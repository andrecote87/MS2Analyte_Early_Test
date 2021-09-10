#!/usr/bin/env python3

"""Tool to convert Waters source data in to MS2Analyte input format"""

import pandas as pd
import pymzml


import ms2analyte.config as config
from ms2analyte.file_handling import data_import


def func(file_path):
    """Import MS data from func files derived from MSeXpress"""
    input_data = pd.read_csv(file_path)
    reformatted_data = input_data[["ScanIndex", "ScanTimeMin", "Mz", "Drift", "Intensity"]]
    reformatted_data = reformatted_data.rename(index=str, columns={"ScanIndex": "scan", "ScanTimeMin": "rt", "Mz": "mz",
                                                "Drift": "drift", "Intensity": "intensity"})

    return reformatted_data


def mzml(file_path, input_structure):
    """Import mzML files derived from applying MSConvert to .raw files.
    Tested with both DIA and DDA data.

    """

    # Waters data includes the lockspray internal calibrant scans as 'MS1' data. These are differentiated from true
    # MS1 data by the 'function' attribute in the spectrum element. Data MS1 scans are function 1. Lockspray scans are
    # assigned the highest possible function number (floating, depends on how many DDA scans were permitted during
    # acquisition setup). Commonly lockspray function=5. This is always 3 for MSe (DIA) data.
    # In order to filter out the lockspray data it is therefore necessary to filter the ms_level 1 scans to only
    # retain those where the function value is also 1.
    # NOTE: For MS2 data this is not an issue, because all MS2 data have ms level = 2, and are therefore all
    # legitimate for inclusion.

    run = pymzml.run.Reader(str(file_path))
    ms1_input_data = []
    ms2_input_data = []
    dda_index = []
    current_ms1_scan = 0
    for spec in run:
        # Handle MS1 and MS2 scans separately
        if spec.ms_level == 1:
            # Skip lockspray or other functions if there are any
            # If not, this is probably not Waters data and should be fine...
            fn = spec.id_dict.get("function")
            if fn is not None:
                if fn != 1:
                    continue
            ms1_input_data.append(extract_scan_data(spec, input_structure))
            current_ms1_scan = spec.ID
        elif spec.ms_level == 2:
            ms2_input_data.append(extract_scan_data(spec, input_structure))
            if input_structure.ms2_type == "DDA":
                for precursor in spec.selected_precursors:
                    dda_index.append(data_import.DdaIndex(current_ms1_scan, spec.ID, precursor["mz"]))
        else:
            continue
        # Print import progress. Useful because importing large mzML files can be slow.
        if spec.ms_level == 1 and spec.ID % 100 == 0 and spec.ID > 0:
            print("Completed import of scan " + str(spec.ID))

    return pd.concat(ms1_input_data, ignore_index=True), pd.concat(ms2_input_data, ignore_index=True), dda_index


def extract_scan_data(spec, input_structure):
    """Extract data from individual scans from mzML files including drift time"""

    scan_data = pd.DataFrame(spec.mz, columns=["mz"])
    scan_data["intensity"] = spec.i
    scan_data["scan"] = spec.ID
    scan_data["rt"] = round(spec.scan_time_in_minutes(), 2)
    # NOTE: Drift time extraction not complete yet.
    if input_structure.ims_exists:
        scan_data["drift"] = spec.drift * 100
    else:
        scan_data["drift"] = None

    scan_data = scan_data.astype({'scan': 'int64',
                                  'rt': 'float64',
                                  'mz': 'float64',
                                  'drift': "float64",
                                  'intensity': 'int64'})

    # Remove data below intensity cutoff, and reorder columns to standard order
    return scan_data[scan_data["intensity"] >= config.intensity_cutoff][["scan", "rt", "mz", "drift", "intensity"]]
