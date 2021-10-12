#!/usr/bin/env python3

"""Process raw input data for individual files and create analytes from individual runs"""

import os
import sys
import pickle

from ms2analyte.converters import waters, thermo, agilent
from ms2analyte.process import filters, peak_create, analyte_create, ms2_process
from ms2analyte.output import tableau
from ms2analyte import config


def peak_process(input_data, input_file, input_structure, input_type, ms_type):
    """Import raw initial dataframe and process to assign unique peak IDs to all mass features.
    Same script appropriate for MS1 or MS2 data, for samples or blanks, with or without IMS"""
    input_data = peak_create.mass_bin(input_data)
    tableau.full_export(input_file, input_data, input_structure, input_type, subname="massbin", mstype=ms_type)
    input_data = peak_create.remove_small_peaks(input_data)
    input_data = peak_create.remove_no_maxima_peaks(input_data)
    input_data = peak_create.mass_range_check(input_data)
    # tableau.full_export(input_file, input_data, input_structure, input_type, subname="massrange", mstype=ms_type)
    if not input_structure.ims_exists:
        if peak_create.find_scan_duplicates(input_data):
            input_data = peak_create.resolve_scan_duplicates(input_data)
    input_data = peak_create.isotope_filter(input_data)
    tableau.full_export(input_file, input_data, input_structure, input_type, subname="isotopefilter", mstype=ms_type)
    input_data = peak_create.scan_split(input_data)
    # tableau.full_export(input_file, input_data, input_structure, input_type, subname="scansplit", mstype=ms_type)
    if input_structure.ims_exists:
        input_data = peak_create.drift_split(input_data)
        # tableau.full_export(input_file, input_data, input_structure, input_type, subname="driftsplit", mstype=ms_type)
        input_data = peak_create.scan_split(input_data)
    input_data = peak_create.remove_no_maxima_peaks(input_data)
    input_data = peak_create.peak_minima_trim(input_data)
    input_data = peak_create.peak_split_max(input_data)
    # tableau.full_export(input_file, input_data, input_structure, input_type, subname="peaksplit", mstype=ms_type)
    input_data = peak_create.strip_flat_peaks(input_data)
    tableau.full_export(input_file, input_data, input_structure, input_type, subname="flatpeakstrip", mstype=ms_type)
    input_data = peak_create.scan_split(input_data)
    input_data = peak_create.remove_small_peaks(input_data)
    input_data = peak_create.remove_no_maxima_peaks(input_data)

    print("Completed initial processing of " + ms_type + " data for sample " + input_file)

    return input_data


def sample_process(input_file, input_structure, input_type):
    """Create analytes for a single data file
    Either sample or blank OK in this function

    """
    sample_name = input_file[:-(len(input_structure.ms_data_file_suffix) + 1)]

    print("Starting processing for sample " + sample_name)

    ms1_data = None
    ms2_data = None
    dda_index = None

    if input_structure.instrument_manufacturer == "Waters":
        if input_structure.ms_data_file_type == "func001":
            if input_type == "Samples":
                ms1_data = filters.intensity(waters.func(os.path.join(input_structure.sample_directory, input_file)))
                if input_structure.ms2_exists:
                    ms2_data = filters.intensity(waters.func(os.path.join(input_structure.sample_directory, "func002",
                                                                          input_file)))
            elif input_type == "Blanks":
                ms1_data = filters.intensity(waters.func(os.path.join(input_structure.blank_directory, input_file)))
                if input_structure.ms2_exists:
                    ms2_data = filters.intensity(waters.func(os.path.join(input_structure.blank_directory, "func002",
                                                                          input_file)))
            else:
                print("Input type not recognized. Options are 'Samples' or 'Blanks'")
                sys.exit()
        elif input_structure.ms_data_file_type == "mzML":
            if input_type == "Samples":
                ms1_data, ms2_data, dda_index = waters.mzml(os.path.join(input_structure.sample_directory, input_file),
                                                            input_structure)
            elif input_type == "Blanks":
                ms1_data, ms2_data, dda_index = waters.mzml(os.path.join(input_structure.blank_directory, input_file),
                                                            input_structure)
            else:
                print("Input type not recognized. Options are 'Samples' or 'Blanks'")
                sys.exit()
        else:
            print("File type not recognized. Options are: 'func001', 'mzML'")
            sys.exit()

    elif input_structure.instrument_manufacturer == "Thermo":
        if input_structure.ms_data_file_type == "mzML":
            if input_type == "Samples":
                ms1_data = filters.intensity(thermo.mzml(os.path.join(input_structure.sample_directory, input_file)))
            elif input_type == "Blanks":
                ms1_data = filters.intensity(thermo.mzml(os.path.join(input_structure.blank_directory, input_file)))
            else:
                print("Input type not recognized. Options are 'Samples' or 'Blanks'")
                sys.exit()
        elif input_structure.ms_data_file_type == "Proteowizard mzML":
            if input_type == "Samples":
                ms1_data = filters.intensity(thermo.proteowizard_mzml(os.path.join(input_structure.sample_directory,
                                                                                     input_file)))
            elif input_type == "Blanks":
                ms1_data = filters.intensity(thermo.proteowizard_mzml(os.path.join(input_structure.blank_directory,
                                                                                     input_file)))
            else:
                print("Input type not recognized. Options are 'Samples' or 'Blanks'")
                sys.exit()

        else:
            print("File type not recognized. Options are: 'mzML', 'Proteowizard mzML")
            sys.exit()

    elif input_structure.instrument_manufacturer == "Agilent":
        if input_structure.ms_data_file_type == "mzXML":
            if input_type == "Samples":
                ms1_data = filters.intensity(agilent.mzxml_import(os.path.join(input_structure.sample_directory,
                                                                               input_file)))
            elif input_type == "Blanks":
                ms1_data = filters.intensity(agilent.mzxml_import(os.path.join(input_structure.blank_directory,
                                                                               input_file)))
            else:
                print("Input type not recognized. Options are 'Samples' or 'Blanks'")
                sys.exit()
        else:
            print("File type not recognized. Options are: 'mzXML'")
            sys.exit()
    else:
        print("Manufacturer not recognized. Options are: 'Waters', 'Thermo', 'Agilent'")
        sys.exit()

    # Process ms1 data
    print("Started MS1 data analysis for " + sample_name)
    ms1_data = peak_process(ms1_data, input_file, input_structure, input_type, "ms1")
    peak_list = analyte_create.peak_df_to_obj(ms1_data)
    ms1_data = analyte_create.peak_to_analyte_vectorized(peak_list, ms1_data)
    analyte_list = analyte_create.create_analytes_objects(peak_list, ms1_data)
    # analyte_list = analyte_create.analyte_isotope_filter(analyte_list)
    tableau.full_export(input_file, ms1_data, input_structure, input_type, mstype="ms1")

    with open(os.path.join(input_structure.output_directory, input_type,
                           input_file[:-(len(input_structure.ms_data_file_suffix) + 1)] + "_ms1_analytes.pickle"),
              "wb") as g:
        pickle.dump(analyte_list, g)

    with open(os.path.join(input_structure.output_directory, input_type,
                           input_file[:-(len(input_structure.ms_data_file_suffix) + 1)] + "_ms1_dataframe.pickle"),
              "wb") as h:
        pickle.dump(ms1_data, h)

    # process ms2 data
    if input_structure.ms2_exists:
        print("Started MS2 data analysis for " + sample_name)
        if input_structure.ms2_type == "DIA":
            ms2_data = peak_process(ms2_data, input_file, input_structure, input_type, "ms2")
            ms2_data= ms2_process.annotate_ms1_peaks(ms1_data, ms2_data, analyte_list)
            ms2_peak_list = analyte_create.peak_df_to_obj(ms2_data)
            ms2_process.append_ms2_spectra(input_structure, input_type, input_file, ms2_data, ms2_peak_list, analyte_list)
            tableau.full_export(input_file, ms2_data, input_structure, input_type, mstype="ms2")

            with open(os.path.join(input_structure.output_directory, input_type,
                                   input_file[:-(len(input_structure.ms_data_file_suffix) + 1)] + "_ms2_dataframe.pickle"),
                      "wb") as i:
                pickle.dump(ms2_data, i)

            tableau.ms1_ms2_combined_export(input_structure, input_type, sample_name)

        elif input_structure.ms2_type == "DDA":
            with open(os.path.join(input_structure.output_directory,
                                   input_file[:-(len(input_structure.ms_data_file_suffix) + 1)] + "_dda_index.pickle"),
                      "wb") as g:
                pickle.dump(dda_index, g)

            ms2_process.append_dda_data(analyte_list, ms1_data, ms2_data, dda_index)
            # Update saved analyte list
            with open(os.path.join(input_structure.output_directory, input_type,
                                   input_file[:-(len(input_structure.ms_data_file_suffix) + 1)] + "_ms1_analytes.pickle"),
                      "wb") as g:
                pickle.dump(analyte_list, g)
