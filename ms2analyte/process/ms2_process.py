#!/usr/bin/env python3

"""Process ms2 data for individual file and create ms2 spectra for each ms1 analyte"""

from scipy.stats import linregress
from scipy import interpolate
import os
import pickle


import ms2analyte.config as config
from ms2analyte.process import analyte_create
from ms2analyte.calculate import mass_features


def calculate_slope(row_data, peak_df):
    """Calculate slope between MS1 and MS2 data for each MS2 mass peak by plotting intensity vs intensity
    for scan by scan data"""
    if row_data.count() >= config.matched_scan_minimum:
        aligned_intensities = peak_df.merge(row_data.dropna(), on="scan")
        if len(aligned_intensities.index) >= config.matched_scan_minimum:
            slope_value = linregress(aligned_intensities["intensity"], aligned_intensities[row_data.name])

            # Intensities must vary colinearly if masses are from the same analyte, so r2 must be high.
            # Relationship must also be positive, so r must be > 0

            if slope_value.rvalue ** 2 >= config.slope_r2_cuttoff and slope_value.slope > 0:
                return True
            else:
                return False
        else:
            return False
    else:
        return False


def ms2_to_analyte_vectorized(ms2_data, interpolated_ms1_data, analyte_id):
    """Identify MS2 peaks that align with each analyte (from interpolated intensities of largest MS1 peak)"""
    unassigned_peaks = ms2_data.loc[(ms2_data["analyte_id"].isnull()) &
                                    (ms2_data["scan"].between(interpolated_ms1_data["scan"].min(),
                                                              interpolated_ms1_data["scan"].max())), ["peak_id",
                                                                                                      "intensity",
                                                                                                      "scan"]]
    if not unassigned_peaks.empty:
        intensity_grid = unassigned_peaks.pivot(index='peak_id', columns='scan', values='intensity')
        # Now calculate slope, and assign analyte_id to all peaks in input_df with appropriate slope value.
        slope_match = intensity_grid.apply(calculate_slope, args=(interpolated_ms1_data, ), axis=1)
        matched_peaks = slope_match[slope_match].index.values
        ms2_data.loc[ms2_data["peak_id"].isin(matched_peaks), "analyte_id"] = int(analyte_id)

    print("Completed DIA assignment for analyte " + str(analyte_id))

    return ms2_data


def annotate_ms1_peaks(ms1_data, ms2_data, analyte_list):
    """Interpolate MS1 intensities for the time points for the MS2 scans for the largest mass peak in each analyte.
    Use relative changes in intensity between interpolated MS1 data and real MS2 data to find MS2 peaks that go with
    each analyte. """
    ms2_data["analyte_id"] = None
    # Extract list of unique scan numbers and corresponding retention times
    ms2_scans = ms2_data[["scan", "rt"]].drop_duplicates().sort_values(by=["scan"])
    for analyte in analyte_list:
        max_peak_data = ms1_data[ms1_data["peak_id"] == analyte.max_peak_id][["scan", "rt", "intensity"]].sort_values(by=["scan"])
        interpolated_range = ms2_scans[ms2_scans["scan"].between(max_peak_data["scan"].min(), max_peak_data["scan"].max())].copy()
        if len(interpolated_range.index) >= config.matched_scan_minimum:
            if len(max_peak_data.index) > 3:
                tck = interpolate.splrep(max_peak_data["rt"].to_numpy(), max_peak_data["intensity"].to_numpy(), s=0)
            elif len(max_peak_data.index) == 3:
                tck = interpolate.splrep(max_peak_data["rt"].to_numpy(), max_peak_data["intensity"].to_numpy(), s=0, k=2)
            else:
                continue
            interpolated_intensities = interpolate.splev(interpolated_range["rt"].to_numpy(), tck, der=0)
            interpolated_range["intensity"] = interpolated_intensities
            ms2_data = ms2_to_analyte_vectorized(ms2_data,
                                                 interpolated_range[["scan", "intensity"]],
                                                 analyte.analyte_id)
        else:
            continue

    return ms2_data


def append_ms2_spectra(input_structure, input_type, input_file, ms2_data, ms2_peak_list, analyte_list):
    """Append MS2 spectra to existing MS1 analyte objects"""

    ms2_peak_dict = analyte_create.peak_list_to_dict(ms2_peak_list)

    for ms2_analyte_id, ms2_data in ms2_data.groupby("analyte_id"):
        ms2_peak_data = []
        for peak_id in ms2_data["peak_id"].unique().tolist():
            ms2_peak_data.append(ms2_peak_dict[peak_id])

        for analyte in analyte_list:
            if analyte.analyte_id == ms2_analyte_id:
                analyte.ms2_peak_list = ms2_peak_data

    with open(os.path.join(input_structure.output_directory, input_type,
                           input_file[:-(len(input_structure.ms_data_file_suffix) + 1)] + "_analytes.pickle"), "wb") \
            as g:
        pickle.dump(analyte_list, g)


def append_dda_data(analyte_list, ms1_data, ms2_data, dda_index):
    """Append DDa data to each MassPeak in each Analyte for which DDA data exists"""

    ms1_peak_match_list = []
    ms2_dda_ms_data = {}

    for dda_entry in dda_index:
        ms1_peak_data = ms1_data[(ms1_data["scan"] == dda_entry.current_ms1_scan) &
                                 (ms1_data["mz"].between(dda_entry.selected_ion - 0.5,
                                                        dda_entry.selected_ion + 0.5))]
        if len(ms1_peak_data.index) > 0:
            peak_id_list = ms1_peak_data.peak_id.unique().tolist()
            # Assign matching MS1 peak ID. NOTE: This assumes there is only one matching peak. Because of floating
            # point errors it is possible that there could be more than one match on very rare occassions.
            # Currently this script arbitrarily assigns the match to the first peak.
            ms1_peak_match = peak_id_list[0]
            if len(peak_id_list) > 1:
                print("WARNING: More than one MS1 peak match found for DDA scan " + str(dda_entry.ms2_scan)
                      + ". Setting DDA match to peak_id " + str(ms1_peak_match))
            dda_data = ms2_data[ms2_data["scan"] == dda_entry.ms2_scan]
            if len(dda_data.index) > 0:
                ms1_peak_match_list.append(ms1_peak_match)
                ms2_dda_ms_data[ms1_peak_match] = [dda_data.mz.tolist(), dda_data.intensity.tolist()]
            else:
                # print("ERROR: No DDA data for MS2 scan " + str(dda_entry.ms2_scan))
                pass
        else:
            # print("WARNING: No matching MS1 peaks found for scan " +
            #       str(dda_entry.ms2_scan) + " selected_ion " + str(dda_entry.selected_ion))
            pass

    for analyte in analyte_list:
        for peak in analyte.peak_list:
            if peak.peak_id in ms1_peak_match_list:
                if not peak.dda_data:
                    peak.dda_data = ms2_dda_ms_data[peak.peak_id]
                else:
                    if max(peak.dda_data[1]) < max(ms2_dda_ms_data[peak.peak_id][1]):
                        peak.dda = ms2_dda_ms_data[peak.peak_id]


