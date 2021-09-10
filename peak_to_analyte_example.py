import os
import pickle
import numpy as np
from scipy.stats import linregress
import timeit
import functools

import ms2analyte.config as config


class Analyte:
    """Class containing data on peak ids for each assigned analyte"""
    def __init__(self, analyte_id, peak_list, max_peak_id, max_peak_intensity_mass, max_peak_intensity, analyte_rt,
                 replicate_match, replicate_analyte_id, experiment_match, experiment_analyte_id, blank_match):
        self.analyte_id = analyte_id
        self.peak_list = peak_list
        self.max_peak_id = max_peak_id
        self.max_peak_intensity_mass = max_peak_intensity_mass
        self.max_peak_intensity = max_peak_intensity
        self.analyte_rt = analyte_rt
        self.replicate_match = replicate_match
        self.replicate_analyte_id = replicate_analyte_id
        self.experiment_match = experiment_match
        self.experiment_analyte_id = experiment_analyte_id
        self.blank_match = blank_match


def calculate_slope(row_data, peak_df):

    if row_data.count() >= config.matched_scan_minimum:
        aligned_intensities = peak_df.merge(row_data.dropna(), on="scan")
        slope_value = linregress(aligned_intensities["intensity"], aligned_intensities[row_data.name])

    # Intensities must vary colinearly if masses are from the same analyte, so r2 must be high.
    # Relationship must also be positive, so r must be > 0

        if slope_value.rvalue ** 2 >= config.slope_r2_cuttoff and slope_value.slope > 0:
            return True
        else:
            return False
    else:
        return False


def peak_to_analyte_vectorized(input_df):

    analyte_id = 1
    available_peaks = input_df["peak_id"].unique().tolist()
    print(available_peaks)

    completed_peaks = []

    for peak_id_value in available_peaks:
        if peak_id_value not in completed_peaks:
            completed_peaks.append(peak_id_value)
            input_df.loc[input_df["peak_id"] == peak_id_value, "analyte_id"] = int(analyte_id)
            peak_data = input_df.loc[input_df["peak_id"] == peak_id_value, ["intensity", "scan"]]
            unassigned_peaks = input_df.loc[(input_df["analyte_id"].isnull()) &
                                        (input_df["scan"].between(peak_data["scan"].min(), peak_data["scan"].max())),
                                            ["peak_id", "intensity", "scan"]]
            if not unassigned_peaks.empty:
                intensity_grid = unassigned_peaks.pivot(index='peak_id', columns='scan', values='intensity')
                # Now calculate slope, and assign analyte_id to all peaks in input_df with appropriate slope value.
                slope_match = intensity_grid.apply(calculate_slope, args=(peak_data, ), axis=1)
                matched_peaks = slope_match[slope_match].index.values
                completed_peaks = completed_peaks + matched_peaks.tolist()
                input_df.loc[input_df["peak_id"].isin(matched_peaks), "analyte_id"] = int(analyte_id)
            if analyte_id % 50 == 1:
                print("Finished analyte " + str(analyte_id))
            analyte_id += 1


def peak_list_to_dict(peak_list):
    """Create dict of peaks"""
    peak_dict = {}

    for peak in peak_list:
        peak_dict[peak.peak_id] = peak

    return peak_dict


def peak_to_analyte(peak_list, input_data):
    """Create analyte from peaks by comparing slope of intensities on scan by scan basis for peaks with overlapping
    scans

    """
    analyte_list = []
    current_analyte_id = 1

    peak_dict = peak_list_to_dict(peak_list)

    # Look at most intense peak in list. If not already assigned to an analyte, start new analyte and add peak

    for origin_peak in peak_list:
        if not origin_peak.peak_assigned:
            origin_peak.peak_assigned = True
            current_analyte = Analyte(current_analyte_id, [origin_peak], None, None, None, None, None, None, None,
                                      None, None)

            # look iteratively at all peaks and decide if peak shape matches first peak in new analyte. If so, add

            for test_peak in peak_list:
                if not test_peak.peak_assigned:
                    matched_scans = sorted(np.intersect1d(origin_peak.scan_array, test_peak.scan_array))
                    if len(matched_scans) >= config.matched_scan_minimum:
                        origin_peak_intensities = input_data[(input_data["peak_id"] == origin_peak.peak_id) &
                                                             (input_data["scan"].isin(matched_scans))].sort_values(by=["scan"])["intensity"].values
                        test_peak_intensities = input_data[(input_data["peak_id"] == test_peak.peak_id) &
                                                           (input_data["scan"].isin(matched_scans))].sort_values(by=["scan"])["intensity"].values
                        slope_value = linregress(origin_peak_intensities, test_peak_intensities)

                        # Intensities must vary colinearly if masses are from the same analyte, so r2 must be high.
                        # Relationship must also be positive, so r must be > 0

                        if slope_value.rvalue ** 2 >= config.slope_r2_cuttoff and slope_value.slope > 0:
                            test_peak.peak_assigned = True
                            current_analyte.peak_list.append(test_peak)

            if len(current_analyte.peak_list) >= config.analyte_peak_minimum:
                analyte_list.append(current_analyte)
                # print("Completed analysis for analyte " + str(current_analyte_id))
                if current_analyte_id % 50 == 1:
                    print("Finished analyte " + str(current_analyte_id))
                current_analyte_id += 1


if __name__ == "__main__":

    with open(os.path.join("data", "raw_input", "FERM_BP_3421_wt_ACE_12C_R1_dataframe.pickle"), "rb") as f:
        input_df = pickle.load(f)

    with open(os.path.join("data", "raw_input", "temp_peak_list.pickle"), "rb") as g:
        peak_list = pickle.load(g)

    simplified_input_df = input_df.drop(["analyte_id", "replicate_analyte_id"], axis=1)

    t = timeit.Timer(functools.partial(peak_to_analyte_vectorized, simplified_input_df))
    print("Vectorized time: " + str(t.timeit(1)))
    u = timeit.Timer(functools.partial(peak_to_analyte, peak_list, simplified_input_df))
    print("Nested for loops time: " + str(u.timeit(1)))
