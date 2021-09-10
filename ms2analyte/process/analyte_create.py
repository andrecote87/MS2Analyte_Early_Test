#!/usr/bin/env python3

"""Tools to build analytes from peak output from peak_create.py"""

from scipy.stats import linregress
from operator import attrgetter

from ms2analyte.calculate import mass_features, rt_features
import ms2analyte.config as config


class MassPeak:
    """Class containing all relevant data for each individual mass peak from peak_create"""
    def __init__(self, dataframe, peak_id, max_intensity, min_scan, max_scan, scan_array, average_mass, rt,
                 average_drift, peak_assigned):
        self.dataframe = dataframe
        self.peak_id = peak_id
        self.max_intensity = max_intensity
        self.min_scan = min_scan
        self.max_scan = max_scan
        self.scan_array = scan_array
        self.average_mass = average_mass
        self.rt = rt
        self.average_drift = average_drift
        self.peak_assigned = peak_assigned
        self.dda_data = None                # [(mz_array), (intensity_array)]. NOTE: If the same mass is selected for
                                            # DDA multiple times, the scan with the highest DDA intensity is retained


class Analyte:
    """Class containing data on peak ids for each assigned analyte"""
    def __init__(self, analyte_id, peak_list):
        self.analyte_id = analyte_id
        self.peak_list = peak_list              # List of MS1 MassPeak objects
        self.ms2_peak_list = None               # List of MS2 MassPeak objects from DIA experiments
        self.max_peak_id = None                 # peak_id of MS1 peak with max intensity
        self.max_peak_intensity_mass = None     # mz value for MS1 peak with max intensity
        self.max_peak_intensity = None          # intensity of MS1 peak with max intensity
        self.analyte_rt = None
        self.replicate_match = None
        self.replicate_analyte_id = None
        self.experiment_match = None
        self.experiment_analyte_id = None
        self.blank_match = None


def peak_df_to_obj(input_data):
    """Create peak objects from Pandas dataframe"""
    peak_list = []

    for peak, data in input_data.groupby("peak_id"):
        peak_id = peak
        max_intensity = data["intensity"].max()
        min_scan = data["scan"].min()
        max_scan = data["scan"].max()
        scan_array = data["scan"].values
        average_mass = round(data["mz"].mean(), 4)
        rt = rt_features.find_rt(data[["rt", "intensity"]])
        average_drift = round(data["drift"].mean(), 1)
        peak_assigned = False

        peak_list.append(MassPeak(data, peak_id, max_intensity, min_scan, max_scan, scan_array, average_mass, rt,
                                  average_drift, peak_assigned))

    sorted_peak_list = sorted(peak_list, key=attrgetter("max_intensity"), reverse=True)

    return sorted_peak_list


def peak_list_to_dict(peak_list):
    """Create dict of peaks"""
    peak_dict = {}

    for peak in peak_list:
        peak_dict[peak.peak_id] = peak

    return peak_dict


def calculate_slope(row_data, peak_df):

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


def peak_to_analyte_vectorized(peak_list, input_df):

    analyte_id = 1
    completed_peaks = []

    for peak in peak_list:
        peak_id_value = peak.peak_id
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
                analyte_id += 1

    print("Total number of peaks = " + str(len(peak_list)))
    print("Total number of analytes = " + str(input_df["analyte_id"].max()))

    return input_df


def create_analytes_objects(peak_list, input_df):

    peak_dict = peak_list_to_dict(peak_list)

    analyte_list = []
    for analyte_id, analyte_data in input_df.groupby("analyte_id"):
        peak_data = []
        for peak_id in analyte_data["peak_id"].unique().tolist():
            peak_data.append(peak_dict[peak_id])
        analyte = Analyte(analyte_id, peak_data)
        max_peak = analyte_data.loc[analyte_data["intensity"].idxmax()]["peak_id"]
        analyte.max_peak_id = max_peak
        analyte.max_peak_intensity = analyte_data.loc[analyte_data["peak_id"] == max_peak]["intensity"].max()
        analyte.max_peak_intensity_mass = peak_dict[max_peak].average_mass
        analyte.analyte_rt = rt_features.find_rt(peak_dict[max_peak].dataframe)
        analyte_list.append(analyte)

    return analyte_list


def analyte_isotope_filter(analyte_list):
    """Filter analytes, and only keep those that contain at least one peak with a 13C isotope"""
    filtered_analyte_list = []

    for analyte in analyte_list:
        isotope_peak = False
        peak_mass_list = []
        for peak in analyte.peak_list:
            peak_mass_list.append(peak.average_mass)
        sorted_peak_mass_list = sorted(peak_mass_list)
        for index, peak_mass in enumerate(sorted_peak_mass_list[:-1]):
            for peak_mass2 in sorted_peak_mass_list[index + 1:]:
                if mass_features.isotope_match(peak_mass, peak_mass2):
                    isotope_peak = True
                    break
        if isotope_peak:
            filtered_analyte_list.append(analyte)

    print("Total number of analytes with isotopes = " + str(len(filtered_analyte_list)))

    return filtered_analyte_list
