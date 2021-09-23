#!/usr/bin/env python3

"""
Set of scripts for calculating chromatographic peak properties of peaks and analytes from mass spectrometry data

"""

import pandas as pd
import ms2analyte.config as config
from scipy.stats import linregress


def find_rt(peak_data):
    """Determine retention time for a peak, given Pandas dataframe of rt and intensity"""
    max_intensity = peak_data["intensity"].max()
    rt = round(peak_data[peak_data["intensity"] == max_intensity].iloc[0]["rt"], 4)

    return rt


def rt_match(rt1, rt2):
    """Determine whether two analytes have the same rt, within defined rt error"""
    if rt1 - config.rt_error <= rt2 <= rt1 + config.rt_error:
        return True
    else:
        return False


def average_rt(rt_list):
    """Calculate average rt for a set of peaks (e.g. experiment analyte rt based on replicate analyte rts).
    Currently this is a simple average calculation, but is broken out as a separate function for possible future
    changes to rt determination strategy

    """

    return round(sum(rt_list)/len(rt_list), 4)


def rt_peak_shape_match(sample_df, compare_df):
    """Test peak shape between analytes from two different samples to identify peak shape match
    Performs all-on-all comparison for peaks from each replicate from each analyte. Stops if any pair gives a match.
    Used in basketing as part of the process to define experiment analytes from different samples.

    """
    # Ensure there is the minimum scan overlap required to compare peaks
    if sample_df["scan"].min() <= compare_df["scan"].max() - config.matched_scan_minimum - 1 \
            and sample_df["scan"].max() >= compare_df["scan"].min() - config.matched_scan_minimum - 1:
        # Iterate through replicates, comparing max intensity peak from each replicate against the corresponding
        # max peaks from each compare replicate to look for peak shape overlap
        for replicate, sample_data in sample_df.groupby("replicate"):
            max_intensity_peak_id = sample_data.loc[sample_data["intensity"].idxmax(), "peak_id"]
            sample_peak = sample_data[sample_data["peak_id"] == max_intensity_peak_id]
            for compare_replicate, compare_data in compare_df.groupby("replicate"):
                max_intensity_compare_peak_id = compare_data.loc[compare_data["intensity"].idxmax(), "peak_id"]
                compare_peak = compare_data[compare_data["peak_id"] == max_intensity_compare_peak_id]
                merged_data = pd.merge(sample_peak, compare_peak, how='inner', on='scan')
                if len(merged_data.index) > config.matched_scan_minimum:
                    slope_value = linregress(merged_data["intensity_x"], merged_data["intensity_y"])
                    if slope_value.rvalue ** 2 >= config.experiment_analyte_slope_r2_cuttoff and slope_value.slope > 0:
                        return True
