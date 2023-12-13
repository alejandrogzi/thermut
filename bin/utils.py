#!/usr/bin/env python3


import numpy as np
from scipy.interpolate import interp1d


MIN_DEPTH = 5


def mask_timepoints(
    times,
    alts,
    depths,
    var_type,
    cutoff_idx,
    depth_fold_change,
    depth_change_pvalue,
    min_depth=MIN_DEPTH,
):
    # first make a copy of alts and depths
    # so that we can modify in place without worrying
    masked_alts = np.copy(alts)
    masked_depths = np.copy(depths)

    # zero out timepoints that don't pass depth threshold
    masked_alts[masked_depths < min_depth] = 0
    masked_depths[masked_depths < min_depth] = 0

    # masked_alts -= masked_alts*(masked_depths < min_depth)
    # masked_depths -= masked_depths*(masked_depths < min_depth)

    # did we infer that a deletion happened?
    # if (var_type=='sv' and depth_fold_change < -2 and depth_change_pvalue < 1e-04) or (var_type!='sv' and depth_fold_change < -1 and depth_change_pvalue < 1e-03):
    if depth_change_pvalue < 1e-02:
        # deletion nearby, trim timecourse
        masked_alts[cutoff_idx:] = 0
        masked_depths[cutoff_idx:] = 0

    good_idxs = np.nonzero(masked_depths > 0.5)[0]

    return good_idxs, masked_alts, masked_depths


###########
#
# Naive frequency estimator, # alts / # depths (0 if no depth)
#
###########
def estimate_frequencies(alts, depths):
    return alts * 1.0 / (depths + (depths == 0))


def create_interpolation_function(times, freqs, tmax=100000, kind="linear"):
    # can create it for anything!

    padded_times = np.zeros(len(times) + 1)
    padded_freqs = np.zeros(len(times) + 1)
    padded_times[0 : len(times)] = times
    padded_freqs[0 : len(times)] = freqs
    padded_times[-1] = tmax
    padded_freqs[-1] = freqs[-1]

    interpolating_function = interp1d(
        padded_times, padded_freqs, kind=kind, bounds_error=True
    )

    return interpolating_function


###########
#
# Naive frequency estimator for clones, restricted to clones with >=min_depth
#
###########
def estimate_clone_frequencies(times, alts, depths, min_depth=20, allowed_times=None):
    masked_times = times[depths > min_depth]
    masked_alts = alts[depths > min_depth]
    masked_depths = depths[depths > min_depth]
    masked_freqs = masked_alts * 1.0 / masked_depths

    return masked_times, masked_freqs


###########
#
# Estimate depth fold change (log2) relative to genome-wide median
#
###########
def estimate_depth_fold_changes(avg_depths, depths, min_depth=20):
    if (depths <= 0).all():
        return np.zeros_like(depths)

    normalization = depths[depths > min_depth][0] / avg_depths[depths > min_depth][0]

    return np.log2(
        (depths + (depths == 0)) / (avg_depths + (avg_depths == 0)) / normalization
    )
