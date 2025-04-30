import pytest
from reduce_event.filters.decluster_hits import decluster_hits
import numpy as np

def test_single_hit_cluster():
    detectorIDs = [10]
    elementIDs = [5]
    driftDistances = [0.5]
    tdcTimes = [1000]
    keep_idx = [0]

    result = decluster_hits(detectorIDs, elementIDs, driftDistances, tdcTimes, keep_idx)
    assert result == [0]

def test_two_hit_cluster_keep_one_based_on_drift():
    detectorIDs = [0, 0]  # D0 detector
    elementIDs = [5, 6]
    driftDistances = [0.4, 0.2]  # good drift distances
    tdcTimes = [1000, 1010]
    keep_idx = [0, 1]

    result = decluster_hits(detectorIDs, elementIDs, driftDistances, tdcTimes, keep_idx)
    # Should keep hit 1 because drift[0] > w_max and drift[1] > w_min
    assert result == [1]

def test_two_hit_cluster_remove_due_to_noise_D3p():
    detectorIDs = [20, 20]  # D3p range
    elementIDs = [10, 11]
    driftDistances = [0.5, 0.5]
    tdcTimes = [1000, 1005]  # close TDC difference (<8)
    keep_idx = [0, 1]

    result = decluster_hits(detectorIDs, elementIDs, driftDistances, tdcTimes, keep_idx)
    # Should remove both
    assert result == []

def test_multi_hit_cluster_noise_remove_all():
    detectorIDs = [10, 10, 10]
    elementIDs = [5, 6, 7]
    driftDistances = [0.5, 0.5, 0.5]
    tdcTimes = [1000, 1005, 1010]  # dt_mean = 5 < 10
    keep_idx = [0, 1, 2]

    result = decluster_hits(detectorIDs, elementIDs, driftDistances, tdcTimes, keep_idx)
    # Should remove all
    assert result == []

def test_multi_hit_cluster_keep_edges():
    detectorIDs = [10, 10, 10]
    elementIDs = [5, 6, 7]
    driftDistances = [0.5, 0.5, 0.5]
    tdcTimes = [1000, 1200, 1400]  # dt_mean = 200 > 10
    keep_idx = [0, 1, 2]

    result = decluster_hits(detectorIDs, elementIDs, driftDistances, tdcTimes, keep_idx)
    # Should keep only first and last
    assert result == [0, 2]
