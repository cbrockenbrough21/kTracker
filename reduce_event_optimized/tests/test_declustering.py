import pytest
from reduce_event.filters.decluster_hits import decluster_hits
import numpy as np

def test_single_hit_cluster():
    detectorIDs = np.array([10], dtype=np.int32)
    elementIDs = np.array([5], dtype=np.int32)
    driftDistances = np.array([0.5], dtype=np.float64)
    tdcTimes = np.array([1000], dtype=np.float64)
    keep_idx = np.array([0], dtype=np.int32)

    result = decluster_hits(detectorIDs, elementIDs, driftDistances, tdcTimes, keep_idx)
    assert list(result) == [0]

def test_two_hit_cluster_keep_one_based_on_drift():
    detectorIDs = np.array([0, 0], dtype=np.int32)
    elementIDs = np.array([5, 6], dtype=np.int32)
    driftDistances = np.array([0.4, 0.2], dtype=np.float64)
    tdcTimes = np.array([1000, 1010], dtype=np.float64)
    keep_idx = np.array([0, 1], dtype=np.int32)

    result = decluster_hits(detectorIDs, elementIDs, driftDistances, tdcTimes, keep_idx)
    assert list(result) == [1]

def test_two_hit_cluster_remove_due_to_noise_D3p():
    detectorIDs = np.array([20, 20], dtype=np.int32)  # D3p group
    elementIDs = np.array([10, 11], dtype=np.int32)
    driftDistances = np.array([0.5, 0.5], dtype=np.float64)
    tdcTimes = np.array([1000, 1005], dtype=np.float64)  # dt < 8
    keep_idx = np.array([0, 1], dtype=np.int32)

    result = decluster_hits(detectorIDs, elementIDs, driftDistances, tdcTimes, keep_idx)
    assert list(result) == []

def test_multi_hit_cluster_noise_remove_all():
    detectorIDs = np.array([10, 10, 10], dtype=np.int32)
    elementIDs = np.array([5, 6, 7], dtype=np.int32)
    driftDistances = np.array([0.5, 0.5, 0.5], dtype=np.float64)
    tdcTimes = np.array([1000, 1005, 1010], dtype=np.float64)  # dt_mean = 5 < 10
    keep_idx = np.array([0, 1, 2], dtype=np.int32)

    result = decluster_hits(detectorIDs, elementIDs, driftDistances, tdcTimes, keep_idx)
    assert list(result) == []

def test_multi_hit_cluster_keep_edges():
    detectorIDs = np.array([10, 10, 10], dtype=np.int32)
    elementIDs = np.array([5, 6, 7], dtype=np.int32)
    driftDistances = np.array([0.5, 0.5, 0.5], dtype=np.float64)
    tdcTimes = np.array([1000, 1200, 1400], dtype=np.float64)  # dt_mean = 200 > 10
    keep_idx = np.array([0, 1, 2], dtype=np.int32)

    result = decluster_hits(detectorIDs, elementIDs, driftDistances, tdcTimes, keep_idx)
    assert list(result) == [0, 2]
