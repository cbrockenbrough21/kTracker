import numpy as np
from numba import njit, types
from numba.typed import List


def decluster_hits(detectorIDs, elementIDs, driftDistances, tdcTimes, keep_idx):
   """
   Python wrapper: sorts keep_idx and passes to JIT declustering.
   """
   sort_order = np.lexsort((elementIDs[keep_idx], detectorIDs[keep_idx]))
   sorted_keep_idx = keep_idx[sort_order].astype(np.int32)

   return decluster_hits_sorted(
       detectorIDs.astype(np.int32),
       elementIDs.astype(np.int32),
       driftDistances.astype(np.float64),
       tdcTimes.astype(np.float64),
       sorted_keep_idx,
   )


@njit
def decluster_hits_sorted(detectorIDs, elementIDs, driftDistances, tdcTimes, sorted_hits):
   result = List.empty_list(types.int32)
   cluster = List.empty_list(types.int32)


   for idx in sorted_hits:
       if len(cluster) == 0:
           cluster.append(idx)
           continue


       last = cluster[-1]
       if detectorIDs[idx] == detectorIDs[last] and abs(elementIDs[idx] - elementIDs[last]) <= 1:
           cluster.append(idx)
       else:
           retained = process_cluster(cluster, detectorIDs, driftDistances, tdcTimes)
           for i in retained:
               result.append(i)
           cluster = List.empty_list(types.int32)
           cluster.append(idx)


   if len(cluster) > 0:
       retained = process_cluster(cluster, detectorIDs, driftDistances, tdcTimes)
       for i in retained:
           result.append(i)


   # Manual conversion to np.array (Numba-compatible)
   final_result = np.empty(len(result), dtype=np.int32)
   for i in range(len(result)):
       final_result[i] = result[i]
   return final_result


@njit
def process_cluster(cluster, detectorIDs, driftDistances, tdcTimes):
   result = List.empty_list(types.int32)
   n = len(cluster)


   if n == 1:
       result.append(cluster[0])
       return result


   det_id = detectorIDs[cluster[0]]
   for i in cluster:
       if detectorIDs[i] != det_id:
           return cluster  # already a typed.List[int32]


   hw = 0.635 / 2
   group_id = det_id // 5
   if group_id == 1:
       hw = 2.083 / 2
   elif group_id in (2, 3):
       hw = 2.021 / 2
   elif group_id in (4, 5):
       hw = 2.000 / 2


   if n == 2:
       i0, i1 = cluster
       drift0, drift1 = driftDistances[i0], driftDistances[i1]
       tdc0, tdc1 = tdcTimes[i0], tdcTimes[i1]


       if 19 <= det_id <= 24 and abs(tdc0 - tdc1) < 8:
           return result  # return empty


       w_max = 0.9 * hw
       w_min = 0.4 * hw


       if drift0 > w_max and drift1 > w_min:
           result.append(i1)
       elif drift1 > w_max and drift0 > w_min:
           result.append(i0)
       else:
           result.append(i0)
           result.append(i1)
       return result


   # n >= 3
   tdcs = [tdcTimes[i] for i in cluster]
   tdcs.sort()


   dt_mean = 0.0
   for i in range(1, len(tdcs)):
       dt_mean += abs(tdcs[i] - tdcs[i - 1])
   dt_mean /= max(len(tdcs) - 1, 1)


   if dt_mean < 10:
       return result
   else:
       result.append(cluster[0])
       result.append(cluster[-1])
       return result

