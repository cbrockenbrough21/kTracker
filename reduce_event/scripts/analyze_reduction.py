import uproot
import numpy as np

def analyze_declustering(original_file, noisy_file, reduced_file):
    tree_orig = uproot.open(original_file)["tree"]
    tree_noisy = uproot.open(noisy_file)["tree"]
    tree_reduced = uproot.open(reduced_file)["tree"]

    n_events = min(len(tree_orig), len(tree_noisy), len(tree_reduced))

    total_real_hits = 0
    total_noise_hits = 0
    total_final_hits = 0
    noise_removed = 0
    real_lost = 0

    for i in range(n_events):
        orig_det = tree_orig["detectorID"].array(entry_start=i, entry_stop=i+1)[0]
        orig_elem = tree_orig["elementID"].array(entry_start=i, entry_stop=i+1)[0]

        noisy_det = tree_noisy["detectorID"].array(entry_start=i, entry_stop=i+1)[0]
        noisy_elem = tree_noisy["elementID"].array(entry_start=i, entry_stop=i+1)[0]

        red_det = tree_reduced["detectorID"].array(entry_start=i, entry_stop=i+1)[0]
        red_elem = tree_reduced["elementID"].array(entry_start=i, entry_stop=i+1)[0]

        # Sets of hits
        real_hits = set(zip(orig_det, orig_elem))
        all_noisy_hits = set(zip(noisy_det, noisy_elem))
        final_hits = set(zip(red_det, red_elem))

        noise_hits = all_noisy_hits - real_hits

        total_real_hits += len(real_hits)
        total_noise_hits += len(noise_hits)
        total_final_hits += len(final_hits)

        noise_removed += len(noise_hits - final_hits)
        real_lost += len(real_hits - final_hits)

    print(f"Events analyzed: {n_events}")
    print(f"Total real hits in original: {total_real_hits}")
    print(f"Total noise hits injected: {total_noise_hits}")
    print(f"Total hits in reduced file: {total_final_hits}")
    print(f"Noise hits removed: {noise_removed}")
    print(f"Real hits lost: {real_lost}")
    
    if total_noise_hits > 0:
        print(f"Noise removal efficiency: {noise_removed / total_noise_hits:.2%}")
    else:
        print("Noise removal efficiency: N/A (no noise hits)")
        
    if total_real_hits > 0:
        print(f"Real hit preservation rate: {(total_real_hits - real_lost) / total_real_hits:.2%}")
    else:
        print("Real hit preservation rate: N/A")




if __name__ == "__main__":
    original_file = "./../noisy_data_gen/MC_negMuon_Dump_Feb21.root"
    noisy_file = "./../noisy_data_gen/noisy_output.root"
    reduced_file = "./cleaned_output_no_electronic_noise.root"
    analyze_declustering(original_file, noisy_file, reduced_file)
  
  
# import uproot
# import numpy as np
# from collections import Counter

# def extract_hits(tree, i):
#     det = np.array(tree["detectorID"].array(entry_start=i, entry_stop=i+1)[0]).astype(np.int32)
#     elem = np.array(tree["elementID"].array(entry_start=i, entry_stop=i+1)[0]).astype(np.int32)
#     return Counter(zip(det, elem))

# def report_duplicates(name, hits, event_index):
#     duplicates = [k for k, v in hits.items() if v > 1]
#     if duplicates:
#         print(f"\n🔁 Duplicates in {name} file, event {event_index}:")
#         for hit in duplicates:
#             print(f"  Hit {hit} appears {hits[hit]} times")

# def analyze_declustering(original_file, noisy_file, reduced_file):
#     tree_orig = uproot.open(original_file)["tree"]
#     tree_noisy = uproot.open(noisy_file)["tree"]
#     tree_reduced = uproot.open(reduced_file)["tree"]

#     n_events = min(len(tree_orig), len(tree_noisy), len(tree_reduced))

#     total_real_hits = 0
#     total_noise_hits = 0
#     total_final_hits = 0
#     noise_removed = 0
#     real_lost = 0

#     for i in range(n_events):
#         real_hits = extract_hits(tree_orig, i)
#         noisy_hits = extract_hits(tree_noisy, i)
#         reduced_hits = extract_hits(tree_reduced, i)
        
#         if i == 1:
#             print("\n🔍 Manual inspection of event 1:")

#             orig_det = np.array(tree_orig["detectorID"].array(entry_start=i, entry_stop=i+1)[0], dtype=np.int32)
#             orig_elem = np.array(tree_orig["elementID"].array(entry_start=i, entry_stop=i+1)[0], dtype=np.int32)
#             real_set = set(zip(orig_det, orig_elem))
#             print("\nOriginal file hits:")
#             for hit in sorted(real_set):
#                 print(f"  {hit}")

#             noisy_det = np.array(tree_noisy["detectorID"].array(entry_start=i, entry_stop=i+1)[0]).astype(np.int32)
#             noisy_elem = np.array(tree_noisy["elementID"].array(entry_start=i, entry_stop=i+1)[0]).astype(np.int32)
#             noisy_set = set(zip(noisy_det, noisy_elem))
#             print("\nNoisy file hits:")
#             for hit in sorted(noisy_set):
#                 print(f"  {hit}")

#             reduced_det = np.array(tree_reduced["detectorID"].array(entry_start=i, entry_stop=i+1)[0]).astype(np.int32)
#             reduced_elem = np.array(tree_reduced["elementID"].array(entry_start=i, entry_stop=i+1)[0]).astype(np.int32)
#             reduced_set = set(zip(reduced_det, reduced_elem))
#             print("\nReduced file hits:")
#             for hit in sorted(reduced_set):
#                 print(f"  {hit}")

#             # Extra line to see what was flagged as unexpected
#             unexpected_hits = reduced_set - (real_set | (noisy_set - real_set))
#             print("\n⚠️ Unexpected hits in reduced file (not in original or noisy):")
#             for hit in sorted(unexpected_hits):
#                 print(f"  {hit}")


#         # Duplicate checks
#         report_duplicates("original", real_hits, i)
#         report_duplicates("noisy", noisy_hits, i)
#         report_duplicates("reduced", reduced_hits, i)

#         injected_noise = noisy_hits - real_hits
#         removed_noise = injected_noise - reduced_hits
#         lost_real = real_hits - reduced_hits

#         # Diagnose: any hit in reduced that's not from real OR noise?
#         all_possible = real_hits + injected_noise
#         unexpected_hits = reduced_hits - all_possible
#         if sum(unexpected_hits.values()) > 0:
#             print(f"\n⚠️ Unexpected hits in event {i}:")
#             for hit, count in unexpected_hits.items():
#                 print(f"  Hit {hit} appears {count} extra time(s) in reduced file")

#         # Count totals
#         total_real_hits += sum(real_hits.values())
#         total_noise_hits += sum(injected_noise.values())
#         total_final_hits += sum(reduced_hits.values())
#         noise_removed += sum(removed_noise.values())
#         real_lost += sum(lost_real.values())

#         # Optional: diagnostics for early events
#         if i < 5:
#             print(f"\nEvent {i}:")
#             print(f"  Real hits: {sum(real_hits.values())}")
#             print(f"  Noise hits: {sum(injected_noise.values())}")
#             print(f"  Final hits: {sum(reduced_hits.values())}")
#             print(f"  Noise removed: {sum(removed_noise.values())}")
#             print(f"  Real lost: {sum(lost_real.values())}")

#     print(f"\nEvents analyzed: {n_events}")
#     print(f"Total real hits in original: {total_real_hits}")
#     print(f"Total noise hits injected: {total_noise_hits}")
#     print(f"Total hits in reduced file: {total_final_hits}")
#     print(f"Noise hits removed: {noise_removed}")
#     print(f"Real hits lost: {real_lost}")
    
#     if total_noise_hits > 0:
#         print(f"Noise removal efficiency: {noise_removed / total_noise_hits:.2%}")
#     else:
#         print("Noise removal efficiency: N/A (no noise hits)")
        
#     if total_real_hits > 0:
#         print(f"Real hit preservation rate: {(total_real_hits - real_lost) / total_real_hits:.2%}")
#     else:
#         print("Real hit preservation rate: N/A")



# if __name__ == "__main__":
#     original_file = "./../noisy_data_gen/MC_negMuon_Dump_Feb21.root"
#     noisy_file = "./../noisy_data_gen/noisy_output.root"
#     reduced_file = "./cleaned_output_no_electronic_noise.root"
#     analyze_declustering(original_file, noisy_file, reduced_file)
