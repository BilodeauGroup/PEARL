import pickle
from MDAnalysis import Universe
from pi_pi_analysis_complete.analysis.universe_loader import load_universes
from pi_pi_analysis_complete.analysis.ring_utils import find_all_rings
from pi_pi_analysis_complete.analysis.pi_pi_detector import find_pi_stacking_rings
from pi_pi_analysis_complete.config import sequences, replicates, base_dir, save_start_frame, save_end_frame, step
import os
os.makedirs("data", exist_ok=True)
def run_analysis():
    universes = load_universes(sequences, replicates, base_dir)
    pi_pi_data = {}

    for label, psf_file, dcd_file in universes:
        print(f"Processing: {label}")
        try:
            u = Universe(psf_file, dcd_file, topology_format="PARM7")
            ag = u.select_atoms("all")

            all_counts_parallel = []
            all_counts_t_shaped = []
            all_times_ns = []

            for i, ts in enumerate(u.trajectory[save_start_frame:save_end_frame:step]):
                current_frame_index = save_start_frame + i * step
                current_time_fs = current_frame_index * u.trajectory.dt
                current_time_ns = current_time_fs / 1000.0
                all_times_ns.append(current_time_ns)

                aromatic_selection = u.select_atoms("resname TYR PHE DPHE DTYR")
                ring_list = find_all_rings(aromatic_selection)

                if len(ring_list) < 2:
                    all_counts_parallel.append(0)
                    all_counts_t_shaped.append(0)
                    continue

                parallel, t_shaped = find_pi_stacking_rings(ring_list, ag)
                all_counts_parallel.append(len(parallel))
                all_counts_t_shaped.append(len(t_shaped))

            pi_pi_data[label] = {
                "times_ns": all_times_ns,
                "parallel_counts": all_counts_parallel,
                "t_shaped_counts": all_counts_t_shaped,
                "step": step,
                "start_frame": save_start_frame,
                "end_frame": save_end_frame,
                "dt_fs": u.trajectory.dt
            }

        except Exception as e:
            print(f"Error processing {label}: {e}")

    with open("data/pi_pi_interactions_data.pkl", "wb") as f:
        pickle.dump(pi_pi_data, f)

    print("✅ Saved π–π interaction data to 'data/pi_pi_interactions_data.pkl'")

if __name__ == "__main__":
    run_analysis()
