from Bio import Phylo
import h5py
import pickle
import matplotlib.pyplot as plt
import numpy as np

HDF5_FILE = "./out/sim-psyduck/dataset-psyduck.hdf5"
RECORD_ID = "record_000003"

with h5py.File(HDF5_FILE, "r") as db_conn:
    rec = db_conn[RECORD_ID]
    demo_tree = pickle.loads(rec["input/tree"][...].tobytes())
    tree_height = rec["input/tree_height"][()]
    present = rec["input/present"][()]
    epidemic_duration = rec["output/parameters/epidemic_duration"][()]
    sampling_prop_change_times = rec["output/parameters/sampling_prop/change_times"][()]

mr_ca_time = present - tree_height
change_times_tree = sampling_prop_change_times - mr_ca_time

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
Phylo.draw(demo_tree, do_show=False, axes=ax1)

for change_time in change_times_tree:
    ax1.axvline(change_time, color="tomato", linestyle="--", linewidth=1.0, alpha=0.8)

ax1.set_title("Sampling proportion change times")
ax1.set_xlim((-mr_ca_time) * 1.1, (epidemic_duration - mr_ca_time) * 1.1)
ax1.axvline(-mr_ca_time, color="red", linestyle="-", linewidth=1.2, alpha=0.9)
ax1.axvline(epidemic_duration - mr_ca_time, color="red", linestyle="-", linewidth=1.2, alpha=0.9)

x = np.linspace(0.0, 1.0, 200)
uniform_pdf = np.where((x >= 0.3) & (x <= 0.7), 1.0 / 0.6, 0.0)
ax2.plot(x, uniform_pdf, color="blue")
ax2.set_title("Uniform(0.3, 0.7) distribution")
ax2.set_xlabel("x")
ax2.set_ylabel("density")
if sampling_prop_change_times.size > 0 and epidemic_duration != 0:
    normalized_change_time = sampling_prop_change_times[0] / epidemic_duration
    ax2.axvline(
        normalized_change_time, color="blue", linestyle="--", linewidth=1.2, alpha=0.9
    )

fig.tight_layout()
fig.savefig("./out/sim-psyduck/plots/demo-tree-change-times.png")
