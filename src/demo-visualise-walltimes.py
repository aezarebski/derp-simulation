import h5py
import matplotlib.pyplot as plt
import squarify

example_name = "charizard"

def wall_time_and_label(db, path):
    label = str(int(path.split("_")[-1]))
    wall_time = db[path].attrs["simulation_wall_time"].item()
    return (label, wall_time)

with h5py.File(f"out/sim-{example_name}/dataset-{example_name}.hdf5", "r") as db:
    times_and_labels = [wall_time_and_label(db, path) for path in db.keys()]
    times_and_labels.sort(key=lambda x: x[1])
    labels, times = zip(*times_and_labels)

plt.figure(figsize=(8, 7), dpi=96)
squarify.plot(sizes=times, color=len(times)*["#1b9e77"], pad=True)
plt.axis("off")
plt.title("Simulation Wall Times")
plt.savefig(f"out/sim-{example_name}/plots/walltimes.png")

# make a histogram of the walltimes with a logarithmic x axis
plt.figure(figsize=(8, 5), dpi=96)
plt.hist(times, bins=20, color="#1b9e77", edgecolor="black")
plt.xscale("log")
plt.xlabel("Wall Time (seconds, log scale)")
plt.ylabel("Number of Simulations")
plt.title("Histogram of Simulation Wall Times")
plt.savefig(f"out/sim-{example_name}/plots/walltime_histogram.png")
