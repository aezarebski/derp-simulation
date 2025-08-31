import h5py
import matplotlib.pyplot as plt
import squarify

def wall_time_and_label(db, path):
    label = str(int(path.split("_")[-1]))
    wall_time = db[path].attrs["simulation_wall_time"].item()
    return (label, wall_time)

with h5py.File("out/sim-charmander/dataset-charmander.hdf5", "r") as db:
    times_and_labels = [wall_time_and_label(db, path) for path in db.keys()]
    times_and_labels.sort(key=lambda x: x[1])
    labels, times = zip(*times_and_labels)

plt.figure(figsize=(8, 7), dpi=96)
squarify.plot(sizes=times, color=len(times)*["#1b9e77"], pad=True)
plt.axis("off")
plt.title("Simulation Wall Times")
plt.savefig("out/sim-charmander/plots/walltimes.png")
