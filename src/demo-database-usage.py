from Bio import Phylo
import h5py
import pickle
import matplotlib.pyplot as plt
import numpy as np

hdf5_file = "../out/sim-charmander/dataset-charmander.hdf5"

db_conn = h5py.File(hdf5_file)

demo_tree = pickle.loads(db_conn['record_000001/input/tree'][...].tobytes())
fig, ax = plt.subplots()
Phylo.draw(demo_tree, do_show=False, axes=ax)
fig.savefig('../out/sim-charmander/plots/demo-tree.png')

measurements = db_conn['record_000001/output/parameters/temporal_measurements'][...]
column_names = measurements.dtype.names
np.savetxt('../out/sim-charmander/demo-measurements.csv',
           measurements, delimiter=',',
           header=','.join(column_names))

db_conn.close()
