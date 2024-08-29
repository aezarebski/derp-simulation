# packages
import h5py
import json
import os
import pandas as pd
import pickle
import plotly.io as pio
import plotnine as p9
import numpy as np
import random



###################################################
# USER INPUT: config to use


# CONFIG_JSON = "config/simulation-charmander.json"
# CONFIG_JSON = "config/simulation-charmeleon.json"
# CONFIG_JSON = "config/simulation-charizard.json"
# CONFIG_JSON = "config/simulation-charmeleon-temporal.json"
CONFIG_JSON = "config/debugging-measurement-times.json"

###################################################




# constants
with open(CONFIG_JSON, "r") as file:
    CONFIG = json.load(file)
SIM_DIR = f"out/{CONFIG['simulation_name']}/simulation/remaster"
SIM_PICKLE_DIR = f"out/{CONFIG['simulation_name']}/simulation/pickle"
DB_PATH = f"out/{CONFIG['simulation_name']}/{CONFIG['output_hdf5']}"
PLOT_DIR = f"out/{CONFIG['simulation_name']}/plots"
if not os.path.exists(PLOT_DIR):
    os.makedirs(PLOT_DIR)
    
MAX_NUM_SIMS_TO_PLOT = 50





# Read all temporal data out

def _read_temporal_data(key, db_conn):
    return {
        "key": key,
        "key_num": int(key.split("_")[1]),
        "temporal_data": db_conn[f"{key}/output/parameters/temporal_measurements"][()],
    }
    
DB_CONN = h5py.File(DB_PATH, "r")
all_temporal_data = [
    _read_temporal_data(key, DB_CONN) for key in DB_CONN.keys() if key.startswith("record")
]
DB_CONN.close()


# Pack a sample of simulations into a dataframe in order to plot
inds_of_sims_to_plot = random.sample(range(len(all_temporal_data)),
        (MAX_NUM_SIMS_TO_PLOT if len(all_temporal_data) > MAX_NUM_SIMS_TO_PLOT else len(all_temporal_data)))
temporal_data_sample_df = pd.concat([pd.DataFrame(
                                      {"key_num": all_temporal_data[sim_ind]["key_num"],
                                      "measurement_times": all_temporal_data[sim_ind]["temporal_data"]["measurement_times"],
                                      "prevalence": all_temporal_data[sim_ind]["temporal_data"]["prevalence"],
                                      "cumulative": all_temporal_data[sim_ind]["temporal_data"]["cumulative"],
                                      "reproductive_number": all_temporal_data[sim_ind]["temporal_data"]["reproductive_number"]})
                                      for sim_ind in inds_of_sims_to_plot])



# Plot prevalence curves over time for the sampled simulations
prevalence_p9 = (
    p9.ggplot()
    + p9.geom_point(data = temporal_data_sample_df,
    mapping=p9.aes(x="measurement_times", y="prevalence", group = "key_num"), color = "cornflowerblue")
    + p9.geom_line(data = temporal_data_sample_df,
    mapping=p9.aes(x="measurement_times", y="prevalence", group = "key_num"), size = 0.2, color = "cornflowerblue")
    + p9.labs(x="Forward time since epidemic start",y="Prevalence")
    + p9.scale_y_log10()
)

prevalence_p9.save(f"{PLOT_DIR}/prevalence_over_time.png", width=10, height=10, dpi=300)
prevalence_p9.save(f"{PLOT_DIR}/prevalence_over_time.svg", width=10, height=10, dpi=300)




# Plot cumulative infection curves over time for the sampled simulations
cumulative_p9 = (
    p9.ggplot()
    + p9.geom_point(data = temporal_data_sample_df,
    mapping=p9.aes(x="measurement_times", y="cumulative", group = "key_num"), color = "darkolivegreen")
    + p9.geom_line(data = temporal_data_sample_df,
    mapping=p9.aes(x="measurement_times", y="cumulative", group = "key_num"), size = 0.2, color = "darkolivegreen")
    + p9.labs(x="Forward time since epidemic start",y="Cumulative infections")
    + p9.scale_y_log10()
)

cumulative_p9.save(f"{PLOT_DIR}/cumulative_over_time.png", width=10, height=10, dpi=300)
cumulative_p9.save(f"{PLOT_DIR}/cumulative_over_time.svg", width=10, height=10, dpi=300)



# Plot reproduction number curves over time for the sampled simulations
rzero_p9 = (
    p9.ggplot()
    + p9.geom_point(data = temporal_data_sample_df,
    mapping=p9.aes(x="measurement_times", y="reproductive_number", group = "key_num"), color = "darkorange")
    + p9.geom_line(data = temporal_data_sample_df,
    mapping=p9.aes(x="measurement_times", y="reproductive_number", group = "key_num"), size = 0.2, color = "darkorange")
    + p9.labs(x="Forward time since epidemic start",y="Reproduction number")
)

rzero_p9.save(f"{PLOT_DIR}/rzero_over_time.png", width=10, height=10, dpi=300)
rzero_p9.save(f"{PLOT_DIR}/rzero_over_time.svg", width=10, height=10, dpi=300)
