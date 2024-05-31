# [[file:vis.org::*Packages][Packages:1]]
import h5py
import json
import os
import pandas as pd
import pickle
import plotly.io as pio
import plotnine as p9
# Packages:1 ends here

# [[file:vis.org::*Read the configuration][Read the configuration:1]]
CONFIG_JSON = "config/debugging.json"

with open(CONFIG_JSON, "r") as file:
    CONFIG = json.load(file)
SIM_DIR = f"out/{CONFIG['simulation-name']}/simulation/remaster"
SIM_PICKLE_DIR = f"out/{CONFIG['simulation-name']}/simulation/pickle"
DB_PATH = f"out/{CONFIG['simulation-name']}/{CONFIG['output-hdf5']}"
PLOT_DIR = f"out/{CONFIG['simulation-name']}/plots"
if not os.path.exists(PLOT_DIR):
    os.makedirs(PLOT_DIR)
# Read the configuration:1 ends here

# [[file:vis.org::*Setting up dataframes from the simulated data][Setting up dataframes from the simulated data:1]]
def _record_summary(key, db_conn):
    return {
        "key": key,
        "key_num": int(key.split("_")[1]),
        "present_time": db_conn[f"{key}/input/present"][()],
        "cumulative_infections": db_conn[f"{key}/output/present_cumulative"][()],
        "prevalence": db_conn[f"{key}/output/present_prevalence"][()],
        "r0_change_times": db_conn[f"{key}/output/parameters/r0/change_times"][()],
        "r0_values": db_conn[f"{key}/output/parameters/r0/values"][()],
        "tree_height": db_conn[f"{key}/input/tree_height"][()],
    }


DB_CONN = h5py.File(DB_PATH, "r")
data_dicts = [
    _record_summary(key, DB_CONN) for key in DB_CONN.keys() if key.startswith("record")
]
DB_CONN.close()

foo = [
    [
        {"key_num": dd["key_num"], "change_ix": ix, "time": t}
        for ix, t in enumerate(dd["r0_change_times"].tolist())
    ]
    for dd in data_dicts
]
bar = []
for f in foo:
    bar.extend(f)
change_time_df = pd.DataFrame(bar)

tree_times_df = pd.DataFrame(
    [
        {
            "key_num": dd["key_num"],
            "present": dd["present_time"],
            "tmrca": dd["present_time"] - dd["tree_height"],
        }
        for dd in data_dicts
    ]
)
# Setting up dataframes from the simulated data:1 ends here

# [[file:vis.org::*Plot: random selection of R0 functions][Plot: random selection of R0 functions:1]]

# Take a random subset of 10 of the elements from data_dicts
tmp = pd.DataFrame(data_dicts).sample((50 if len(data_dicts) > 50 else len(data_dicts)))

def _r0_plot_df(subset_data_dicts_df, key_num):
    global CONFIG
    max_sim_duration = CONFIG["simulation-hyperparameters"]["duration-range"][-1]
    foo = tmp[tmp.key_num == key_num].r0_change_times.item().tolist()
    foo.insert(0, 0)
    foo.insert(len(foo), max_sim_duration)
    bar = tmp[tmp.key_num == key_num].r0_values.item().tolist()
    bar.insert(len(bar), bar[-1])
    return pd.DataFrame({"time": foo, "r0": bar, "key_num": key_num})

r0_plot_df = pd.concat([_r0_plot_df(tmp, k) for k in tmp.key_num.tolist()])

r0_trajectories_p9 = (p9.ggplot() + p9.geom_step(data=r0_plot_df, mapping=p9.aes(x="time", y="r0", group="key_num"), alpha = 0.5) + p9.theme_bw())
r0_trajectories_p9.save(f"{PLOT_DIR}/r0_trajectories.png", width=10, height=10, dpi=300)
r0_trajectories_p9.save(f"{PLOT_DIR}/r0_trajectories.svg", width=10, height=10, dpi=300)
# Plot: random selection of R0 functions:1 ends here

# [[file:vis.org::*Simulation timelines][Simulation timelines:1]]
timelines_p9 = (
    p9.ggplot()
    + p9.geom_hline(
        data=change_time_df,
        mapping=p9.aes(yintercept="key_num"),
        color="gray",
        linetype="dashed",
    )
    + p9.geom_point(data=change_time_df, mapping=p9.aes(x="time", y="key_num"))
    + p9.geom_point(
        data=tree_times_df, mapping=p9.aes(x="present", y="key_num"), color="red"
    )
    + p9.geom_point(
        data=tree_times_df, mapping=p9.aes(x="tmrca", y="key_num"), color="blue"
    )
    + p9.theme_bw()
)

timelines_p9.save("out/plots/timelines.png", width=10, height=10, dpi=300)
# Simulation timelines:1 ends here

# [[file:vis.org::*Distribution of last sequence times][Distribution of last sequence times:1]]
last_seq_hist_p9 = (
    p9.ggplot()
    + p9.geom_histogram(
        data=tree_times_df,
        mapping=p9.aes(x="present"),
        bins=20,
    )
    + p9.geom_vline(xintercept=[20, 40], linetype="dashed", color="red")
    + p9.scale_x_continuous(limits=(0, 42), name="Time of last sequence")
    + p9.theme_bw()
    + p9.theme(axis_title_y=p9.element_blank())
)

# Save the plot
last_seq_hist_p9.save("out/plots/last_seq_hist.png", width=10, height=10, dpi=300)
# Distribution of last sequence times:1 ends here
