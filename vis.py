# [[file:vis.org::*Packages and global variables][Packages and global variables:1]]
import h5py
import json
import os
import pandas as pd
import pickle
import plotly.io as pio
import plotnine as p9

# DB_PATH = "out/debugging/dataset-demo.hdf5"
DB_PATH = "out/medium-size/dataset-medium-size.hdf5"
if not os.path.exists("out/plots"):
    os.makedirs("out/plots")
# Packages and global variables:1 ends here

# [[file:vis.org::*Setting up dataframes from the simulated data][Setting up dataframes from the simulated data:1]]
def _record_summary(key, db_conn):
    return {
        "key": key,
        "key_num": int(key.split("_")[1]),
        "present_time": db_conn[f"{key}/input/present"][()],
        "cumulative_infections": db_conn[f"{key}/output/present_cumulative"][()],
        "prevalence": db_conn[f"{key}/output/present_prevalence"][()],
        "r0_change_times": db_conn[f"{key}/output/parameters/r0/change_times"][()],
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
