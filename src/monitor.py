import os
import re
import json
import time
import tqdm

# CONFIG_JSON = "config/simulation-bulbasaur.json"
CONFIG_JSON = os.sys.argv[1]


with open(CONFIG_JSON, "r") as file:
    CONFIG = json.load(file)

NUM_SIMS = CONFIG["num_simulations"]
SIM_DIR = f"out/{CONFIG['simulation_name']}/simulation/remaster"
SIM_PICKLE_DIR = f"out/{CONFIG['simulation_name']}/simulation/pickle"

xml_files = [f for f in os.listdir(SIM_DIR) if re.match(r"[0-9]{6}.xml", f)]
assert (
    len(xml_files) == NUM_SIMS
), f"Expected {NUM_SIMS} XML files, found {len(xml_files)}"

print("=============================")
print("Monitoring BEAST2 simulations")
print("=============================\n")
pbar = tqdm.tqdm(
    total=NUM_SIMS,
    ncols=100,
    miniters=1,
    mininterval=0.1,
    smoothing=0,
    bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt} [Elapsed: {elapsed} Remaining: {remaining}]",
)
old_num_tree_files = 0
new_num_tree_files = 0
while new_num_tree_files < NUM_SIMS:
    new_num_tree_files = len(
        [f for f in os.listdir(SIM_DIR) if re.match(r"[0-9]{6}.tree", f)]
    )
    pbar.update(new_num_tree_files - old_num_tree_files)
    old_num_tree_files = new_num_tree_files
    time.sleep(0.1)
pbar.close()

print("========================")
print("Monitoring tree pickling")
print("========================\n")
pbar = tqdm.tqdm(
    total=NUM_SIMS,
    ncols=100,
    miniters=1,
    mininterval=0.1,
    smoothing=0,
    bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt} [Elapsed: {elapsed} Remaining: {remaining}]",
)
old_num_pickle_files = 0
new_num_pickle_files = 0
while new_num_pickle_files < NUM_SIMS:
    new_num_pickle_files = len(
        [f for f in os.listdir(SIM_PICKLE_DIR) if re.match(r"[0-9]{6}.pickle", f)]
    )
    pbar.update(new_num_pickle_files - old_num_pickle_files)
    old_num_pickle_files = new_num_pickle_files
    time.sleep(0.1)
pbar.close()
