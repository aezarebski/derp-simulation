from Bio import Phylo
from concurrent.futures import ThreadPoolExecutor, as_completed
import datetime
import h5py
from lxml import etree
import numpy as np
import os
import pandas as pd
import pickle
import re
import subprocess


np.random.seed(42)


REMASTER_XML = "src/remaster-template.xml"
NUM_WORKERS = 3
NUM_SIMS = 20
SIM_DIR = "out/simulation/remaster"
SIM_PICKLE_DIR = "out/simulation/pickle"
DB_PATH = "out/dataset-demo.hdf5"


if os.path.exists(DB_PATH):
    # To prevent overwriting existing data throw an exception.
    raise Exception(f"File {DB_PATH} already exists.")


def _update_attr(root, xpath: str, attr: str, val) -> None:
    """Update the attribute of an element in an XML tree."""
    tmp = root.xpath(xpath)[0]
    tmp.attrib[attr] = str(val)
    return None


def random_remaster_parameters():
    """
    Generate random parameters for the remaster model.

    This is where the simulation is configured, so if you want to
    change the model used this is the function you want to edit.

    Note that we are using a Dirichlet distribution to generate the
    change times. This is to avoid the change times being too close
    together, which is biological implausible. Also, to reduce the
    variability in the parameter values, we shrink the values towards
    their mean. This leads to smoother parameter trajectories but
    maintains the average value.
    """
    p = {}
    p["epidemic_duration"] = np.random.randint(20, 40 + 1)
    p["num_changes"] = np.random.randint(1, 2 + 1)
    alpha_param = 3
    # cts = np.sort(np.random.rand(p["num_changes"]) * p["epidemic_duration"])
    cts = (
        p["epidemic_duration"]
        * np.cumsum(np.random.dirichlet([alpha_param] * (p["num_changes"] + 1)))[0:-1]
    )
    p["change_times"] = cts
    # Epidemic parameterisation
    shrink = lambda x: 0.5 * x + (1 - 0.5) * x.mean()
    p["r0"] = {
        "values": shrink(np.random.uniform(1.0, 2.0, size=p["num_changes"] + 1)),
        "change_times": cts,
    }
    p["net_removal_rate"] = {
        "values": shrink(1 / np.random.uniform(2.0, 14.0)),
        "change_times": [],
    }
    p["sampling_prop"] = {
        "values": shrink(np.random.uniform(0.05, 0.50, size=p["num_changes"] + 1)),
        "change_times": cts,
    }
    # Rate parameterisation
    p["birth_rate"] = {
        "values": p["r0"]["values"] * p["net_removal_rate"]["values"],
        "change_times": cts,
    }
    p["death_rate"] = {
        "values": p["net_removal_rate"]["values"] * (1 - p["sampling_prop"]["values"]),
        "change_times": cts,
    }
    p["sampling_rate"] = {
        "values": p["net_removal_rate"]["values"] * p["sampling_prop"]["values"],
        "change_times": cts,
    }
    return p


def write_simulation_xml(simulation_xml, parameters):
    remaster_xml_obj = etree.parse(REMASTER_XML)
    b = remaster_xml_obj.getroot()
    _update_attr(
        b,
        ".//reaction[@id='lambdaReaction']",
        "rate",
        " ".join([str(br) for br in parameters["birth_rate"]["values"]]),
    )
    if parameters["birth_rate"]["change_times"].shape[0] > 0:
        _update_attr(
            b,
            ".//reaction[@id='lambdaReaction']",
            "changeTimes",
            " ".join([str(ct) for ct in parameters["birth_rate"]["change_times"]]),
        )
    _update_attr(
        b,
        ".//reaction[@id='muReaction']",
        "rate",
        " ".join([str(dr) for dr in parameters["death_rate"]["values"]]),
    )
    if parameters["death_rate"]["change_times"].shape[0] > 0:
        _update_attr(
            b,
            ".//reaction[@id='muReaction']",
            "changeTimes",
            " ".join([str(ct) for ct in parameters["death_rate"]["change_times"]]),
        )
    _update_attr(
        b,
        ".//reaction[@id='psiReaction']",
        "rate",
        " ".join([str(sr) for sr in parameters["sampling_rate"]["values"]]),
    )
    if parameters["sampling_rate"]["change_times"].shape[0] > 0:
        _update_attr(
            b,
            ".//reaction[@id='psiReaction']",
            "changeTimes",
            " ".join([str(ct) for ct in parameters["sampling_rate"]["change_times"]]),
        )
    _update_attr(b, ".//trajectory", "maxTime", parameters["epidemic_duration"])
    _update_attr(
        b,
        ".//logger[@mode='tree']",
        "fileName",
        simulation_xml.replace(".xml", ".tree"),
    )
    _update_attr(
        b, ".//logger[not(@mode)]", "fileName", simulation_xml.replace(".xml", ".traj")
    )
    remaster_xml_obj.write(simulation_xml, pretty_print=True)


def run_beast2_simulations_parallel(simulation_xml_list, num_jobs):
    def run_beast2(simulation_xml):
        """
        $ ./lib/beast/bin/beast -seed 1 -overwrite <simulation_xml>
        """
        print(f"Running simulation: {simulation_xml}")
        command = ["./lib/beast/bin/beast", "-seed", "1", "-overwrite", simulation_xml]
        try:
            result = subprocess.run(command, check=True, capture_output=True, text=True)
            return result.stdout
        except subprocess.CalledProcessError as e:
            return f"Error occurred while running BEAST2 simulation for {simulation_xml}: {e.stderr}"

    with ThreadPoolExecutor(max_workers=num_jobs) as executor:
        future_to_xml = {
            executor.submit(run_beast2, xml): xml for xml in simulation_xml_list
        }
        for future in as_completed(future_to_xml):
            xml = future_to_xml[future]
            try:
                data = future.result()
                print(f"Simulation completed for {xml}: {data}")
            except Exception as exc:
                print(f"Simulation generated an exception for {xml}: {exc}")


def read_simulation_results(simulation_xml):
    tree_generator = Phylo.parse(simulation_xml.replace("xml", "tree"), "nexus")
    tree = next(tree_generator).root
    traj_file = simulation_xml.replace("xml", "traj")
    # read the time of the last sequenced sample and the prevalence at
    # that time.
    traj_df = pd.read_csv(traj_file, sep="\t")
    psi_df = traj_df[traj_df["population"] == "Psi"]
    psi_df = psi_df[psi_df["value"] == psi_df["value"].max()]
    psi_df = psi_df[psi_df["t"] == psi_df["t"].min()]
    last_psi_time = psi_df["t"].values[0]
    last_rows = traj_df[traj_df["t"] == last_psi_time]
    last_X = last_rows[last_rows["population"] == "X"]["value"].values[0]
    last_Psi = last_rows[last_rows["population"] == "Psi"]["value"].values[0]
    last_Mu = last_rows[last_rows["population"] == "Mu"]["value"].values[0]
    return {
        "tree": tree,
        "tree_height": max(tree.depths().values()),
        "present": last_psi_time,
        "present_prevalence": last_X,
        "present_cumulative": last_Psi + last_Mu + last_X,
    }


def pickle_simulation_result(sim_pickle, sim_xml, params):
    tree_file = os.path.basename(sim_xml).replace(".xml", ".tree")
    if os.path.exists(f"{SIM_DIR}/{tree_file}"):
        result = {
            "parameters": params,
            "simulation_xml": sim_xml,
            "simulation_results": read_simulation_results(sim_xml),
        }
        with open(sim_pickle, "wb") as f:
            pickle.dump(result, f)
        return sim_pickle
    return None


def run_pickling_parallel(pickle_files, sim_xml_list, params_list):
    with ThreadPoolExecutor(max_workers=NUM_WORKERS) as executor:
        futures = [
            executor.submit(pickle_simulation_result, sim_pickle, sim_xml, params)
            for sim_pickle, sim_xml, params in zip(
                pickle_files, sim_xml_list, params_list
            )
        ]
        completed_files = [future.result() for future in as_completed(futures)]
    return [f for f in completed_files if f is not None]


def run_simulations(num_sims):
    params_list = [random_remaster_parameters() for _ in range(num_sims)]
    if not os.path.exists(SIM_DIR):
        os.makedirs(SIM_DIR)
    sim_xml_list = [
        f"{SIM_DIR}/{sim_num:06d}.xml" for sim_num in range(1, num_sims + 1)
    ]
    for sim_xml, params in zip(sim_xml_list, params_list):
        write_simulation_xml(sim_xml, params)

    run_beast2_simulations_parallel(sim_xml_list, num_jobs=NUM_WORKERS)
    if not os.path.exists(SIM_PICKLE_DIR):
        os.makedirs(SIM_PICKLE_DIR)
    pickle_files = [
        f"{SIM_PICKLE_DIR}/{ix:06d}.pickle" for ix in range(1, num_sims + 1)
    ]
    return run_pickling_parallel(pickle_files, sim_xml_list, params_list)


def _tree_to_uint8(tree):
    return np.frombuffer(pickle.dumps(tree), dtype="uint8")


def create_database(pickle_files):
    db_conn = h5py.File(DB_PATH, "w")
    parameter_keys = [
        "birth_rate",
        "death_rate",
        "sampling_rate",
        "r0",
        "net_removal_rate",
        "sampling_prop",
    ]
    num_sims = 0
    for pf in pickle_files:
        if not os.path.exists(pf):
            continue
        num_sims += 1
        ix_str = re.search(r"\d{6}", pf).group(0)
        print(f"Processing record {ix_str}")
        with open(pf, "rb") as f:
            foobar = pickle.load(f)
            rec_grp = db_conn.create_group(f"record_{ix_str}")
            rec_grp.attrs["simulation_xml"] = foobar["simulation_xml"]
            in_grp = rec_grp.create_group("input")
            in_grp.create_dataset(
                "tree", data=_tree_to_uint8(foobar["simulation_results"]["tree"])
            )
            in_grp.create_dataset(
                "tree_height", data=foobar["simulation_results"]["tree_height"]
            )
            in_grp.create_dataset(
                "present", data=foobar["simulation_results"]["present"]
            )
            out_grp = rec_grp.create_group("output")
            params_grp = out_grp.create_group("parameters")
            params_grp.create_dataset(
                "epidemic_duration", data=foobar["parameters"]["epidemic_duration"]
            )
            for key in parameter_keys:
                param_grp = params_grp.create_group(key)
                param_grp.create_dataset(
                    "values", data=foobar["parameters"][key]["values"]
                )
                param_grp.create_dataset(
                    "change_times", data=foobar["parameters"][key]["change_times"]
                )
            out_grp.create_dataset(
                "present_prevalence",
                data=foobar["simulation_results"]["present_prevalence"],
            )
            out_grp.create_dataset(
                "present_cumulative",
                data=foobar["simulation_results"]["present_cumulative"],
            )
    db_conn.attrs["num_simulations"] = num_sims
    db_conn.attrs["creation_date"] = datetime.datetime.now().isoformat()
    db_conn.close()


def main():
    sim_pickles = run_simulations(NUM_SIMS)
    create_database(sim_pickles)


if __name__ == "__main__":
    main()
