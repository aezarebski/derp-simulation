from Bio import Phylo
from concurrent.futures import ThreadPoolExecutor, as_completed
import datetime
import h5py
import json
from lxml import etree
import numpy as np
import os
import pandas as pd
import pickle
import re
import subprocess
import glob

if len(os.sys.argv) < 2:
    raise Exception("Please provide the path to the configuration file. For example ./config/simulation-charmander.json")

with open(os.sys.argv[1], "r") as file:
    CONFIG = json.load(file)

np.random.seed(CONFIG["seed"])


REMASTER_XML = CONFIG["remaster_xml"]
NUM_WORKERS = CONFIG["num_workers"]
NUM_SIMS = CONFIG["num_simulations"]
SIM_DIR = f"out/{CONFIG['simulation_name']}/simulation/remaster"
SIM_PICKLE_DIR = f"out/{CONFIG['simulation_name']}/simulation/pickle"
DB_PATH = f"out/{CONFIG['simulation_name']}/{CONFIG['output_hdf5']}"
NUM_TEMP_MEASUREMENTS = CONFIG["simulation_hyperparameters"]["num_temp_measurements"]
LIMITED_TIME_SAMPLING = CONFIG["simulation_hyperparameters"].get("limited_time_sampling", False)


def prompt_user(message):
    while True:
        response = input(message).lower()
        if response in ["y", "n"]:
            return response
        else:
            print("Invalid input. Please enter 'y' or 'n'.")


if os.path.exists(DB_PATH):
    response = prompt_user(
        f"File {DB_PATH} already exists. Do you want to delete it and proceed? [y/n]: "
    )
    if response == "y":
        confirmation = prompt_user(
            "Are you sure you want to delete the existing database file? This action cannot be undone. [y/n]: "
        )
        if confirmation == "y":
            os.remove(DB_PATH)
            print(f"File {DB_PATH} has been deleted.")
        else:
            raise Exception("Deletion cancelled.")
    else:
        raise Exception(f"File {DB_PATH} already exists.")


def _update_attr(root, xpath: str, attr: str, val) -> None:
    """Update the attribute of an element in an XML tree."""
    tmp = root.xpath(xpath)[0]
    tmp.attrib[attr] = str(val)
    return None

def random_remaster_parameters():
    """
    Generate random parameters for the remaster model.

    NOTE that this makes use of the global `CONFIG` variable.

    NOTE that we are using a Dirichlet distribution to generate the
    change times (hard coded). This is to avoid the change times being 
    too close together, which is biologically implausible.
    """
    hyperparams = CONFIG["simulation_hyperparameters"]
    p = {}
    if hyperparams["duration_range"]["dist"] == "uniform_int":
        p["epidemic_duration"] = np.random.randint(
            hyperparams["duration_range"]["lower_bound"], 
            hyperparams["duration_range"]["upper_bound"] + 1
        )
    else:
        raise NotImplementedError("Currently, only the uniform (integer) distribution is supported for epidemic duration")
    if hyperparams["num_changes"]["dist"] == "uniform_int":
        p["num_changes"] = np.random.randint(
            hyperparams["num_changes"]["lower_bound"], 
            hyperparams["num_changes"]["upper_bound"] + 1
        )
    else:
        raise NotImplementedError("Currently, only the uniform (integer) distribution is supported for number of changes")
    alpha_param = 3
    cts = (
        p["epidemic_duration"]
        * np.cumsum(np.random.dirichlet([alpha_param] * (p["num_changes"] + 1)))[0:-1]
    )
    p["change_times"] = cts
    # Epidemic parameterisation
    if hyperparams["r0"]["dist"] == "lognormal":
        p["r0"] = {
            "values": np.random.lognormal(
                mean = hyperparams["r0"]["LN_mean"],
                sigma = hyperparams["r0"]["LN_sigma"],
                size=p["num_changes"] + 1
            ),
            "change_times": cts,
        }
    else:
        raise NotImplementedError("Currently, only the lognormal distribution is supported for r0")
    # The following sets up the remaining parameters which depend upon
    # whether there is contemporaneous sampling or not.
    if not hyperparams.get("contemporaneous_sample", False):
        p["contemporaneous_sample"] = False
        return _rand_remaster_params_serial(p, hyperparams)
    else:
        p["contemporaneous_sample"] = True
        return _rand_remaster_params_contemporaneous(p, hyperparams)


def _rand_remaster_params_serial(p, hyperparams):
    # The net removal rate is the total rate at which infected
    # individuals cease to be infectious. I.e. it is the reciprocal of
    # the mean infectious duration. In the canonical parameterisation
    # this is the sum of the sampling and death rate; in the
    # epidemiological parameterisation it is the "net becoming
    # uninfectious rate"
    if hyperparams["net_removal_rate"]["dist"] == "lognormal":
        p["net_removal_rate"] = {
            "values": np.random.lognormal(
                mean = hyperparams["net_removal_rate"]["LN_mean"],
                sigma = hyperparams["net_removal_rate"]["LN_sigma"],
                size = 1
            ),
            "change_times": [],
        }
    else:
        raise NotImplementedError("Currently, only the lognormal distribution is supported for net removal rate")
    if LIMITED_TIME_SAMPLING:
        if hyperparams["sampling_prop"]["dist"] == "uniform":
            p["sampling_prop"] = {
                "values": np.array(
                    [
                        0.0,
                        np.random.uniform(
                            hyperparams["sampling_prop"]["lower_bound"],
                            hyperparams["sampling_prop"]["upper_bound"]
                        ),
                    ]
                ),
                # TODO: this just randomly selects ANY time uniformly - should be more specific
                "change_times": np.array([p["epidemic_duration"]*np.random.uniform()]),
            }
        else:
            raise NotImplementedError("Currently, only the uniform distribution is supported for sampling prop")
    else:
        if hyperparams["sampling_prop"]["dist"] == "uniform":
            p["sampling_prop"] = {
                "values": np.random.uniform(
                    hyperparams["sampling_prop"]["lower_bound"],
                    hyperparams["sampling_prop"]["upper_bound"],
                    size=p["num_changes"] + 1,
                ),
                "change_times": p["change_times"],
            }
        else:
            raise NotImplementedError("Currently, only the uniform distribution is supported for sampling prop")
    # Rate parameterisation
    p["birth_rate"] = {
        "values": p["r0"]["values"] * p["net_removal_rate"]["values"],
        "change_times": p["change_times"],
    }
    p["death_rate"] = {
        "values": p["net_removal_rate"]["values"] * (1 - p["sampling_prop"]["values"]),
        "change_times": p["sampling_prop"]["change_times"],
    }
    p["sampling_rate"] = {
        "values": p["net_removal_rate"]["values"] * p["sampling_prop"]["values"],
        "change_times": p["sampling_prop"]["change_times"]
    }
    return p


def _rand_remaster_params_contemporaneous(p, hyperparams):
    if hyperparams["net_removal_rate"]["dist"] == "lognormal":
        p["net_removal_rate"] = {
            "values": np.random.lognormal(
                mean = hyperparams["net_removal_rate"]["LN_mean"],
                sigma = hyperparams["net_removal_rate"]["LN_sigma"],
                size = 1
            ),
            "change_times": [],
        }
    else:
        raise NotImplementedError("Currently, only the lognormal distribution is supported for net removal rate")
    p["sampling_prop"] = {
        "values": np.array([0]),
        "change_times": np.array([]),
    }
    p["birth_rate"] = {
        "values": p["r0"]["values"] * p["net_removal_rate"]["values"],
        "change_times": p["change_times"],
    }
    p["death_rate"] = {
        "values": p["net_removal_rate"]["values"],
        "change_times": np.array([])
        if not p["net_removal_rate"]["change_times"]
        else p["net_removal_rate"]["change_times"],
    }
    p["sampling_rate"] = {
        "values": np.array([0]),
        "change_times": np.array([]),
    }
    if hyperparams["sampling_prop"]["dist"] == "uniform":
        p["rho"] = {
            "values": np.random.uniform(
                hyperparams["sampling_prop"]["lower_bound"],
                hyperparams["sampling_prop"]["upper_bound"],
                size=1,
            ),
            "change_times": None,
        }
    else:
        raise NotImplementedError("Currently, only the uniform distribution is supported for sampling prop")
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

    if parameters["contemporaneous_sample"]:
        _update_attr(
            b,
            ".//reaction[@id='rhoReaction']",
            "p",
            parameters["rho"]["values"][0],
        )
        _update_attr(
            b,
            ".//reaction[@id='rhoReaction']",
            "times",
            parameters["epidemic_duration"],
        )
    else:
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
                " ".join(
                    [str(ct) for ct in parameters["sampling_rate"]["change_times"]]
                ),
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
        Run a BEAST2 simulation using the provided XML file.

        If the simulation does not finish within 5 minutes, it is
        considered to have timed out.

        $ ./lib/beast/bin/beast -seed 1 -overwrite <simulation_xml>
        """
        print(f"Running simulation: {simulation_xml}")

        # Find BEAST executable - first checks for an installation
        # in this directory (as in a linux system) or in the
        # Applications folder (on a mac)
        beast_executable_linux = "./lib/beast/bin/beast"
        beast_folder_mac = '/Applications/BEAST*'
        beast_fname_mac = 'bin/beast'
        if os.path.exists(beast_executable_linux):
            beast_executable = beast_executable_linux
        elif glob.glob(beast_folder_mac):
            latest_beast_ver_mac = sorted(glob.glob(beast_folder_mac))[-1]
            beast_executable = os.path.join(latest_beast_ver_mac, beast_fname_mac)
        else:
            raise Exception("BEAST2 executable not found. You might find the src/setupbeast2.py script helpful.")
        command = [beast_executable, "-seed", "1", "-overwrite", simulation_xml]

        # If there is a local packages directory, then the command
        # should be ammended to use that as the package directory.
        maybe_package_dir = "./lib/packages"
        if os.path.exists(maybe_package_dir):
            command.insert(1, "-packagedir")
            command.insert(2, maybe_package_dir)

        try:
            result = subprocess.run(
                command,
                check=True,
                capture_output=True,
                text=True,
                timeout=300,
            )
            return result.stdout
        except subprocess.TimeoutExpired:
            return f"BEAST2 simulation for {simulation_xml} timed out."
        except subprocess.CalledProcessError as e:
            return f"""
            ========================================
            Error occurred while running BEAST2 simulation for {simulation_xml}: {e.stderr}
            Command that was attempted is: {command}
            Double check that you have remaster installed where beast can find it.
            ========================================
            """

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


def read_simulation_results(simulation_xml, params):
    sim_xml_obj = etree.parse(simulation_xml)
    sx = sim_xml_obj.getroot()
    tree_file = sx.xpath(".//logger[@mode='tree']")[0].attrib["fileName"]
    traj_file = sx.xpath(".//logger[not(@mode)]")[0].attrib["fileName"]
    is_serial = sx.find(".//reaction[@spec='PunctualReaction']") is None
    tree_generator = Phylo.parse(tree_file, "nexus")
    tree = next(tree_generator).root
    # read the time of the last sequenced sample and the prevalence at
    # that time.
    traj_df = pd.read_csv(traj_file, sep="\t")
    if is_serial:
        psi_df = traj_df[traj_df["population"] == "Psi"]
        psi_df = psi_df[psi_df["value"] == psi_df["value"].max()]
        psi_df = psi_df[psi_df["t"] == psi_df["t"].min()]
        last_psi_time = psi_df["t"].values[0]
        last_rows = traj_df[traj_df["t"] == last_psi_time]
    else:
        last_psi_time = traj_df["t"].max()
        last_rows = traj_df[traj_df["t"] == last_psi_time]
    last_X = last_rows[last_rows["population"] == "X"]["value"].values[0]
    last_Psi = last_rows[last_rows["population"] == "Psi"]["value"].values[0]
    last_Mu = last_rows[last_rows["population"] == "Mu"]["value"].values[0]
    sim_result_dict = {
        "tree": tree,
        "tree_height": max(tree.depths().values()),
        "present": last_psi_time,
        "present_prevalence": last_X,
        "present_cumulative": last_Psi + last_Mu + last_X,
    }

    meas_times = np.sort(
        np.random.uniform(
            low=0.0, high=sim_result_dict["present"], size=NUM_TEMP_MEASUREMENTS
        )
    )
    r0_change_times = pd.Series(params["r0"]["change_times"])

    temp_data_headers = ",".join(
        ["measurement_times", "prevalence", "cumulative", "reproductive_number"]
    )
    temp_data = []

    for time_ind in range(NUM_TEMP_MEASUREMENTS):
        this_meas_time = meas_times[time_ind]

        most_recent_change_time = traj_df[traj_df["t"] <= this_meas_time]["t"].max()
        rows_this_time = traj_df[traj_df["t"] == most_recent_change_time]
        this_X = rows_this_time[rows_this_time["population"] == "X"][
            "value"
        ].values[0]
        this_Psi = rows_this_time[rows_this_time["population"] == "Psi"][
            "value"
        ].values[0]
        this_Mu = rows_this_time[rows_this_time["population"] == "Mu"][
            "value"
        ].values[0]

        prev_meas_this_time = this_X
        cumul_meas_this_time = this_X + this_Mu + this_Psi

        num_r0_changes_so_far = len(
            r0_change_times[r0_change_times <= this_meas_time]
        )
        r0_meas_this_time = params["r0"]["values"][num_r0_changes_so_far]

        temp_data.append(
            (
                this_meas_time,
                prev_meas_this_time,
                cumul_meas_this_time,
                r0_meas_this_time,
            )
        )

        sim_result_dict["temporal_measurements"] = np.rec.fromrecords(
            temp_data, names=temp_data_headers
        )

    return sim_result_dict


def pickle_simulation_result(sim_pickle, sim_xml, params):
    tree_file = os.path.basename(sim_xml).replace(".xml", ".tree")
    if os.path.exists(f"{SIM_DIR}/{tree_file}"):
        result = {
            "parameters": params,
            "simulation_xml": sim_xml,
            "simulation_results": read_simulation_results(sim_xml, params),
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
        # NOTE We skip over cases where the simulation failed to
        # generate a pickle file containing the reconstructed tree.
        if not os.path.exists(pf):
            continue
        num_sims += 1
        ix_str = re.search(r"\d{6}", pf).group(0)
        print(f"Processing record {ix_str}")
        with open(pf, "rb") as f:
            sim = pickle.load(f)
            rec_grp = db_conn.create_group(f"record_{ix_str}")
            rec_grp.attrs["simulation_xml"] = sim["simulation_xml"]
            in_grp = rec_grp.create_group("input")
            in_grp.create_dataset(
                "tree", data=_tree_to_uint8(sim["simulation_results"]["tree"])
            )
            in_grp.create_dataset(
                "tree_height", data=sim["simulation_results"]["tree_height"]
            )
            in_grp.create_dataset("present", data=sim["simulation_results"]["present"])
            out_grp = rec_grp.create_group("output")
            params_grp = out_grp.create_group("parameters")
            params_grp.create_dataset(
                "epidemic_duration", data=sim["parameters"]["epidemic_duration"]
            )
            params_grp.create_dataset(
                "temporal_measurements",
                data=sim["simulation_results"]["temporal_measurements"],
            )

            for key in parameter_keys:
                param_grp = params_grp.create_group(key)
                param_grp.create_dataset(
                    "values", data=sim["parameters"][key]["values"]
                )
                param_grp.create_dataset(
                    "change_times", data=sim["parameters"][key]["change_times"]
                )
            out_grp.create_dataset(
                "present_prevalence",
                data=sim["simulation_results"]["present_prevalence"],
            )
            out_grp.create_dataset(
                "present_cumulative",
                data=sim["simulation_results"]["present_cumulative"],
            )
    db_conn.attrs["num_simulations"] = num_sims
    db_conn.attrs["creation_date"] = datetime.datetime.now().isoformat()
    db_conn.close()


def main():
    sim_pickles = run_simulations(NUM_SIMS)
    create_database(sim_pickles)


if __name__ == "__main__":
    main()
