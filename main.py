from Bio import Phylo
from concurrent.futures import ThreadPoolExecutor, as_completed
import h5py
from lxml import etree
import numpy as np
import os
import pickle
import re
import subprocess


np.random.seed(42)


REMASTER_XML = "src/remaster-template.xml"
NUM_WORKERS = 1
NUM_SIMS = 10
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
    """
    p = {}
    p["epidemic_duration"] = np.random.randint(20, 40 + 1)
    p["num_changes"] = np.random.randint(2, 5 + 1)
    cts = np.sort(np.random.rand(p["num_changes"]) * p["epidemic_duration"])
    p["change_times"] = cts
    # Epidemic parameterisation
    p["r0"] = {
        "values": np.random.uniform(1.0, 2.0, size=p["num_changes"] + 1),
        "change_times": cts
    }
    p["net_removal_rate"] = {
        "values": 1 / np.random.uniform(2.0, 14.0),
        "change_times": []
    }
    p["sampling_prop"] = {
        "values": np.random.uniform(0.05, 0.95, size=p["num_changes"] + 1),
        "change_times": cts
    }
    # Rate parameterisation
    p["birth_rate"] = {
        "values": p["r0"]["values"] * p["net_removal_rate"]["values"],
        "change_times": cts
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
    return {"tree": tree}


def run_simulations(num_sims):
    params_list = [random_remaster_parameters() for _ in range(num_sims)]
    if not os.path.exists(SIM_DIR):
        os.makedirs(SIM_DIR)
    simulation_xml_list = [
        f"{SIM_DIR}/{sim_num:06d}.xml" for sim_num in range(1, num_sims + 1)
    ]
    for sim_xml, params in zip(simulation_xml_list, params_list):
        write_simulation_xml(sim_xml, params)

    run_beast2_simulations_parallel(simulation_xml_list, num_jobs=NUM_WORKERS)
    if not os.path.exists(SIM_PICKLE_DIR):
        os.makedirs(SIM_PICKLE_DIR)
    # NOTE we only record a single simulation per iteration to avoid
    # memory issues with large trees.
    pickle_files = [
        f"{SIM_PICKLE_DIR}/{ix:06d}.pickle" for ix in range(1, num_sims + 1)
    ]
    for sim_pickle, sim_xml, params in zip(
        pickle_files, simulation_xml_list, params_list
    ):
        result = {
            "parameters": params,
            "simulation_xml": sim_xml,
            "simulation_results": read_simulation_results(sim_xml),
        }
        with open(sim_pickle, "wb") as f:
            pickle.dump(result, f)
    return pickle_files


def _tree_to_uint8(tree):
    return np.frombuffer(pickle.dumps(tree), dtype="uint8")


def create_database(pickle_files):
    db_conn = h5py.File(DB_PATH, "w")
    parameter_keys = ["birth_rate", "death_rate", "sampling_rate",
                      "r0", "net_removal_rate", "sampling_prop"]
    for pf in pickle_files:
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

    db_conn.close()


def main():
    create_database(run_simulations(NUM_SIMS))


if __name__ == "__main__":
    main()


# foo = h5py.File(DB_PATH, "r")
# pickle.loads(foo['record_000001']['input']['tree'][...].tobytes())
# foo['record_000001']['output']['parameters']['epidemic_duration'][()]
# foo['record_000001']['output']['parameters']['birth_rate']['values'][()]
# foo['record_000001']['output']['parameters']['birth_rate']['change_times'][()]
# foo.close()
