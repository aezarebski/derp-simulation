from Bio import Phylo
from concurrent.futures import ThreadPoolExecutor, as_completed
from lxml import etree
import os
import pickle
import subprocess


REMASTER_XML = "src/remaster-template.xml"
NUM_WORKERS = 4
NUM_SIMS = 2
SIM_DIR = "out/simulation/remaster"
SIM_PICKLE_DIR = "out/simulation/pickle"


def _update_attr(root, xpath: str, attr: str, val) -> None:
    """Update the attribute of an element in an XML tree."""
    tmp = root.xpath(xpath)[0]
    tmp.attrib[attr] = str(val)
    return None


def random_remaster_parameters():
    parameters = {}
    parameters["epidemic_duration"] = 10
    parameters["num_changes"] = 3
    parameters["birth_rate"] = {"values": [0.3, 0.3, 0.3], "change_times": [3.0, 6.0]}
    parameters["death_rate"] = {"values": [0.1, 0.1, 0.1], "change_times": [3.0, 6.0]}
    parameters["sampling_rate"] = {"values": [0.1], "change_times": []}
    return parameters


def write_simulation_xml(simulation_xml, parameters):
    remaster_xml_obj = etree.parse(REMASTER_XML)
    b = remaster_xml_obj.getroot()
    _update_attr(
        b,
        ".//reaction[@id='lambdaReaction']",
        "rate",
        " ".join([str(br) for br in parameters["birth_rate"]["values"]]),
    )
    if parameters["birth_rate"]["change_times"]:
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
    if parameters["death_rate"]["change_times"]:
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
    if parameters["sampling_rate"]["change_times"]:
        _update_attr(
            b,
            ".//reaction[@id='psiReaction']",
            "changeTimes",
            " ".join([str(ct) for ct in parameters["sampling_rate"]["change_times"]]),
        )
    _update_attr(b, ".//trajectory", "maxTime", parameters["epidemic_duration"])
    _update_attr(
        b, ".//logger[@mode='tree']", "fileName", simulation_xml.replace(".xml", ".tree")
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
        f"{SIM_DIR}/{sim_num:06d}.xml"
        for sim_num in range(1, num_sims + 1)
    ]
    for sim_xml, params in zip(simulation_xml_list, params_list):
        write_simulation_xml(sim_xml, params)

    run_beast2_simulations_parallel(simulation_xml_list, num_jobs=NUM_WORKERS)
    if not os.path.exists(SIM_PICKLE_DIR):
        os.makedirs(SIM_PICKLE_DIR)

    # NOTE we only record a single simulation per iteration to avoid
    # memory issues with large trees.
    for ix, sim_xml, params in zip(
        range(1, num_sims + 1), simulation_xml_list, params_list
    ):
        result = {
            "parameters": params,
            "simulation_xml": sim_xml,
            "simulation_results": read_simulation_results(sim_xml),
        }
        simulation_pickle = f"{SIM_PICKLE_DIR}/{ix:06d}.pickle"
        with open(simulation_pickle, "wb") as f:
            pickle.dump(result, f)


def main():
    run_simulations(NUM_SIMS)
    # record_simulations(NUM_SIMS)


if __name__ == "__main__":
    main()


# import json
# foo = open("out/simulation/pickle/000001.pickle", "rb")
# foobar = pickle.load(foo)
# foobar["simulation_results"]["tree"]  = str(pickle.dumps( foobar["simulation_results"]["tree"] ))
# with open("demo-obj.json", "w") as f:
#     # pretty print the dictionary called foobar
#     json.dump(foobar, f, indent=4)
