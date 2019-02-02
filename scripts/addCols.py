
import argparse, os, copy, errno, csv

cohort_configs = {
    "CN_128__CS_4": "cn128:cs4",
    "CN_16__CS_32": "cn16:cs32",
    "CN_1__CS_512": "cn1:cs512",
    "CN_256__CS_2": "cn256:cs2",
    "CN_2__CS_256": "cn2:cs256",
    "CN_32__CS_16": "cn32:cs16",
    "CN_4__CS_128": "cn4:cs128",
    "CN_64__CS_8": "cn64:cs8",
    "CN_8__CS_64": "cn8:cs64"
}

sel_modes = {
    "SEL_COHORT_LEX": "cohort lex",
    "SEL_PROG_ONLY_COHORT_LEX": "prog-only cohorts",
    "SEL_DOWN_SAMPLE_TESTS": "sample tests",
    "SEL_TRUNCATED": "truncated lex"
}

parser = argparse.ArgumentParser(description="Data aggregation script.")
parser.add_argument("data_file", type=str, help="Target data file")
args = parser.parse_args()

fpath = args.data_file

file_content = None
with open(fpath, "r") as fp:
    file_content = fp.read().strip().split("\n")

header = file_content[0].split(",")
header_lu = {header[i].strip():i for i in range(0, len(header))}
file_content = file_content[1:]

solutions = [l for l in csv.reader(file_content, quotechar='"', delimiter=',', quoting=csv.QUOTE_ALL, skipinitialspace=True)]

modified_content = ",".join(header) + ",cohort_config,sel_mode\n"
for sol in solutions:
    treatment = sol[header_lu["treatment"]]
    test_mode = treatment.split("__")[1].replace("TESTS_", "")
    
    cohort_config = None
    for thing in cohort_configs:
        if thing in treatment: cohort_config = cohort_configs[thing]
    if cohort_config == None: 
        print("Unrecognized cohort config! Exiting.")
        exit()

    sel_mode = None
    for thing in sel_modes:
        if thing in treatment: sel_mode = sel_modes[thing]
    if sel_mode == None: 
        print("Unrecognized selection mode! Exiting.")
        exit()

    new_line = ",".join([sol[header_lu["treatment"]],sol[header_lu["run_id"]],sol[header_lu["problem"]],sol[header_lu["uses_cohorts"]],sol[header_lu["solution_found"]],sol[header_lu["solution_length"]],sol[header_lu["update_found"]],sol[header_lu["evaluation_found"]],sol[header_lu["update_first_solution_found"]],sol[header_lu["evaluation_first_solution_found"]], "\"" + sol[header_lu["program"]] + "\"", cohort_config, sel_mode]) + "\n"
    modified_content += new_line
    
with open(fpath.replace(".csv", "__modified.csv"), "w") as fp:
    fp.write(modified_content)