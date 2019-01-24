
import argparse, os, copy, errno, csv

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

info_by_treatment = {}
for sol in solutions:
    treatment = sol[header_lu["treatment"]]
    if not treatment in info_by_treatment:
        info_by_treatment[treatment] = {"solutions_found":0, "total_runs":0, "min_solution_length":"NONE"}
    if (sol[header_lu["solution_found"]] == "True"):
        info_by_treatment[treatment]["solutions_found"] += 1
        if info_by_treatment[treatment]["min_solution_length"] == "NONE":
            info_by_treatment[treatment]["min_solution_length"] = sol[header_lu["solution_length"]]
        elif info_by_treatment[treatment]["min_solution_length"] > sol[header_lu["solution_length"]]:
            info_by_treatment[treatment]["min_solution_length"] = sol[header_lu["solution_length"]]

    info_by_treatment[treatment]["total_runs"] += 1
    info_by_treatment[treatment]["uses_cohorts"] = sol[header_lu["uses_cohorts"]]
    info_by_treatment[treatment]["problem"] = sol[header_lu["problem"]]

solutions_summary = "treatment,test_mode,problem,uses_cohorts,solutions_found,total_runs,min_solution_size\n"
for treatment in info_by_treatment:
    info = info_by_treatment[treatment]
    print(treatment)
    test_mode = treatment.split("__")[1].replace("TESTS_", "")
    solutions_summary += ",".join(map(str, [treatment, test_mode, info["problem"], info["uses_cohorts"], info["solutions_found"], info["total_runs"], info["min_solution_length"]])) + "\n"

new_name = fpath.split("/")[-1].strip(".csv") + "__solutions_summary.csv"
with open(new_name, "w") as fp:
    fp.write(solutions_summary)   
    
