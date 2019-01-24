'''
Script: agg_min_correct_networks.py

For each run, grab smallest correct solution network. If run has none, report none.

'''

import argparse, os, copy, errno, csv


def main():
    parser = argparse.ArgumentParser(description="Data aggregation script.")
    parser.add_argument("data_directory", type=str, help="Target experiment directory.")
    parser.add_argument("-r", "--replicates", type=int, help="Expected number of replicates.")

    args = parser.parse_args()

    data_directory = args.data_directory
    
    # Get a list of all runs
    runs = [d for d in os.listdir(data_directory) if "__" in d]
    runs.sort()

    run_ids_by_treatment = {}

    undone_content = "treatment,run_id,target_gen,final_update\n"
    cnt = 0
    total = 0
    for run in runs:
        print("Run: {}".format(run))
        run_dir = os.path.join(data_directory, run)
        run_id = run.split("__")[-1]
        run_name = "__".join(run.split("__")[:-1])
        run_log = None
        with open(os.path.join(run_dir, "run.log"), "r") as fp:
            run_log = fp.read().strip().split("\n")
        target_gen = None
        for line in run_log:
            if "set GENERATIONS" in line:
                target_gen = line.split(" ")[2]
                break
        if target_gen == None:
            print("Failed to find target generations in run log!")
            exit(-1)
        # Did this run finish?
        final_line = run_log[-1]
        finished = False
        if "Update: {};".format(target_gen) in final_line:
            finished = True
            print("  ==> Finished!")
        else:
            finished = False
            print("  ==> Not Finished!")
            cnt+=1
        final_update = final_line.split(";")[0].split(" ")[-1]
        
        if not finished:
            undone_content += ",".join([run_name, run_id, target_gen, final_update]) + "\n"

        if not (run_name in run_ids_by_treatment): run_ids_by_treatment[run_name] = []
        run_ids_by_treatment[run_name].append(int(run_id))
            
        total += 1

    with open("undone.csv", "w") as fp:
        fp.write(undone_content)
    print ("Runs not done: " + str(cnt))

    # Try to figure out if there are any jobs that failed to get submitted.
    if args.replicates != None:

        treatments_with_missing_runs = []
        for treatment in run_ids_by_treatment:
            print("Runs from treatment {}: {}".format(treatment, len(run_ids_by_treatment[treatment])))
            if len(run_ids_by_treatment[treatment]) < args.replicates:
                treatments_with_missing_runs.append(treatment)
    
        for treatment in treatments_with_missing_runs:
            run_ids = run_ids_by_treatment[treatment]
            run_ids.sort()
            min_run = min(run_ids)
            max_run = max(run_ids)
            expected_sequence = [i for i in range(min_run, max_run+1)]
            missing_runs = [rid for rid in expected_sequence if not (rid in run_ids)]
            print("== Treatment: {} ==".format(treatment))
            print("  MIN RID: {}".format(min_run))
            print("  MAX RID: {}".format(max_run))
            print("  Best guess at missing runs: {}".format(str(missing_runs)))

if __name__ == "__main__":
    main()