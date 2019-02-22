'''
Script: agg_min_correct_networks.py

For each run, grab smallest correct solution network. If run has none, report none.

'''

import argparse, os, copy, errno, csv

default_update = 10000

problem_whitelist = ["grade", "number-io", "median", "smallest", "small-or-large", "compare-string-lengths"]

cohort_configs = {
    "CN_128__CS_4": "cn128:cs4",
    "CN_16__CS_32": "cn16:cs32",
    "CN_1__CS_512": "cn1:cs512",
    "CN_256__CS_2": "cn256:cs2",
    "CN_2__CS_256": "cn2:cs256",
    "CN_32__CS_16": "cn32:cs16",
    "CN_4__CS_128": "cn4:cs128",
    "CN_64__CS_8": "cn64:cs8",
    "CN_8__CS_64": "cn8:cs64",
    "CN_1__CS_100": "cn1:cs100",
    "CN_2__CS_50": "cn2:cs50",
    "CN_4__CS_25": "cn4:cs25",
    "CN_10__CS_10": "cn10:cs10",
    "CN_20__CS_5": "cn20:cs5"
}

def mkdir_p(path):
    """
    This is functionally equivalent to the mkdir -p [fname] bash command
    """
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise

def FixDownSampledEval(treatment, e, u):
    cohort_config = None
    for thing in cohort_configs:
        if thing in treatment: cohort_config = cohort_configs[thing]
    if cohort_config == None: 
        print("Unrecognized cohort config! Exiting.")
        exit()
   
    # fixing it now.
    if ("SEL_DOWN_SAMPLE_TESTS" in treatment):
        cn = float(cohort_config.split(":")[0][2:])
        cs = float(cohort_config.split(":")[1][2:])
        # u_eval_found = int(sol[header_lu["update_found"]])
        # u_eval_first_sol_found = int(sol[header_lu["update_first_solution_found"]])
        
        # evaluation_found = (cn * cs * cs)*u_eval_found
        # evaluation_first_solution_found = (cn * cs * cs)*u_eval_first_sol_found
        return (cn * cs * cs)*u
    else:
        return e

def main():
    parser = argparse.ArgumentParser(description="Data aggregation script.")
    parser.add_argument("data_directory", type=str, help="Target experiment directory.")
    parser.add_argument("dump_directory", type=str, help="Where to dump this?")
    parser.add_argument("-u", "--update", type=int, help="max update to look for solutions")
    parser.add_argument("-e", "--evaluations", type=int, help="max evaluation to look for solutions")

    args = parser.parse_args()

    data_directory = args.data_directory
    dump = args.dump_directory
    
    print("Pulling smallest network solutions from all runs in {}".format(data_directory))

    mkdir_p(dump)

    # Get a list of all runs
    runs = [d for d in os.listdir(data_directory) if d.strip("PROBLEM_").split("__")[0] in problem_whitelist]
    runs.sort()

    if args.update != None:
        update = args.update   
        print("Looking for best solutions from update {} or earlier.".format(update)) 
        
        solutions_content = "treatment,run_id,problem,uses_cohorts,solution_found,solution_length,update_found,evaluation_found,update_first_solution_found,evaluation_first_solution_found,program\n"
        
        for run in runs:
            print("Run: {}".format(run))
            run_dir = os.path.join(data_directory, run)
            run_id = run.split("__")[-1]
            run = "__".join(run.split("__")[:-1])
            treatment = run
            run_sols = os.path.join(run_dir, "output", "solutions.csv")

            uses_cohorts = "1" if not "SEL_LEX" in treatment else "0"
            problem = run.strip("PROBLEM_").split("__")[0]
            
            file_content = None
            with open(run_sols, "r") as fp:
                file_content = fp.read().strip().split("\n")

            header = file_content[0].split(",")
            header_lu = {header[i].strip():i for i in range(0, len(header))}
            file_content = file_content[1:]
            
            solutions = [l for l in csv.reader(file_content, quotechar='"', delimiter=',', quoting=csv.QUOTE_ALL, skipinitialspace=True)]
            for i in range(0, len(solutions)):
                sol_evaluation = float(solutions[i][header_lu["evaluations"]])
                sol_update = float(solutions[i][header_lu["update"]])
                sol_evaluation = FixDownSampledEval(treatment, e=sol_evaluation, u=sol_update)
                solutions[i][header_lu["evaluations"]] = str(sol_evaluation)
            # Add smallest solution to smallest solution doc
            min_program = None
            sol_found = False
            if len(solutions) > 0:
                # Find the smallest program
                for i in range(0, len(solutions)):
                    sol_update = int(solutions[i][header_lu["update"]])
                    if sol_update > update: continue
                    if min_program == None:
                        min_program = i
                        sol_found = True
                    elif float(solutions[i][header_lu["program_len"]]) < float(solutions[min_program][header_lu["program_len"]]):
                        min_program = i
                        sol_found = True
            
            if sol_found:
                # Record timing info about first solution
                update_first_sol = solutions[0][header_lu["update"]]
                eval_first_sol = solutions[0][header_lu["evaluations"]]
                # Record info about smallest solution
                min_sol = solutions[min_program]
                program_len = min_sol[header_lu["program_len"]]
                update_found = min_sol[header_lu["update"]]
                evaluation_found = min_sol[header_lu["evaluations"]]
                program = min_sol[header_lu["program"]]
            else:
                update_first_sol = "NONE"
                eval_first_sol = "NONE"
                program_len = "NONE"
                update_found = "NONE"
                evaluation_found = "NONE"
                program = "NONE"
            # "treatment,run_id,problem,uses_cohorts,solution_found,solution_length,update_found,program\n"
            solutions_content += ",".join(map(str,[treatment, run_id, problem, uses_cohorts, sol_found, program_len, update_found, evaluation_found, update_first_sol, eval_first_sol, '"{}"'.format(program)])) + "\n"
        with open(os.path.join(dump, "min_programs__update_{}.csv".format(update)), "w") as fp:
            fp.write(solutions_content)


    if args.evaluations != None:
        evaluations = args.evaluations   
        print("Looking for best solutions from evaluations {} or earlier.".format(evaluations)) 
        
        solutions_content = "treatment,run_id,problem,uses_cohorts,solution_found,solution_length,update_found,evaluation_found,update_first_solution_found,evaluation_first_solution_found,program\n"
        
        for run in runs:
            print("Run: {}".format(run))
            run_dir = os.path.join(data_directory, run)
            run_id = run.split("__")[-1]
            run = "__".join(run.split("__")[:-1])
            treatment = run
            run_sols = os.path.join(run_dir, "output", "solutions.csv")

            uses_cohorts = "1" if not "SEL_LEX" in treatment else "0"
            problem = run.strip("PROBLEM_").split("__")[0]
            
            file_content = None
            with open(run_sols, "r") as fp:
                file_content = fp.read().strip().split("\n")

            header = file_content[0].split(",")
            header_lu = {header[i].strip():i for i in range(0, len(header))}
            file_content = file_content[1:]
            
            solutions = [l for l in csv.reader(file_content, quotechar='"', delimiter=',', quoting=csv.QUOTE_ALL, skipinitialspace=True)]
            for i in range(0, len(solutions)):
                sol_evaluation = float(solutions[i][header_lu["evaluations"]])
                sol_update = float(solutions[i][header_lu["update"]])
                sol_evaluation = FixDownSampledEval(treatment, e=sol_evaluation, u=sol_update)
                solutions[i][header_lu["evaluations"]] = str(sol_evaluation)

            # Add smallest solution to smallest solution doc
            min_program = None
            sol_found = False
            
            if len(solutions) > 0:
                # Find the smallest program
                for i in range(0, len(solutions)):
                    sol_evaluation = float(solutions[i][header_lu["evaluations"]])
                    sol_update = float(solutions[i][header_lu["update"]])
                    if sol_evaluation > evaluations: continue
                    if min_program == None:
                        min_program = i
                        sol_found = True
                    elif float(solutions[i][header_lu["program_len"]]) < float(solutions[min_program][header_lu["program_len"]]):
                        min_program = i
                        sol_found = True
            
            if sol_found:
                # Record timing info about first solution
                update_first_sol = solutions[0][header_lu["update"]]
                eval_first_sol = solutions[0][header_lu["evaluations"]]
                # Record info about smallest solution
                min_sol = solutions[min_program]
                program_len = min_sol[header_lu["program_len"]]
                update_found = min_sol[header_lu["update"]]
                evaluation_found = min_sol[header_lu["evaluations"]]
                program = min_sol[header_lu["program"]]
            else:
                update_first_sol = "NONE"
                eval_first_sol = "NONE"
                program_len = "NONE"
                update_found = "NONE"
                evaluation_found = "NONE"
                program = "NONE"
            # "treatment,run_id,problem,uses_cohorts,solution_found,solution_length,update_found,program\n"
            solutions_content += ",".join(map(str,[treatment, run_id, problem, uses_cohorts, sol_found, program_len, update_found, evaluation_found, update_first_sol, eval_first_sol, '"{}"'.format(program)])) + "\n"
        with open(os.path.join(dump, "min_programs__eval_{}.csv".format(evaluations)), "w") as fp:
            fp.write(solutions_content)

if __name__ == "__main__":
    main()