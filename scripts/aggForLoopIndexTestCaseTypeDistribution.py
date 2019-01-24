import argparse, os, copy, errno, csv


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

def main():
    parser = argparse.ArgumentParser(description="Data aggregation script.")
    parser.add_argument("data_directory", type=str, help="Target experiment directory.")
    parser.add_argument("dump_directory", type=str, help="Where to dump this?")
    parser.add_argument("update", type=int, help="max update to look at distribution")

    args = parser.parse_args()

    data_directory = args.data_directory
    dump = args.dump_directory
    
    print("Pulling test case distributions in for-loop-index problem from {}".format(data_directory))

    mkdir_p(dump)

    # Get a list of all runs
    runs = [d for d in os.listdir(data_directory) if d.strip("PROBLEM_").split("__")[0] == "for-loop-index"]
    runs.sort()

    out_content = "treatment,run_id,solution_found,test_out_len_1,total_tests\n"

    for run in runs:
        print("Run: {}".format(run))
        run_dir = os.path.join(data_directory, run)
        run_id = run.split("__")[-1]
        run = "__".join(run.split("__")[:-1])
        treatment = run

        pop_dir = os.path.join(run_dir, "output", "pop_{}".format(args.update))
        if not os.path.isdir(pop_dir):
            print("  pop_{} dir not found".format(args.update))
            continue
        
        pop_fpath = os.path.join(pop_dir, "test_pop_{}.csv".format(args.update))
        file_content = None
        with open(pop_fpath, "r") as fp:
            csvreader = csv.reader(fp, delimiter=',', quotechar='"')
            next(csvreader, None)
            tests = [list(map(int,row[-1].split(","))) for row in csvreader]
        
        # Categorize test type distribution
        test_type__out_len_1 = 0
        total_tests = len(tests)

        for test in tests:
            start = test[0]
            end = test[1]
            step = test[2]
            if ((start + step) >= end): test_type__out_len_1 += 1

        # Figure out if this run produced a solution
        file_content = None
        run_sols = os.path.join(run_dir, "output", "solutions.csv")
        with open(run_sols, "r") as fp:
            file_content = fp.read().strip().split("\n")

        header = file_content[0].split(",")
        header_lu = {header[i].strip():i for i in range(0, len(header))}
        file_content = file_content[1:]    
        solutions = [l for l in csv.reader(file_content, quotechar='"', delimiter=',', quoting=csv.QUOTE_ALL, skipinitialspace=True)]
        
        sol_found = "1" if len(solutions) > 0 else "0"

        #out_content = "treatment,run_id,solution,test_out_len_1,total_tests\n"
        out_content += ",".join(list(map(str, [treatment, run_id, sol_found, test_type__out_len_1, total_tests]))) + "\n"
    
    with open(os.path.join(dump, "for-loop-index_test_distribution__update_{}.csv".format(args.update)), "w") as fp:
        fp.write(out_content)





        

        
            

if __name__ == "__main__":
    main()



