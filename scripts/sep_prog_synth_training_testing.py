import argparse, os, copy, errno, csv

aggregator_dump = "./aggregated_data"

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
    parser = argparse.ArgumentParser(description="Script to separate training and testing data from original programming synthesis benchmark example csv files.")
    parser.add_argument("examples_directory", type=str, help="Examples directory.")
    parser.add_argument("dump_directory", type=str, help="Dump directory.")

    args = parser.parse_args()
    
    ex_dir = args.examples_directory
    dump_dir = args.dump_directory

    mkdir_p(dump_dir)

    example_sets = [f for f in os.listdir(ex_dir) if ".csv" in f]
    for fname in example_sets:
        problem = "-".join(fname.replace(".csv", "").split("-")[1:])
        print ("Processing examples for {} (problem: {})".format(fname, problem))

        training_content = []
        testing_content = []
        testing = False
        
        ex_fpath = os.path.join(ex_dir, fname)
        content = None
        with open(ex_fpath, "r") as fp:
            content = fp.read().split("\n")
            #content = [l for l in csv.reader(fp, quotechar='"', delimiter=',', quoting=csv.QUOTE_ALL, skipinitialspace=True)]
        
        for line in content:
            if len(line) >= 13:
                if line[:13] == "test_input_1,":
                    print(">> Found test set.")
                    testing = True
            if testing:
                testing_content.append(line)
            else:
                training_content.append(line)

        ########################################################################
        ########################################################################
        # Extra, problem-specific processing.
        ########################################################################
        if problem == "for-loop-index":
            processed_testing_content = []
            for line in testing_content:
                if "," in line: 
                    processed_testing_content.append(line)
                    continue
                processed_testing_content[-1] += "," + line
            testing_content = processed_testing_content

            processed_training_content = []
            for line in training_content:
                if "," in line: 
                    processed_training_content.append(line)
                    continue
                processed_training_content[-1] += "," + line
            training_content = processed_training_content
        ########################################################################
        
        testing_fname = "testing-examples-{}.csv".format(problem)
        training_fname = "training-examples-{}.csv".format(problem)
        with open(os.path.join(dump_dir, testing_fname), "w") as fp:
            fp.write("\n".join(testing_content).strip()) 
        with open(os.path.join(dump_dir, training_fname), "w") as fp:
            fp.write("\n".join(training_content).strip()) 
        

if __name__ == "__main__":
    main()