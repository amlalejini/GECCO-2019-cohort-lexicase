#Will list all of the incomplete id's that need to finish running
#SEED_ORIGINAL & TREATMENTS  will need to be adjusted based on the problems, treatments, and seeds that the project requires. 
#Will also need to handle RANGE if different from the expected results!
#
#Input 1: Directory where all the data is located, expecting to see undone.csv and undone.txt | cntUndone.py > undone.txt
#Input 2: Directory where the data will be placed
#Input 3: Number of pairings of CN,CS we are seeing, expecting n, where we have [2^0, 2^1, ..., 2^n]
#
#Output : Will create a text file with all the missing IDs sepereated by problem-treatment!
#
#python3

import argparse, os, copy, errno
import pandas as pd

#Will hold all original seeds in the scipts for original run!
SEED_ORIGINAL = {'for-loop-index':{'COHORT_LEX':5000, 'PROG_ONLY_COHORT_LEX':6000, 'DOWN_SAMPLE_TESTS':7000, 'TRUNCATED':8000},
                'median':{'COHORT_LEX':9000, 'PROG_ONLY_COHORT_LEX':10000, 'DOWN_SAMPLE_TESTS':11000, 'TRUNCATED':12000}, 
                'smallest':{'COHORT_LEX':13000, 'PROG_ONLY_COHORT_LEX':14000, 'DOWN_SAMPLE_TESTS':15000, 'TRUNCATED':16000},
                'small-or-large':{'COHORT_LEX':1000, 'PROG_ONLY_COHORT_LEX':2000, 'DOWN_SAMPLE_TESTS':3000, 'TRUNCATED':4000},
                'sum-of-squares':{'COHORT_LEX':17000, 'PROG_ONLY_COHORT_LEX':18000, 'DOWN_SAMPLE_TESTS':19000, 'TRUNCATED':20000}}

#Will hold all missing/incomplete ids!
TREATMENTS = {'for-loop-index':{'COHORT_LEX':[], 'PROG_ONLY_COHORT_LEX':[], 'DOWN_SAMPLE_TESTS':[], 'TRUNCATED':[]},
                'median':{'COHORT_LEX':[], 'PROG_ONLY_COHORT_LEX':[], 'DOWN_SAMPLE_TESTS':[], 'TRUNCATED':[]}, 
                'smallest':{'COHORT_LEX':[], 'PROG_ONLY_COHORT_LEX':[], 'DOWN_SAMPLE_TESTS':[], 'TRUNCATED':[]},
                'small-or-large':{'COHORT_LEX':[], 'PROG_ONLY_COHORT_LEX':[], 'DOWN_SAMPLE_TESTS':[], 'TRUNCATED':[]},
                'sum-of-squares':{'COHORT_LEX':[], 'PROG_ONLY_COHORT_LEX':[], 'DOWN_SAMPLE_TESTS':[], 'TRUNCATED':[]}}

#Will hold all of the original ranges for number of cohorts!
RANGE={}

#Will generate dictionary with boundaries!
#Expecting that they are 50 apart and powers of 2!
REPLICATES = 50
def create(n):
    #Will creat the min and max ranges between replicates!
    for i in range(n-1):
        minn = 1 + (REPLICATES * i)
        maxx = REPLICATES * (i+1)
        RANGE[2 ** i] = tuple([minn, maxx])

def main():
    parser = argparse.ArgumentParser(description="Data aggregation script.")
    parser.add_argument("data_directory", type=str, help="Target experiment directory.")
    parser.add_argument("dump_directory", type=str, help="Target dump directory")
    parser.add_argument("pairings", type=int, help="Number of CN,C# pairings there are")

    args = parser.parse_args()
    data_directory = args.data_directory
    write_dir = args.dump_directory
    create(args.pairings)
     
    if not os.path.exists(data_directory):
        print('Data directory does not exist!')
        return 0

#####Will start looking into the incomplete directories in the csv file!
    print('****Looking for all incomplete directories from undone.csv!**** \n')
    df = pd.read_csv(data_directory+'undone.csv')
    df = df.values.tolist()
    missing = {}
    for row in df:
        treat = row[0].split('__') 
        prob = treat[0][len('PROBLEM_'):]
        sel = treat[1][4:]
        cn = int(treat[2].strip('CN_'))
        seed = SEED_ORIGINAL[prob][sel]

        if prob not in missing:
            missing[prob] = {}
            if sel not in missing[prob]:
                missing[prob][sel] = [row[1] - seed]
            else:
                missing[prob][sel].append(row[1] - seed)
        else:
            if sel not in missing[prob]:
                missing[prob][sel] = [row[1] - seed]
            else:
                missing[prob][sel].append(row[1] - seed)

    for prob in missing.keys():
        for key,val in missing[prob].items():
            missing1 = str(val).strip('[]')
            final = ''
            for char in missing1:
                if char != ' ':
                    final += char
            
            print(prob+'__'+key+': ')
            print('Incomplete: ', final)
            print()
            TREATMENTS[prob][key] += val

#####Will start looking at all the missing directories within the missing text file!
    print('****Looking for all missing directories from undone.txt!**** \n')
    f = open(data_directory+'undone.txt')
    missing = {}
    for line in f:
        l = line.split(':')
        if l[0] == '== Treatment':
            treat = line.split(':')[1].strip(' =\n')
            minn = int(f.readline().strip().split(':')[1])
            maxx = int(f.readline().strip().split(':')[1])
            guess = f.readline().strip().split(':')[1]
            guess = guess.strip('[]\n ').split(',')

            if guess == ['']:
                missing[treat] = [minn, maxx, []]

            else:
                guess = [int(x) for x in guess]
                missing[treat] = [minn, maxx, guess]

    for k,vals in missing.items():
        treat = k.split('__')
        prob = treat[0][len('PROBLEM_'):]
        sel = treat[1][4:]
        cn = int(treat[2].strip('CN_'))
        minn = RANGE[cn][0]
        maxx = RANGE[cn][1]
        seed = SEED_ORIGINAL[prob][sel]

        print(prob + '__' + sel + '__CN_' + str(cn) + ': ')
        lower_b = vals[0] - SEED_ORIGINAL[prob][sel]
        lower_l = []
        if minn != lower_b:
            lower_l = list(range(minn,lower_b+1))
        upper_b = vals[1] - SEED_ORIGINAL[prob][sel]
        upper_l = []
        if maxx != upper_b:
            upper_l = list(range(upper_b, maxx+1))

        
        missing1 = lower_l + upper_l + [(x-seed) for x in vals[2]]
        miss = missing1
        final = ''
        for char in str(missing1).strip(' []'):
            if char != ' ':
                final += char

        TREATMENTS[prob][sel] += miss
        print('Missing: ' + final)
        print()


####Will write the file and output all the required ids that need rerun####
    fp = open(write_dir+'undone_list.txt', 'w')
    print()
    fp.write("Id's that need to be reran per problem and selection \n \n ")
    print('****Printing out the final list of all directories that need a rerun!****\n')
    for prob in TREATMENTS.keys():
        fp.write(prob+': \n \n')
        print(prob+': \n \n')
        final = ''
        for k,v in TREATMENTS[prob].items():
            v.sort()
            final = ''
            for char in str(v).strip('[]'):
                if char != ' ':
                    final += char

            fp.write(k + ': '+ final + '\n')
            print(k + ': ', final + '\n')

        fp.write('\n')
    

if __name__ == "__main__":
    main()