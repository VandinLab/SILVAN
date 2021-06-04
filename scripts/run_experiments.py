import os
import time
import numpy as np
import math

bcest_max = {
    "cit-HepTh.txt" : 0.00001 ,
    "cit-HepPh.txt" : 0.00001 ,
    "dip20090126_MAX.nde" : 0.495560 ,
    "in_2004.nde" : 0.177073  ,
    "p2p-Gnutella31.txt" : 0.007749,
    "cnr_2000.nde" : 0.258061 ,
    "com-amazon.ungraph.txt" : 0.011748,
    "email-Enron.txt" : 0.065,
    "ca-GrQc.txt" : 0.00148,
    "oregon1_010526.txt" : 0.0097,
    "soc-Epinions1.txt" : 0.072,
    "wiki-Vote.txt" : 0.045
}
datasets = bcest_max.keys()

num_runs = 20
db_path = "../datasets/"

# parameters to choose epsilon
max_eps = 10**-2
min_eps = 10**-6
numeps_ = 7

choose_random_eps = 1
if choose_random_eps == 1:
    epsilons = [1.0]
else:
    epsilons = np.logspace(math.log(min_eps,10) , math.log(max_eps,10),numeps_)
    print("epsilons",epsilons)


for epsilon in epsilons:
    for run_id in range(num_runs):
        for index , dataset in enumerate(datasets):
            if choose_random_eps == 1:
                bcest_loc = bcest_max[dataset]
                # set minimum value of eps to choose
                min_eps_est = max(10**-3*bcest_loc , min_eps)
                min_eps_est = min(min_eps_est , 10**-4)
                # set maximum value of eps to choose
                max_eps_loc = max_eps
                # choose epsilon
                epsilon_rnd = np.exp(np.random.uniform(math.log(min_eps_est), math.log(max_eps_loc)))
                print(epsilon_rnd)
                epsilon = epsilon_rnd
            # run the experiment
            cmd = "python3 run_kadabra.py -db "+db_path+dataset+" -e "+str(epsilon)
            print("run ",run_id," ",cmd)
            os.system(cmd)
            time.sleep(1)
