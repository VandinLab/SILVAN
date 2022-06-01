import os
import time
import numpy as np
import math
import pandas as pd

graph_experiments_paths = "graphs_experiments.csv"
graphs = pd.read_csv(graph_experiments_paths,sep=";")
graphs = graphs.set_index('graph_name').T.to_dict('list')
print(graphs)
#graphs = graphs["graph_name"].values

num_runs = 10
db_path = "../datasets/"

# parameters to choose epsilon
max_eps = 10**-2
min_eps = 10**-4
numeps_ = 5
run_kadabra = 1

choose_random_eps = 0
if choose_random_eps == 1:
    epsilons = [1.0]
else:
    epsilons = list(np.logspace(math.log(min_eps,10) , math.log(max_eps,10),numeps_))
    epsilons.reverse()
    epsilons = [0.01 , 0.005 , 0.0025 , 0.001 , 0.0005]
    print("epsilons",epsilons)


for run_id in range(num_runs):
    for epsilon in epsilons:
        for graph_name in graphs:
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
            directed = graphs[graph_name]
            directed = directed[0] == 'D'
            directed_flag = ""
            if directed == True:
                directed_flag = " -t 1"
            # run the experiment
            graph_name = graph_name.replace(".txt","_pre.txt")
            cmd = "python3 run_kadabra.py -db "+db_path+graph_name+" -e "+str(epsilon)+directed_flag
            if run_kadabra == 1:
                print("run ",run_id," ",cmd)
                os.system(cmd)
                time.sleep(1)
            cmd = "python3 run_silvan.py -db "+db_path+graph_name+" -e "+str(epsilon)+directed_flag
            print("run ",run_id," ",cmd)
            os.system(cmd)
            time.sleep(1)
