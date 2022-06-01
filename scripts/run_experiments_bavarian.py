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

num_runs = 5
db_path = "../datasets/"
approx_modes = ["rk" , "ab" , "bp"]

# parameters to choose epsilon
max_eps = 10**-2
min_eps = 10**-4
numeps_ = 5
expertoskip = set()

choose_random_eps = 0
if choose_random_eps == 1:
    epsilons = [1.0]
else:
    epsilons = list(np.logspace(math.log(min_eps,10) , math.log(max_eps,10),numeps_))
    epsilons.reverse()
    epsilons = [0.01 , 0.005 , 0.0025 , 0.001 , 0.0005]
    print("epsilons",epsilons)


def check_experiment(eps , db , appx):
    #print("checking ",eps," ",db," ",appx)
    df = pd.read_csv("results_bavarian.csv",sep=";")
    df = df[ df["graph_name"]==db ]
    #print(df)
    max_eps = max(epsilons)
    #print("max_eps",max_eps)
    df_ = df[ df["epsilon"]==max_eps ]
    #print(df)
    if df_.shape[0] == 0:
        return 0
    if not df_["time"].max() > 0.:
        return 1
    df_ = df[ df["estimator"]==appx ]
    if not df_["time"].min() > 0.:
        return 1

    df = df[ df["epsilon"]==eps ]
    df = df[ df["estimator"]==appx ]
    if df.shape[0] > 0:
        return 1
    return 0


for run_id in range(num_runs):
    for epsilon in epsilons:
        for graph_name in graphs:
            directed = graphs[graph_name]
            directed = directed[0] == 'D'
            directed_flag = ""
            if directed == True:
                directed_flag = " -t 1"
            # run the experiment
            graph_name = graph_name.replace(".txt","_bav.txt")
            for approx_ in approx_modes:
                already_run = check_experiment(epsilon , db_path+graph_name , approx_)
                maxtime_ = 21600
                cmd = "python3 run_bavarian.py -db "+db_path+graph_name+" -e "+str(epsilon)+directed_flag+" -a "+str(approx_)+" -maxtime "+str(maxtime_)
                if cmd not in expertoskip and already_run == 0:
                    print("run ",run_id," ",cmd)
                    start = time.time()
                    os.system(cmd)
                    end = time.time()
                    if end - start > 21000:
                        expertoskip.add(cmd)
                    time.sleep(1)
