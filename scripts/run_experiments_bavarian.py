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
approx_modes = ["rk" , "ab" , "bp"]

expertoskip = set()
epsilons = [0.01 , 0.005 , 0.0025 , 0.001 , 0.0005]
print("epsilons",epsilons)


def check_experiment(eps , db , appx):
    if os.path.isfile("results_bavarian.csv") == False:
        return 0
    #print("checking ",eps," ",db," ",appx)
    df_bav_res = pd.read_csv("results_bavarian.csv",sep=";")

    #check if this run was already done
    df_g_e = df_bav_res.copy()
    df_g_e = df_g_e[ df_g_e["graph_name"]==db ]
    df_g_e = df_g_e[ df_g_e["estimator"]==appx ]

    max_eps = max(epsilons)
    if max_eps == eps:
        df_max_eps = df_g_e[ df_g_e["epsilon"]==max_eps ]
        # consider the larget epsilon first
        # if no experiments was run with the largest epsilon, run it
        if df_max_eps.shape[0] == 0:
            return 0

    # check if current experiment will not terminate
    df_g_e_term = df_g_e[ df_g_e["epsilon"]>=eps ]
    if df_g_e_term.shape[0] > 0:
        if not df_g_e_term["terminated"].min() > 0:
            # this experiment will not terminate
            print("skipping exp (",eps,db,appx,") for reason 1")
            return 1

    # consider the current epsilon
    df_curr_eps = df_g_e[ df_g_e["epsilon"]==eps ]
    # do not repeat if number of runs is enough
    if df_curr_eps.shape[0] >= num_runs:
        print("skipping exp (",eps,db,appx,") for reason 2")
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
