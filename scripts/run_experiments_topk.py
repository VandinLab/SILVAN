import os
import time
import numpy as np
import math
import pandas as pd
from datetime import datetime

graph_experiments_paths = "graphs_experiments.csv"
graphs = pd.read_csv(graph_experiments_paths,sep=";")
graphs = graphs.set_index('graph_name').T.to_dict('list')
print(graphs)

num_runs = 10
db_path = "../datasets/"
run_silvan = 1

# parameters for topk experiments
ks = [5 , 10 , 25]
epsrels = [0.25 , 0.1 , 0.05]

if run_silvan == 1:
    for run_id in range(num_runs):
        for k in ks:
            for epsrel in epsrels:
                for graph_name in graphs:
                    directed = graphs[graph_name]
                    directed = directed[0] == 'D'
                    directed_flag = ""
                    if directed == True:
                        directed_flag = " -t 1"
                    # run the experiment
                    graph_name = graph_name.replace(".txt","_pre.txt")
                    epsilon = epsrel
                    cmd = "python3 run_silvan.py -k "+str(k)+" -db "+db_path+graph_name+" -e "+str(epsilon)+directed_flag
                    now = datetime.now()
                    print("now =", now)
                    print("run ",run_id," ",cmd)
                    os.system(cmd)
                    time.sleep(1)

# when finished, we have the parameters to use for kadabra
df = pd.read_csv("results_silvan_topk.csv",sep=";")
for index, row in df.iterrows():
    k = row['k']
    epsrel_ = row['epsilon']
    graph_type = row['type']
    try:
        k = int(k)
        epsrel_ = float(epsrel_)
        graph_type = int(graph_type)
    except ValueError:
        pass
    epsilon_kadabra = row['eps_final_topk']
    time_silvan = row['time']
    max_time_kadabra = 43200
    try:
        max_time_kadabra = min(max_time_kadabra,100*float(time_silvan))
    except ValueError:
        pass
    graph_path = row['graph_name']
    directed_flag = ""
    directed = graph_type == 1
    if directed == True:
        directed_flag = " -t 1"
    # first run kadabra using the -k parameter
    cmd = "python3 run_kadabra.py -k "+str(k)+" -db "+graph_path+" -e "+str(epsilon_kadabra)+directed_flag
    cmd = cmd+" -maxtime "+str(int(max_time_kadabra))
    #print("run ",run_id," ",cmd)
    now = datetime.now()
    print("now =", now)
    print(cmd)
    os.system(cmd)
    time.sleep(1)
