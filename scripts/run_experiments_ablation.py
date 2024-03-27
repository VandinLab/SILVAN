import os
import time
import numpy as np
import math
import pandas as pd

graph_experiments_paths = "graphs_experiments_abl.csv"
graphs = pd.read_csv(graph_experiments_paths,sep=";")
graphs = graphs.set_index('graph_name').T.to_dict('list')
print(graphs)
#graphs = graphs["graph_name"].values

num_runs = 10
db_path = "../datasets/"

epsilons = [0.01 , 0.005 , 0.0025 , 0.001 , 0.0005]
print("epsilons",epsilons)


for run_id in range(num_runs):
    for epsilon in epsilons:
        for graph_name in graphs:
            directed = graphs[graph_name]
            directed = directed[0] == 'D'
            directed_flag = ""
            if directed == True:
                directed_flag = " -t 1"
            # run the experiment
            graph_name = graph_name.replace(".txt","_pre.txt")
            cmd = "python3 run_silvan.py -db "+db_path+graph_name+" -mh 0 -eempp 0 -o results_silvan_abl.csv -e "+str(epsilon)+directed_flag
            print("run ",run_id," ",cmd)
            os.system(cmd)
            time.sleep(1)
            cmd = "python3 run_silvan.py -db "+db_path+graph_name+" -mh 0 -o results_silvan_abl.csv -e "+str(epsilon)+directed_flag
            print("run ",run_id," ",cmd)
            os.system(cmd)
            time.sleep(1)
            cmd = "python3 run_silvan.py -db "+db_path+graph_name+" -o results_silvan_abl.csv -e "+str(epsilon)+directed_flag
            print("run ",run_id," ",cmd)
            os.system(cmd)
            time.sleep(1)
