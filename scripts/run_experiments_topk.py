import os
import time
import numpy as np
import math

bcest_max = {
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

num_runs = 10
epsilons = [0.25 , 0.1 , 0.05]
ks = [5 , 10 , 25]
db_path = "../datasets/"

print("rel epsilons",epsilons)

for run_id in range(num_runs):
    for epsilon in epsilons:
        for k in ks:
            for index , dataset in enumerate(datasets):
                cmd = "python3 run_kadabra_topk.py -k "+str(k)+" -db "+db_path+dataset+" -e "+str(epsilon)
                print("run ",run_id," ",cmd)
                os.system(cmd)
                time.sleep(1)
