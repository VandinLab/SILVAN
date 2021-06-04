import os
import math
import time

in_mcrade_exp = open("results_topk_new.csv","r")
path_to_datasets = "../datasets/"
path_to_kadabra = "../kadabra/kadabra"

min_time_kadabra = 5*60.



def takeresultskadabra():

    kadabra_runout = 1
    kadabra_samples = 0
    kadabra_time = 0.
    file_out = open("out_k.txt","r")
    for line in file_out:
        if "Finished after " in line:
            #this means kadabra finished in reasonable time
            kadabra_runout = 0
            # take number of samples
            line_ = line.replace("\n","")
            line_ = line_.replace("Finished after ","")
            line_ = line_.replace(" iterations.","")
            kadabra_samples = int(line_)
        if kadabra_runout == 0:
            if "Total time: " in line:
                line_ = line.replace("\n","")
                kadabra_time = line_.replace("Total time: ","")
                kadabra_time = float(kadabra_time)
    file_out.close()
    if kadabra_runout == 1:
        # kadabra did not finish in reasonable time
        # take max num of samples and max time
        str_iterations = "Situation after "
        file_out = open("out_k.txt","r")
        for line in file_out:
            if str_iterations in line:
                line_ = line.replace("\n","")
                line_ = line_.replace(str_iterations,"")
                line_ = line_.replace(" iterations.","")
                kadabra_samples = max(kadabra_samples , int(line_))

    return kadabra_samples , kadabra_runout , kadabra_time


def run_kadabra(db , eps , mcradetime , k_exp , strict):

    db = db.replace("../datasets/","")

    print("run kadabra experiment: ",db,eps,mcradetime)

    mcradetime = float(mcradetime)
    kadabra_max_time = int(math.ceil(mcradetime*10. + 1.))
    kadabra_max_time = max(kadabra_max_time,min_time_kadabra)
    timeout_str = "timeout "+str(kadabra_max_time)+"s "
    if strict == 0:
        cmd = timeout_str+path_to_kadabra+" -k "+str(k_exp)+" "+str(eps)+" 0.05 "+path_to_datasets+db+" > out_k.txt 2>&1"
    else:
        cmd = timeout_str+path_to_kadabra+" "+str(eps)+" 0.05 "+path_to_datasets+db+" > out_k.txt 2>&1"
    print(cmd)
    os.system(cmd)
    time.sleep(1)

    kadabra_samples , kadabra_runout , kadabra_time = takeresultskadabra()
    if kadabra_runout == 1:
        kadabra_time = kadabra_max_time

    time.sleep(1)

    return kadabra_samples , kadabra_runout , kadabra_time


diameters={ '../datasets/dip20090126_MAX.nde': 6,
'../datasets/in_2004.nde': 8,
 '../datasets/p2p-Gnutella31.txt': 16,
 '../datasets/cnr_2000.nde': 8,
 '../datasets/com-amazon.ungraph.txt': 56,
 '../datasets/email-Enron.txt': 16,
 '../datasets/ca-GrQc.txt' : 20,
 '../datasets/oregon1_010526.txt': 12,
 '../datasets/soc-Epinions1.txt' : 16,
 '../datasets/wiki-Vote.txt': 8}


def compute_eps_kadabra(delta , m , diam , epsrel):

    eps_omega = math.sqrt((math.log(2*(diam-1),2)+math.log(1./delta))/(2*m) )
    # correct eps to get same guarantees of SILVAN
    eps_omega = eps_omega*epsrel*((1-epsrel)/(1+epsrel))**2
    return eps_omega


for line in in_mcrade_exp:
    line_ = line.replace("\n","")
    items = line_.split(";")
    db_exp = items[0]
    k_exp = items[1]
    epsrel_exp = float(items[2])
    delta = float(items[3])
    mcradesamples_exp = items[5]
    mcradetime_exp = items[6]
    omega_eps = float(items[7])
    eps_kadabra = compute_eps_kadabra(delta , omega_eps , diameters[db_exp] , epsrel_exp)
    # experiment with top-k guarantees from KADABRA
    numsamples_kadabra , runout_kadabra , kadkadabra_time = run_kadabra(db_exp , eps_kadabra , mcradetime_exp , k_exp , 0)
    out_kadabra_exp = open("kadabra_experiments_out.csv","a")
    out_kadabra_exp.write(line_+str(numsamples_kadabra)+";"+str(runout_kadabra)+";"+str(kadkadabra_time)+"\n")
    out_kadabra_exp.close()

    # experiment with no relaxation of guarantees
    numsamples_kadabra , runout_kadabra , kadkadabra_time = run_kadabra(db_exp , eps_kadabra , mcradetime_exp , k_exp , 1)
    out_kadabra_exp = open("kadabra_strict_experiments_out.csv","a")
    out_kadabra_exp.write(line_+str(numsamples_kadabra)+";"+str(runout_kadabra)+";"+str(kadkadabra_time)+"\n")
    out_kadabra_exp.close()
