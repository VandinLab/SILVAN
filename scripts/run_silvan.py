import os
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-db", help="path to graph input file")
parser.add_argument("-e","--epsilon", type=float ,help="approximation accuracy parameter (in (0,1))",default=0.01)
parser.add_argument("-d","--delta", type=float ,help="approximation confidence parameter (in (0,1), def. 0.1)",default=0.05)
parser.add_argument("-a","--aemppeeling", type=float ,help="parameter a for empirical peeling (def. = 2)",default=2)
parser.add_argument("-s","--alpha", type=float ,help="parameter alpha for sampling shortest paths (def. = 2.3)",default=2.3)
parser.add_argument("-k", type=int ,help="parameter for top-k approximation",default=0)
parser.add_argument("-mh", type=int ,help="enable computation of m_hat (def.=1)",default=1)
parser.add_argument("-eempp", type=int ,help="enable empirical peeling (def.=1)",default=1)
parser.add_argument("-o", type=str, help="output path (def.=results_silvan.csv or results_silvan_topk.csv)",default="results_silvan.csv")
parser.add_argument("-t","--type", type=int ,help="type of graph. 1 means directed, 0 undirected (def. undirected)",default=0)
args = parser.parse_args()

delta = args.delta
file_output_path = "output.txt"
path_executable = "../silvan/silvan"
output_path = args.o
directed_flag = ""
if args.type == 1:
    directed_flag = "-d "
topk_flag = ""
if args.epsilon <= 0.:
    print("the parameter epsilon should be positive!")
    exit()
algname = "SILVAN"
if args.k > 0:
    topk_flag = "-k "+str(args.k)+" "
    output_path = "results_silvan_topk.csv"
    file_output_path = "output_topk.txt"
    algname = "SILVAN-TOPK"
optional_parameters = ""
if args.alpha != 2.3:
    if args.alpha > 1:
        optional_parameters = "-s "+str(args.alpha)+" "
    else:
        print("The input parameter args.alpha should be > 1!")
if args.aemppeeling != 2:
    if args.aemppeeling > 1 and args.eempp == 1:
        optional_parameters = optional_parameters+"-a "+str(args.aemppeeling)+" "
    else:
        print("The input parameter args.aemppeeling should be > 1!")
if args.eempp == 0:
    optional_parameters = optional_parameters+"-a 0 "
if args.mh == 0:
    optional_parameters = optional_parameters+"-m "

print("Running",algname,"with eps =",args.epsilon)
cmd = path_executable+" "+directed_flag+topk_flag+optional_parameters+str(args.epsilon)+" "+str(delta)+" "+str(args.db)+" > "+file_output_path
print(cmd)
retval = os.system(cmd)
if retval > 0:
    print(algname,"terminated with errors")
    exit()



if os.path.isfile(output_path) == False:
    # write the header for the results
    out = "graph_name;epsilon;delta;k;type;m_hat_enabled;emp_peel_enabled;emp_peel_a;alpha;m;time;diameter;first_samples;first_time;eps_mcera;iter_index_mcera;top1_est;top1_ub;avg_spl_est;avg_spl_ub;samples_bound;time_bfs;time_critical;eps_final_topk\n"
    out_file = open(output_path,"w")
    out_file.write(out)
    print(out)
    out_file.close()

def get_result(pattern , path ,  verbose=1):
    fin = open(path,'r')
    to_return = ""
    line_to_print = ""
    for line in fin:
        if pattern in line:
            line = line.replace('\n','')
            line_to_print = line
            to_return = line[len(pattern):]
    fin.close()
    if verbose == 1 and len(line_to_print) > 0:
        print(line_to_print)
    return to_return

# results to gather from output file
results_strings = ["Finished after ","Total time: ","estimated diameter of the graph: ","First pass finished after ", "time for first pass ", "MCRADE STOPS WITH eps ","MCRADE STOPS at iteration ","top-1 bc first pass: ","    top1bc_upperbound: ","avg_diam_firstpass: ", "avg_diam_upperbound: ", "max_num_samples: ","Time bfs: ","Time critical: ","eps_final_topk "]

out_file = open(output_path,"a")
line_out = str(args.db)+";"+str(args.epsilon)+";"+str(args.delta)+";"+str(args.k)+";"+str(args.type)+";"+str(args.mh)+";"+str(args.eempp)+";"+str(args.aemppeeling)+";"+str(args.alpha)+";"
for idx, result_string in enumerate(results_strings):
    value_result = get_result(result_string , file_output_path ,  1)
    if " iterations." in value_result:
        value_result = value_result.replace(" iterations.","")
    if " iteration" in value_result:
        value_result = value_result.replace(" iteration","")
    #print(value_result)
    #print(float(value_result))
    line_out = line_out+value_result
    if idx < len(results_strings)-1:
        line_out = line_out+";"
line_out = line_out+"\n"
print(line_out)
out_file.write(line_out)
