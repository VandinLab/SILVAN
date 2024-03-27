import os
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-db", help="path to graph input file")
parser.add_argument("-e","--epsilon", type=float ,help="approximation accuracy parameter (in (0,1))",default=0.01)
parser.add_argument("-d","--delta", type=float ,help="approximation confidence parameter (in (0,1), def. 0.1)",default=0.05)
parser.add_argument("-k", type=int ,help="parameter for top-k approximation",default=0)
parser.add_argument("-t","--type", type=int ,help="type of graph. 1 means directed, 0 undirected (def. undirected)",default=0)
parser.add_argument("-maxtime", type=int ,help="maximum seconds before killing kadabra",default=-1)
args = parser.parse_args()

delta = args.delta
file_output_path = "output.txt"
path_executable = "../kadabra/kadabra"
if args.maxtime > 0:
    path_executable = "timeout "+str(args.maxtime)+" "+path_executable
output_path = "results_kadabra.csv"
directed_flag = ""
if args.type == 1:
    directed_flag = "-d "
topk_flag = ""
if args.k > 0:
    topk_flag = "-k "+str(args.k)+" "
    output_path = "results_kadabra_topk.csv"
    file_output_path = "output_topk.txt"
cmd = path_executable+" "+directed_flag+topk_flag+str(args.epsilon)+" "+str(delta)+" "+str(args.db)+" > "+file_output_path
print(cmd)
retval = os.system(cmd)
if retval > 0:
    print("KADABRA terminated with errors")
    exit()


if os.path.isfile(output_path) == False:
    # write the header for the results
    out = "graph_name;epsilon;delta;k;type;m;time;diameter;first_samples;first_time;eps_mcera;iter_index_mcera;top1_est;top1_ub;avg_spl_est;avg_spl_ub;samples_bound;time_bfs;time_critical;eps_final_topk\n"
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
#results_strings = ["Finished after ","Total time: ", "sup_mcera_all: ","sup bc est: ","estimated diameter of the graph: ","Maximum confidence interval: ",
#"sup_bcest_high_freq: ","sup_mcera_high_freq: ","sup_bcest_low_freq: ","sup_mcera_low_freq: ","high_freq_elements: ","MCRADE STOPS WITH eps ","MCRADE STOPS at iteration ",
results_strings = ["Finished after ","Total time: ","estimated diameter of the graph: ","First pass finished after ", "time for first pass ", "MCRADE STOPS WITH eps ","MCRADE STOPS at iteration ","top-1 bc first pass: ","top1bc_upperbound: ","avg_diam_firstpass: ", "avg_diam_upperbound: ", "max_num_samples: ","Time bfs: ","Time critical: ","eps_final_topk "]

out_file = open(output_path,"a")
line_out = str(args.db)+";"+str(args.epsilon)+";"+str(args.delta)+";"+str(args.k)+";"+str(args.type)+";"
for idx, result_string in enumerate(results_strings):
    value_result = get_result(result_string , file_output_path ,  1)
    if " iterations." in value_result:
        value_result = value_result.replace(" iterations.","")
    if " iteration" in value_result:
        value_result = value_result.replace(" iteration","")
    #print(float(value_result))
    line_out = line_out+value_result
    if idx < len(results_strings)-1:
        line_out = line_out+";"
line_out = line_out+"\n"
print(line_out)
out_file.write(line_out)
