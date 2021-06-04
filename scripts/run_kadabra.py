import os
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-db", help="path to graph input file")
parser.add_argument("-e","--epsilon", type=float ,help="approximation accuracy parameter (in (0,1))",default=0.01)
parser.add_argument("-d","--delta", type=float ,help="approximation confidence parameter (in (0,1), def. 0.1)",default=0.05)
args = parser.parse_args()

delta = args.delta
kadabra_output_path = "output.txt"
cmd = "../mcrade-kadabra/kadabra "+str(args.epsilon)+" "+str(delta)+" "+str(args.db)+" > "+kadabra_output_path
print(cmd)
os.system(cmd)

output_path = "results.csv"

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
    if verbose == 1:
        print(line_to_print)
    return to_return

# results to gather from output file
results_strings = ["Finished after ","Total time: ", "sup_mcera_all: ","sup bc est: ","estimated diameter of the graph: ","Maximum confidence interval: ",
"sup_bcest_high_freq: ","sup_mcera_high_freq: ","sup_bcest_low_freq: ","sup_mcera_low_freq: ","high_freq_elements: ","MCRADE STOPS WITH eps ","MCRADE STOPS at iteration ",
"MCRADE STOPS at sample size ","MCRADE STOPS after seconds "]

out_file = open(output_path,"a")
line_out = str(args.db)+";"
for result_string in results_strings:
    value_result = get_result(result_string , kadabra_output_path ,  1)
    if " iterations." in value_result:
        value_result = value_result.replace(" iterations.","")
    print(value_result)
    print(float(value_result))
    line_out = line_out+value_result+";"
line_out = line_out+"\n"
print(line_out)
out_file.write(line_out)
