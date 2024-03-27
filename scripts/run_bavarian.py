
import os
import time
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-db", help="path to graph input file")
parser.add_argument("-e","--epsilon", type=float ,help="approximation accuracy parameter (in (0,1))",default=0.01)
parser.add_argument("-d","--delta", type=float ,help="approximation confidence parameter (in (0,1), def. 0.1)",default=0.05)
parser.add_argument("-t","--type", type=int ,help="type of graph. 1 means directed, 0 undirected (def. undirected)",default=0)
parser.add_argument("-a", "--approxmode", help="type of estimator (rk, ab, or bp)",default="rk")
parser.add_argument("-mc", "--mctrials", type=int, help="number of mc trials",default=25)
parser.add_argument("-maxtime", type=int ,help="maximum seconds before killing bavarian",default=-1)
args = parser.parse_args()

delta = args.delta
start_t = time.time()
file_output_path = "output_bavarian.txt"
path_executable = "../bavarian/bin/progrbavarian"
if args.maxtime > 0:
    path_executable = "timeout "+str(args.maxtime)+" "+path_executable
output_path = "results_bavarian.csv"
directed_flag = ""
if args.type == 1:
    directed_flag = "-d "
cmd = path_executable+" "+directed_flag+"-m 1.2 -v -v 1 "+str(args.epsilon)+" "+str(args.mctrials)+" "+str(args.approxmode)+" "+str(delta)+" "+str(args.db)+" > "+file_output_path+" 2>&1"
print(cmd)
retval = os.system(cmd)
end_t = time.time()
true_time_bavarian = end_t - start_t
if retval > 0:
    print("Bavarian terminated with errors")
    exit()

if os.path.isfile(output_path) == False:
    # write the header for the results
    out = "graph_name;epsilon;delta;type;m;mctrials;time;true_time;estimator;terminated\n"
    out_file = open(output_path,"w")
    out_file.write(out)
    print(out)
    out_file.close()

# get running time and number of samples
time_bavarian = 0.
samples_bavarian = 0.
bavarian_terminated = 0
fin = open(file_output_path,"r")
for line in fin:
    if "[INFO ]: Iteration " in line:
        items_ = line.split(" ")
        for idx,item_ in enumerate(items_):
            if "taking" in item_:
                samples_bavarian = max(samples_bavarian , int(items_[idx+1]))
    if line[0:4] == "done":
        line = line.replace("\n","")
        line = line.replace(" millisecs)","")
        line = line.replace("done (","")
        time_bavarian = float(line)/1000.
    if "INFO: Getting results...done" in line:
        bavarian_terminated = 1
if bavarian_terminated == 0:
    time_bavarian = 0.
    samples_bavarian = 0.
print("bavarian_terminated",bavarian_terminated,"time_bavarian",time_bavarian,"true_time_bavarian",true_time_bavarian,"samples_bavarian",samples_bavarian)

out_file = open(output_path,"a")
line_out = str(args.db)+";"+str(args.epsilon)+";"+str(args.delta)+";"+str(args.type)+";"+str(samples_bavarian)+";"+str(args.mctrials)+";"+str(time_bavarian)+";"+str(true_time_bavarian)+";"+str(args.approxmode)+";"+str(bavarian_terminated)
line_out = line_out+"\n"
print(line_out)
out_file.write(line_out)
