# preprocess graphs by shifting node ids to consecutive values
# this is to make sure KADABRA and SILVAN work correctly
# convert space separated edges to tab separated for BAVARIAN

import os
import os.path
import pandas as pd

graph_experiments_paths = "graphs_experiments.csv"
graphs = pd.read_csv(graph_experiments_paths,sep=";")

graphs = graphs["graph_name"].values

for graph in graphs:

    graph_path = "../datasets/"+graph
    graph_path_out = graph_path
    graph_path_out = graph_path_out.replace(".txt","")
    graph_path_out = graph_path_out.replace(".nde","")
    graph_path_out_bav = graph_path_out+"_bav.txt"
    graph_path_out = graph_path_out+"_pre.txt"
    if not os.path.isfile(graph_path_out) or not os.path.isfile(graph_path_out_bav):
        fin = open(graph_path,"r")
        print("reading ",graph_path,"...")
        max_node_id = 0
        nodes_set = set()
        for line in fin:
            if "#" not in line and "%" not in line:
                line = line.replace(" \n","")
                line = line.replace("\n","")
                if " " in line or "\t" in line:
                    if " " in line:
                        nodes = line.split(" ")
                    else:
                        if "\t" in line:
                            nodes = line.split("\t")
                    try:
                        nodes = [int(n) for n in nodes]
                    except ValueError:
                        print("ValueError encountered!")
                        print(nodes)
                        print(line)
                        exit()
                    if len(nodes) >= 2:
                        for n in nodes:
                            nodes_set.add(n)
                            max_node_id = max(max_node_id , n)
        fin.close()

        node_ids_map = dict()

        new_node_id = 0
        for node_id in nodes_set:
            node_ids_map[node_id] = new_node_id
            new_node_id += 1

        fin = open(graph_path,"r")
        print("writing ",graph_path_out,"...")
        fout = open(graph_path_out,"w")
        fout_bav = open(graph_path_out_bav,"w")
        for line in fin:
            if "#" not in line and "%" not in line:
                line = line.replace(" \n","")
                line = line.replace("\n","")
                if " " in line or "\t" in line:
                    if " " in line:
                        nodes = line.split(" ")
                    else:
                        if "\t" in line:
                            nodes = line.split("\t")
                    try:
                        nodes = [int(n) for n in nodes]
                    except ValueError:
                        print("ValueError encountered!")
                        print(nodes)
                        print(line)
                        os.system("rm "+graph_path_out)
                        exit()
                    if len(nodes) >= 2:
                        new_edge_line = str(node_ids_map[nodes[0]])+" "+str(node_ids_map[nodes[1]])+"\n"
                        new_edge_line_bav = str(node_ids_map[nodes[0]])+"\t"+str(node_ids_map[nodes[1]])+"\n"
                        fout.write(new_edge_line)
                        fout_bav.write(new_edge_line_bav)
        fout.close()
        fout_bav.close()


        print("done graph",graph)
        print(max_node_id)
        print(len(nodes_set))
        print("")
    else:
        print("skipping graph",graph)
