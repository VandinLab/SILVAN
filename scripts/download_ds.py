import os

directed_graphs_links = [
"http://snap.stanford.edu/data/p2p-Gnutella31.txt.gz",
"http://snap.stanford.edu/data/cit-HepTh.txt.gz",
"http://snap.stanford.edu/data/soc-Epinions1.txt.gz",
"http://snap.stanford.edu/data/wiki-Vote.txt.gz",
"http://snap.stanford.edu/data/cit-HepPh.txt.gz",
"http://snap.stanford.edu/data/wiki-topcats.txt.gz",
"http://konect.cc/files/download.tsv.wikipedia_link_en.tar.bz2", # http://konect.cc/networks/wikipedia_link_en/
"http://snap.stanford.edu/data/email-EuAll.txt.gz", # http://snap.stanford.edu/data/email-EuAll.html
"http://snap.stanford.edu/data/wiki-Talk.txt.gz", # http://snap.stanford.edu/data/wiki-Talk.html
"http://snap.stanford.edu/data/soc-LiveJournal1.txt.gz", # http://snap.stanford.edu/data/soc-LiveJournal1.html
"http://snap.stanford.edu/data/soc-pokec-relationships.txt.gz" # http://snap.stanford.edu/data/soc-Pokec.html
]

undirected_graphs_links = [
"http://snap.stanford.edu/data/bigdata/communities/com-amazon.ungraph.txt.gz",
"http://snap.stanford.edu/data/email-Enron.txt.gz",
"http://snap.stanford.edu/data/ca-GrQc.txt.gz",
"http://snap.stanford.edu/data/bigdata/communities/com-youtube.ungraph.txt.gz",
"http://konect.cc/files/download.tsv.actor-collaboration.tar.bz2", # http://konect.cc/networks/actor-collaboration/
"http://snap.stanford.edu/data/bigdata/communities/com-dblp.ungraph.txt.gz", # http://snap.stanford.edu/data/com-DBLP.html
"http://snap.stanford.edu/data/ca-AstroPh.txt.gz" # http://snap.stanford.edu/data/ca-AstroPh.html
]

all_links = directed_graphs_links+undirected_graphs_links

filepaths = []
graph_experiments_paths = "graphs_experiments.csv"
graphs_list_file = open(graph_experiments_paths,"w")
graphs_list_file.write("graph_name;directed\n")
for link in all_links:
    filepath = link.rfind("/")
    filepath = link[filepath+1:]
    print(filepath)
    filepaths.append(filepath)
    filepath = filepath.replace(".gz","")
    filepath = filepath.replace(".zip","")
    filepath = filepath.replace(".tar","")
    filepath = filepath.replace(".bz2","")
    filepath = filepath.replace("download.tsv.","")
    if ".txt" not in filepath:
        filepath = filepath+".txt"
    if link in undirected_graphs_links:
        graphs_list_file.write(filepath+";U\n")
    else:
        graphs_list_file.write(filepath+";D\n")
graphs_list_file.close()

os.system("cd .. && mkdir datasets")

for link in all_links:
    os.system("cd ../datasets && wget -N "+link)

for filepath in filepaths:
    if "zip" in filepath:
        os.system("cd ../datasets && unzip -o "+filepath)
    if ".gz" in filepath:
        os.system("cd ../datasets && gunzip -kf "+filepath)
    if ".tar.bz2" in filepath:
        os.system("cd ../datasets && tar xjf "+filepath)
    if "download.tsv." in filepath:
        filepath_ = filepath.replace("download.tsv.","")
        filepath_ = filepath_.replace(".tar.bz2","")
        file_to_copy_path = "../datasets/"+filepath_+"/out."+filepath_
        file_final_path = "../datasets/"+filepath_+".txt"
        os.system("cp "+file_to_copy_path+" "+file_final_path)

import pandas as pd

graph_experiments_paths = "graphs_experiments.csv"
graphs = pd.read_csv(graph_experiments_paths,sep=";")
graphs = graphs["graph_name"].values
for graph_filename in graphs:
    if os.path.isfile("../datasets/"+graph_filename) == False:
        print("Graph file for ",graph_filename," not found!")
