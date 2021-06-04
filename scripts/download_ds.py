import os

dblinks = [
"https://www.pilucrescenzi.it/lasagne/graphs_dataset/nde/web/cnr_2000.zip",
"https://www.pilucrescenzi.it/lasagne/graphs_dataset/nde/web/in_2004.zip",
"http://snap.stanford.edu/data/p2p-Gnutella31.txt.gz",
"http://snap.stanford.edu/data/bigdata/communities/com-amazon.ungraph.txt.gz",
"http://snap.stanford.edu/data/cit-HepTh.txt.gz",
"https://www.pilucrescenzi.it/lasagne/graphs_dataset/nde/biological/dip20090126_MAX.zip",
"http://snap.stanford.edu/data/email-Enron.txt.gz",
"http://snap.stanford.edu/data/ca-GrQc.txt.gz",
"http://snap.stanford.edu/data/oregon1_010526.txt.gz",
"http://snap.stanford.edu/data/soc-Epinions1.txt.gz",
"http://snap.stanford.edu/data/wiki-Vote.txt.gz",
"http://snap.stanford.edu/data/cit-HepPh.txt.gz"
]

filepaths = [
"cnr_2000.zip",
"in_2004.zip",
"p2p-Gnutella31.txt.gz",
"com-amazon.ungraph.txt.gz",
"cit-HepTh.txt.gz",
"dip20090126_MAX.zip",
"email-Enron.txt.gz",
"ca-GrQc.txt.gz",
"oregon1_010526.txt.gz",
"soc-Epinions1.txt.gz",
"wiki-Vote.txt.gz",
"cit-HepPh.txt.gz"
]


os.system("cd .. && mkdir datasets")

for link in dblinks:
    os.system("cd ../datasets && wget -N "+link)

for filepath in filepaths:
    if "zip" in filepath:
        os.system("cd ../datasets && unzip -o "+filepath)
    if ".gz" in filepath:
        os.system("cd ../datasets && gunzip -f "+filepath)
