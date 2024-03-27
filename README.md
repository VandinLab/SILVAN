# SILVAN: Estimating Betweenness Centralities with Progressive Sampling and Non-uniform Rademacher Bounds #

This repository contains the implementation of SILVAN, supporting our paper ["SILVAN: Estimating Betweenness Centralities with Progressive Sampling and Non-uniform Rademacher Bounds"](https://dl.acm.org/doi/10.1145/3628601).

SILVAN is an algorithm to compute approximations of the Betweenness Centralities from graphs with progressive sampling. More details can be found in the paper, published in ACM Transactions on Knowledge Discovery from Data: https://dl.acm.org/doi/10.1145/3628601

Part of the underlying implementation of SILVAN is based on the sampling algorithm [KADABRA](https://github.com/natema/kadabra) by Michele Borassi and Emanuele Natale. Therefore, it has the same compilation dependencies (described below) and it is distributed with the same license (Apache License 2.0).

## Installation

The software requires the [OpenMP API](http://openmp.org/wp/). After cloning this repository,
build the software by issuing the `make` command inside the silvan folder.

### Running SILVAN ###

The network has to be provided as a file containing two space-separated
integers `u v` per line, for each edge `(u,v)` of the network. The labels of
the vertices are assumed to be consecutive.

To run SILVAN, you can use the `run_silvan.py` python script found in the `scripts` folder. It takes the following input parameters:

```
usage: run_silvan.py [-h] [-db DB] [-e EPSILON] [-d DELTA] [-a AEMPPEELING]
                     [-s ALPHA] [-k K] [-mh MH] [-eempp EEMPP] [-o O]
                     [-t TYPE]
optional arguments:
  -h, --help            show this help message and exit
  -db DB                path to graph input file
  -e EPSILON, --epsilon EPSILON
                        approximation accuracy parameter (in (0,1))
  -d DELTA, --delta DELTA
                        approximation confidence parameter (in (0,1), def. 0.1)
  -a AEMPPEELING, --aemppeeling AEMPPEELING
                        parameter a for empirical peeling (def. = 2)
  -s ALPHA, --alpha ALPHA
                        parameter alpha for sampling shortest paths (def. = 2.3)
  -k K                  parameter for top-k approximation
  -mh MH                enable computation of m_hat (def.=1)
  -eempp EEMPP          enable empirical peeling (def.=1)
  -o O                  output path (def.=results_silvan.csv or results_silvan_topk.csv)
  -t TYPE, --type TYPE  type of graph. 1 means directed, 0 undirected (def. undirected)
```

For example, to approximate the betweenness centrality of all nodes of the undirected graph `graph.txt` with absolute accuracy 0.01 and with probability at least 0.95, you can use

`python run_silvan.py -db graph.txt -e 0.01 -d 0.05`

or to approximate the top-10 most central nodes of the directed graph `digraph.txt` with relative accuracy 0.1, you can use

`python run_silvan.py -db digraph.txt -e 0.1 -k 10 -t 1`

### Reproducing experiments
To reproduce the experiments described in the paper, follow the instructions listed below.

1. First you need to compile the algorithms. You can do so with the command `make` within the folders `kadabra` (found inside the `kadabra.zip` archive) and `silvan`. The folder kadabra contains a copy of the code from https://github.com/natema/kadabra with minimal modifications (mainly to better parse some of its statistics).
2. Then, move to the `script` folder. Download all graphs using `python download_ds.py`. Graphs will be downloaded in the `datasets` folder. Run `python preprocessing.py` to preprocess the downloaded graphs. Then, use the script `python setup_bavarian.py` to automatically setup BAVARIAN.
3. To run experiments described in Section 5.1, run `python run_experiments.py`. Results will be appended in the files `results_silvan.csv` for SILVAN, and `results_kadabra.csv` for KADABRA. The ablation experiments of Section 5.1.4 can be ran with the command `python run_experiments_ablation.py`
4. To run experiments described in Section 5.2, run `python run_experiments_topk.py`. Results will be appended in the file `results_silvan_topk.csv` for SILVAN, and in `results_kadabra_topk.csv` for KADABRA.
5. To run all experiments for BAVARIAN, run `python run_experiments_bavarian.py`. Results will be appended in the file `results_bavarian.csv`.

### Contacts ###
You can contact us at leonardo.pellegrina@unipd.it and fabio.vandin@unipd.it for any questions and for reporting bugs.

### Aknowledgments ###
We would like to thank the authors of KADABRA for making their code freely available; this truly helped us developing our algorithm.
