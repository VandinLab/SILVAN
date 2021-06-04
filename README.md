# README #

This repository contains the implementation of SILVAN, supporting our paper ["SILVAN: Estimating Betweenness Centralities with Progressive Sampling and Non-uniform Rademacher Bounds"](http://www.dei.unipd.it/~pellegri/papers/silvan-ext.pdf).

SILVAN is an algorithm to compute approximations of the Betweenness Centralities from graphs with progressive sampling. More details can be found in the paper, available here: http://www.dei.unipd.it/~pellegri/papers/silvan-ext.pdf

The underlying implementation of SILVAN is based on the sampling algorithm [KADABRA](https://github.com/natema/kadabra) by Michele Borassi and Emanuele Natale. Therefore, it has the same compilation dependencies (see the `kadabra` subfolder for more details) and it is distributed with the same license.

### Running SILVAN ###

Coming soon! (a cleaner and more optimized version of SILVAN is in the process)

### Reproducing experiments
To reproduce experiments described in the paper, follow the instructions listed below.

1. First you need to compile the algorithms. You can do so with the command `make` within the folders `kadabra` and `mcrade-kadabra`. The latter contains the implementation of a modified version of KADABRA that also runs SILVAN (a more optimized version of SILVAN is in the process).
2. Then, move to the `script` folder. Download all graphs using `python download_ds.py`. Graphs will be downloaded in the `datasets` folder.
3. To run experiments described in Section 5.1, run `python run_experiments.py`. Results will be appended in the file `results.csv`
4. To run experiments described in Section 5.2, run `python run_topk_experiments_all.py`. Results will be appended in the file `results_topk_new.csv` for SILVAN, and in files `results_topk_new.csv` and `results_topk_new.csv` for the variants of KADABRA described in the paper (resp., K_k and K_\varepsilon).

### Contacts ###
You can contact us at pellegri@dei.unipd.it for any questions and for reporting bugs.

### Aknowledgments ###
We would like to thank the authors of KADABRA for making their code freely available; this truly helped us developing our algorithm.
