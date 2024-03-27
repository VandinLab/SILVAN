#include <iostream>
#include <dirent.h>
#include <limits.h>
#include <math.h>
#include <ctime>
#include <cstddef>
#include <unistd.h>
#include <ctype.h>
#include <fstream>

#include "Rand_gen.h"
#include "Graph.h"
#include "utilities.h"
#include "Probabilistic.h"

extern char *optarg;
static const std::string ERROR_HEADER = "ERROR: ";
using namespace std;

bool directed = false;
double verb = 60;
double delta;
double err;
char *graph_file;
std::string output_file;
int64_t k = 0;
double sampling_rate = 2.3;
bool alpha_given = false;
double empirical_peeling_param = 2.0;
bool m_hat_enabled = true;

// mcrade
int num_mc = 10;

/**
 * Print usage on stderr.
 */
void usage(const char *binary_name) {
    std::cerr << binary_name
        << ": compute betweenness centrality approximations for all nodes"
        << std::endl;
    std::cerr << "USAGE: " << binary_name << " [-dhm] [-v verbosity] [-k k_value] [-o output] [-a a_emp_peeling] [-s alpha] epsilon delta graph"
        << std::endl;
    std::cerr << "\t-d: consider the graph as directed" << std::endl;
    std::cerr << "\t-k: compute the top-k betweenness centralities (if 0, compute all of them with absolute error) " << std::endl;
    std::cerr << "\t-h: print this help message" << std::endl;
    std::cerr << "\t-v: print additional messages (verbosity is the time in second between subsequent outputs)" << std::endl;
    std::cerr << "\t-o: path for the output file (if empty, do not write the output)" << std::endl;
    std::cerr << "\t-a: parameter a for empirical peeling (def. = 2)" << std::endl;
    std::cerr << "\t-s: parameter alpha for sampling shortest paths (def. = 2.3)" << std::endl;
    std::cerr << "\t-m: disable the computation of m_hat" << std::endl;
    std::cerr << "\terr: accuracy (0 < epsilon < 1), relative accuracy if k > 0" << std::endl;
    std::cerr << "\tdelta: confidence (0 < delta < 1)" << std::endl;
    std::cerr << "\tgraph: graph edge list file" << std::endl;
}

/**
 * Parse command line options.
 * Return 0 if everything went well, 1 if there were errors, 2 if -h was specified.
 */
int parse_command_line(int& argc, char *argv[]) {
    int opt;
    while ((opt = getopt(argc, argv, "dhmk:o:s:a:v:")) != -1) {
        switch (opt) {
        case 'd':
            directed = true;
            break;
        case 'h':
            return 2;
            break;
        case 'm':
            m_hat_enabled = false;
            break;
        case 'o':
            std::cerr << "Writing output to " << optarg << std::endl;
            output_file = optarg;
            break;
        case 's':
            sampling_rate = std::strtod(optarg, NULL);
            alpha_given = true;
            break;
        case 'a':
            empirical_peeling_param = std::strtod(optarg, NULL);
            if (errno == ERANGE || empirical_peeling_param <= 1) {
                std::cerr << ERROR_HEADER
                    << "The value a should be >= 1. Empirical peeling disabled."
                    << std::endl;
            }
            alpha_given = true;
            break;
        case 'k':
            k = std::strtod(optarg, NULL);
            if (errno == ERANGE || k < 0 || k > UINT_MAX) {
                std::cerr << ERROR_HEADER
                    << "The value k should be between 0 and 2^32-1."
                    << std::endl;
                return 1;
            }
            break;
        case 'v':
            verb = std::strtod(optarg, NULL);
            if (errno == ERANGE || verb < 0) {
                std::cerr << ERROR_HEADER
                    << "The verbosity should be a positive number, or 0 to produce no output."
                    << std::endl;
                return 1;
            }
            break;
        }
    }

    if (optind != argc - 3) {
        std::cerr << ERROR_HEADER << "Wrong number of arguments" << std::endl;
        return 1;
    } else {
        err = std::strtod(argv[argc - 3], NULL);
        if (errno == ERANGE || err >= 1.0 || err <= 0.0) {
            std::cerr << ERROR_HEADER <<
                "The error err should be greater than 0 and smaller than 1"
                << std::endl;
            return 1;
        }
        delta = std::strtod(argv[argc - 2], NULL);
        if (errno == ERANGE || delta >= 1.0 || delta <= 0.0) {
            std::cerr << ERROR_HEADER <<
                "Delta should be greater than 0 and smaller than 1"
                << std::endl;
            return 1;
        }
        graph_file = argv[argc - 1];
        // test if input graph file exists
        std::ifstream infile(graph_file);
        if(!infile.good()){
          std::cerr << "Problems with input graph file" << std::endl;
          return 1;
        }
    }

    return 0;
}

int main(int argc, char *argv[]){
    int correct_parse = parse_command_line(argc, argv);

    if (correct_parse != 0) {
        usage(argv[0]);
        return correct_parse!=2;
    }

    Probabilistic G( graph_file, directed, verb , sampling_rate , alpha_given , empirical_peeling_param , m_hat_enabled , output_file);
    G.run((uint32_t) k, delta, err);
    std::cout << "run finished" << std::endl;
    return 0;
}
