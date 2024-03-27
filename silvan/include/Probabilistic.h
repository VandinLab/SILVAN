#ifndef PROBABILISTIC_H
#define PROBABILISTIC_H

#include <Graph.h>
#include <string>
#include <utility>
#include <vector>
#include "Ranking_list.h"
#include "Rand_gen.h"
#include "Sp_sampler.h"

class Status {
    public:
        Status( const uint32_t k );
        virtual ~Status();
        const uint32_t k;
        uint32_t *top_k;
        double *approx_top_k;
        uint64_t n_pairs;
        bool *finished;

        double *bet;
        double *err_l;
        double *err_u;
};

class Probabilistic : public Graph
{
    public:
        Probabilistic( const std::string &filename, bool directed = false, const double verb = 60, const double sampling_rate_ = 2.3, bool alpha_given_ = false, const double empirical_peeling_param_ = 2.0 , const bool enable_m_hat_ = true, const std::string output_file_ = "" );
        virtual ~Probabilistic();
        void run(const uint32_t k, const double delta, const double err = 0,
                 const uint32_t union_sample = 0,
                 const uint32_t start_factor = 100);
        inline double get_centrality(const uint32_t v) const {
            return (double) approx[v] / n_pairs;
        }
        inline long get_n_pairs() const {
            return n_pairs;
        }
        inline long get_vis_edges() const {
            return vis_edges;
        }
        double verbose = 60;
    protected:
    private:
        void print_status(Status *Status, const bool full = false) const;
        void get_status (Status *status) const;
        void one_round(Sp_sampler &sp_sampler);
        bool compute_finished_mcrade(Status *status);
        double get_next_stopping_sample();
        double get_epsilon_mcrade(double sup_emp_wimpy_var_ , double mcera_ , double delta_ , double num_samples_ , double num_nodes_ , bool verbose) const;
        double check_topk_guarantees(bool verbose);
        double getUpperBoundTop1BC(double top1_est_bc , double delta);
        double getUpperBoundAvgDiameter(double delta , bool verbose);
        double computeRelBound(double bc_est , double delta , double rho , double m_samples , bool upper);


        double delta;
        double err;
        double last_output;
        bool absolute;
        uint32_t k;
        uint32_t union_sample;
        Ranking_list *top_k;
        double *approx;
        uint64_t n_pairs;
        uint64_t vis_edges;
        double start_time;
        double *time_bfs;
        double *time_critical;
        double *time_critical_round;
        double *time_comp_finished;
        double *time_mcera;
        double omega;
        // mcrade
        int graph_diameter;
        double last_stopping_samples;
        double first_stopping_samples;
        bool second_phase_started;
        int64_t *mcrade;
        double *emp_wimpy_vars;
        double *approx_toadd;
        Rand_gen *mcrade_randgen;
        double sup_emp_wimpy_var;
        double sup_bcest;
        uint64_t num_samples;
        int *partition_index;
        uint32_t mctrials = 25;
        bool enable_m_hat;
        bool enable_emp_peel;
        double delta_for_progressive_bound;
        uint64_t next_stopping_samples;
        int iteration_index;
        uint64_t distinct_nodes_top_k;
        uint64_t *sup_bcest_partition;
        double *sup_empwvar_partition;
        int64_t *max_mcera_partition;
        double *epsilon_partition;
        double eps_final_topk;
        bool firstpass;
        int64_t void_samples;
        uint64_t *sp_lengths;
        double alpha_sp_sampling;
        bool alpha_sp_given;
        double empirical_peeling_a;
        int number_of_non_empty_partitions;
        std::map<int, int> partitions_ids_map;
        std::string output_file;
        uint32_t numresults_topk;
};


#endif
