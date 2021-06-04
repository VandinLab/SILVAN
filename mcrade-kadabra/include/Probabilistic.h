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
        Probabilistic( const std::string &filename, bool directed = false, const double verb = 60 );
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
        void compute_bet_err(Status *status, double *bet, double *err_l, double *err_u) const;
        void print_status(Status *Status, const bool full = false) const;
        void get_status (Status *status) const;
        void one_round(Sp_sampler &sp_sampler);
        bool compute_finished(Status *status) const;
        bool compute_finished_mcrade(Status *status , bool last_iter);
        double get_next_stopping_sample();
        double get_max_epsilon_mcrade(double maxestbc_high , double maxestbc_low  , double mcrade_low , bool verbose , bool last_iter) const;
        double check_topk_guarantees(double* mcrade_clusters , bool verbose , bool last_iter);
        void compute_delta_guess();
        double compute_f( const double btilde, const uint64_t iter_num, const double delta_l ) const;
        double compute_g( const double btilde, const uint64_t iter_num, const double delta_u ) const;

        double delta;
        double err;
        double last_output;
        bool absolute;
        uint32_t k;
        uint32_t union_sample;
        uint32_t start_factor;
        Ranking_list *top_k;
        double *approx;
        double *delta_l_guess;
        double *delta_u_guess;
        double delta_l_min_guess;
        double delta_u_min_guess;
        uint64_t n_pairs;
        uint64_t vis_edges;
        double start_time;
        double *time_bfs;
        double *time_critical;
        double *time_comp_finished;
        double omega;
        // mcrade
        int64_t *mcrade;
        Rand_gen *mcrade_randgen;
        double sup_mcera;
        uint64_t sup_mcera_node;
        double sup_bcest_supdiscr;
        double sup_bcest;
        uint64_t sup_bcest_node;
        uint64_t num_samples;
        int *highfreq_set;
        uint32_t mctrials = 100;
        int64_t* max_mcera_exact;
        int64_t* max_mcera_exact_low;
        int64_t* max_mcera_exact_high;
        uint64_t next_stopping_samples;
        double delta_progressive_current;
        int iteration_index;
        double num_freq;
        uint64_t distinct_nodes_top_k;
        double theta_top_k;
        int num_clusters_top_k;
        uint64_t *sup_bcest_cluster;
        int64_t *max_mcera_cluster;
        double eps_final_topk;
};


#endif
