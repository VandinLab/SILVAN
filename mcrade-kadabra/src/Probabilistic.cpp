#include <string>
#include <iostream>
#include <math.h>
#include <limits.h>
#include <iomanip>
#include <time.h>
#include <cfloat>
#include <omp.h>

#include "utilities.h"
#include "Probabilistic.h"
#include "Sp_sampler.h"

#define SEED 42


using namespace std;

// The status class contains the data about the k most central vertices.
Status::Status(const uint32_t k) : k(k) {
    approx_top_k = (double *) malloc( k*sizeof(double));
    top_k = (uint32_t *) malloc( k*sizeof(uint32_t) );
    finished = (bool *) malloc( k*sizeof(bool) );
    bet = (double*) malloc(k * sizeof(double));
    err_l = (double*) malloc(k * sizeof(double));
    err_u = (double*) malloc(k * sizeof(double));
}

Status::~Status() {
    free(approx_top_k);
    free(top_k);
    free(finished);
    free(bet);
    free(err_l);
    free(err_u);
}

// Creates the graph for running the approximation algorithm.
// For more information see the graph class.
Probabilistic::Probabilistic( const std::string &filename, const bool directed, const double verb ): Graph( filename, directed ), verbose(verb) {
    approx = (double *) calloc( get_nn(), sizeof(double) );
    delta_l_guess = (double *) calloc( get_nn(), sizeof(double) );
    delta_u_guess = (double *) calloc( get_nn(), sizeof(double) );
    time_bfs = (double *) calloc( omp_get_max_threads(), sizeof(double) );
    time_comp_finished = (double *) calloc( omp_get_max_threads(), sizeof(double) );
    time_critical = (double *) calloc( omp_get_max_threads(), sizeof(double) );
    n_pairs = 0;
    vis_edges = 0;
    if (verbose > 0) {
        print_data();
    }
    // mcrade
    mcrade = (int64_t *) calloc( get_nn()*mctrials, sizeof(int64_t) );
    highfreq_set = (int *) calloc( get_nn(), sizeof(int) );
    mcrade_randgen = new Rand_gen( 2021 );
    max_mcera_exact = (int64_t*) calloc( mctrials , sizeof(int64_t) );
    max_mcera_exact_low = (int64_t*) calloc( mctrials , sizeof(int64_t) );
    max_mcera_exact_high = (int64_t*) calloc( mctrials , sizeof(int64_t) );
    sup_mcera = 0.0;
    sup_bcest = 0.0;
}

// Decides whether the algorithm should terminate
// INPUT: a Status object describing the current status of the algorithm
bool Probabilistic::compute_finished(Status *status) const {

    double *bet = status->bet;
    double *err_l = status->err_l;
    double *err_u = status->err_u;
    bool all_finished = true;

    uint32_t i;
    for (i = 0; i < status->k-1; i++) {
        bet[i] = status->approx_top_k[i] / status->n_pairs;
        err_l[i] = compute_f( bet[i], status->n_pairs, delta_l_guess[status->top_k[i]] );
        err_u[i] = compute_g( bet[i], status->n_pairs, delta_u_guess[status->top_k[i]] );
    }
    bet[i] = status->approx_top_k[i] / status->n_pairs;
    err_l[i] = compute_f( bet[i], status->n_pairs, this->delta_l_min_guess );
    err_u[i] = compute_g( bet[i], status->n_pairs, this->delta_u_min_guess );

    if (absolute) {
        for (uint32_t i = 0; i < status->k; i++) {
            status->finished[i] = (err_l[i] < err && err_u[i] < err);
            all_finished = all_finished && status->finished[i];
        }
    } else {
        for (uint32_t i = 0; i < status->k; i++) {
            if (i == 0) {
                status->finished[i] = (bet[i]-err_l[i] > bet[i+1]+err_u[i+1]);
            } else if (i < k) {
                status->finished[i] = (bet[i-1]-err_l[i-1] > bet[i]+err_u[i]) && (bet[i]-err_l[i] > bet[i+1]+err_u[i+1]);
            } else {
                status->finished[i] = bet[k-1]-err_u[k-1] > bet[i]+err_u[i];
            }
            status->finished[i] = status->finished[i] || (err_l[i] < err && err_u[i] < err);
            all_finished = all_finished && status->finished[i];
        }
    }

    return all_finished;
}


// Decides whether the algorithm should terminate
// according to mcera
// INPUT: a Status object describing the current status of the algorithm
bool Probabilistic::compute_finished_mcrade(Status *status , bool last_iter) {

  std::cout << std::endl;
  cout << "EVALUATING STOPPING CONDITION at iteration " << num_samples << std::endl;

  // first, reset counters for computation of mcera
  if(absolute){
    for(uint32_t j = 0; j < mctrials; j++){
      max_mcera_exact[j] = -(int64_t)status->n_pairs;
      max_mcera_exact_low[j] = -(int64_t)status->n_pairs;
      max_mcera_exact_high[j] = -(int64_t)status->n_pairs;
    }
  }
  else{
    for(int i = 0; i < num_clusters_top_k; i++){
      for(uint32_t j = 0; j < mctrials; j++){
        max_mcera_cluster[i*mctrials+j] = -(int64_t)status->n_pairs;
      }
    }
  }
  //uint64_t max_mcera_node = 0;
  double sup_bcest_high_freq = 0.0;
  double sup_mcera_high_freq = 0.0;
  double sup_bcest_low_freq = 0.0;
  double sup_mcera_low_freq = 0.0;
  // iterate over nodes to update the mcera
  for (uint32_t i = 0; i < status->k; i++) {
      uint32_t v = status->top_k[i];
      uint32_t v_rade_idx = v*mctrials;
      if(absolute){
        for(uint32_t j = 0; j < mctrials; j++){
          max_mcera_exact[j] = max(max_mcera_exact[j] , mcrade[v_rade_idx+j]);
        }
        if(highfreq_set[v] == 1){
          // this node was selected to be high freq
          sup_bcest_high_freq = max(sup_bcest_high_freq , approx[v]);
          //sup_mcera_high_freq = max(sup_mcera_high_freq , mcrade[v]);
          for(uint32_t j = 0; j < mctrials; j++){
            max_mcera_exact_high[j] = max(max_mcera_exact_high[j] , mcrade[v_rade_idx+j]);
          }
        }
        else{
          // this node has low freq
          sup_bcest_low_freq = max(sup_bcest_low_freq , approx[v]);
          //sup_mcera_low_freq = max(sup_mcera_low_freq , mcrade[v]);
          for(uint32_t j = 0; j < mctrials; j++){
            max_mcera_exact_low[j] = max(max_mcera_exact_low[j] , mcrade[v_rade_idx+j]);
          }
        }
      }
      // top-k algorithm
      else{
        int cluster_index = highfreq_set[v];
        sup_bcest_cluster[cluster_index] = max(sup_bcest_cluster[cluster_index] , (uint64_t)approx[v]);
        int mcera_cluster_index = mctrials*cluster_index;
        for(uint32_t j = 0; j < mctrials; j++){
          max_mcera_cluster[mcera_cluster_index+j] = max(max_mcera_cluster[mcera_cluster_index+j] , mcrade[v_rade_idx+j]);
        }
      }
      //max_mcera_exact = max(mcrade[v] , max_mcera_exact);
  }

  // the mcera has been computed , now evaluate the stopping condition
  // stopping condition for absolute (uniform) approximation
  if(absolute){
    cout << "max_mcera_exact: " << max_mcera_exact << std::endl;
    //cout << "max_mcera_node: " << max_mcera_node << std::endl;

    // compute averages of supremum discrepancies
    sup_mcera_high_freq = 0.0;
    cout << "  max_mcera_exact_high \n    ";
    for(uint32_t j = 0; j < mctrials; j++){
      sup_mcera_high_freq += max_mcera_exact_high[j];
      cout << " " << max_mcera_exact_high[j] << " ";
    }
    cout << std::endl;
    sup_mcera_high_freq = sup_mcera_high_freq / mctrials;
    sup_mcera_low_freq = 0.0;
    cout << "  max_mcera_exact_low \n    ";
    for(uint32_t j = 0; j < mctrials; j++){
      sup_mcera_low_freq += max_mcera_exact_low[j];
      cout << " " << max_mcera_exact_low[j] << " ";
    }
    cout << std::endl;
    sup_mcera_low_freq = sup_mcera_low_freq / mctrials;
    double sup_mcera_all = 0.0;
    for(uint32_t j = 0; j < mctrials; j++){
      sup_mcera_all += max_mcera_exact[j];
    }
    sup_mcera_all = sup_mcera_all / mctrials;

    cout << "sup_bcest_high_freq: " << sup_bcest_high_freq << std::endl;
    cout << "sup_mcera_high_freq: " << sup_mcera_high_freq << std::endl;
    cout << "sup_bcest_low_freq: " << sup_bcest_low_freq << std::endl;
    cout << "sup_mcera_low_freq: " << sup_mcera_low_freq << std::endl;
    cout << "sup_mcera_all: " << sup_mcera_all << std::endl;

    double num_samples_d = (double)num_samples;
    double maxestbc_high , maxestbc_low  , mcrade_low;
    maxestbc_high = sup_bcest_high_freq/num_samples_d;
    maxestbc_low = sup_bcest_low_freq/num_samples_d;
    mcrade_low = max(0. , sup_mcera_low_freq/num_samples_d);
    cout << "num_samples_d: " << num_samples_d << std::endl;
    cout << "    maxestbc_low: " << maxestbc_low << std::endl;

    double sup_eps = get_max_epsilon_mcrade(maxestbc_high , maxestbc_low  , mcrade_low , true, last_iter);

    // if we at the last iteration and MCRade has not converged yet
    if(last_iter == true){
      // this means that MCRade has not converged, while kadabra has
      if(sup_eps > err){
        cout << "MCRADE STOPS WITH eps " << err << std::endl;
        cout << "MCRADE STOPS at iteration " << iteration_index << std::endl;
        cout << "MCRADE STOPS at sample size " << (int)max((double)num_samples,omega) << std::endl;
        cout << "MCRADE STOPS after seconds " << get_time_sec() - start_time << std::endl;
      }
    }

    if(sup_eps <= err){
      cout << "MCRADE STOPS WITH eps " << sup_eps << std::endl;
      cout << "MCRADE STOPS at iteration " << iteration_index << std::endl;
      cout << "MCRADE STOPS at sample size " << num_samples << std::endl;
      cout << "MCRADE STOPS after seconds " << get_time_sec() - start_time << std::endl;
    }

    return sup_eps <= err;
  }
  // stopping condition for top-k approximation
  else{
    // this array containes the mcera for every cluster
    double *mcera_clusters_avg = (double*)calloc( (uint32_t)num_clusters_top_k , sizeof(double));
    double mcera_avg_ = 0.0;
    int mcera_cluster_index = 0;
    for(int i = 0; i < num_clusters_top_k; i++){
      mcera_avg_ = 0.0;
      mcera_cluster_index = mctrials*i;
      for(uint32_t j = 0; j < mctrials; j++){
        mcera_avg_ += max_mcera_cluster[mcera_cluster_index+j]/(double)mctrials;
      }
      mcera_clusters_avg[i] = mcera_avg_;
    }

    // check if the guarantees on the top-k are satisfied
    bool top_k_done = check_topk_guarantees(mcera_clusters_avg , true, last_iter);
    // output the results
    if(top_k_done){
        cout << "MCRADE REL STOPS " << std::endl;
        cout << "MCRADE STOPS at iteration " << iteration_index << std::endl;
        cout << "MCRADE STOPS at sample size " << num_samples << std::endl;
        cout << "MCRADE STOPS after seconds " << get_time_sec() - start_time << std::endl;
    }
    else{
        cout << "MCRADE REL DOES NOT STOP " << std::endl;
    }


    free(mcera_clusters_avg);

    return top_k_done;

  }

}

// Computes the function f that bounds the betweenness of a vertex from below.
// For more information, see Borassi, Natale (2016).
double Probabilistic::compute_f( const double btilde, const uint64_t iter_num, const double delta_l ) const {
    double tmp = (((double) omega) / iter_num - 1./3);
    double err_chern = (log(1./delta_l)) * 1./iter_num * (-tmp + sqrt(tmp * tmp + 2 * btilde * omega / (log(1./delta_l))));
    return min(err_chern, btilde);
}

// Computes the function g that bounds the betweenness of a vertex from above.
// For more information, see Borassi, Natale (2016).
double Probabilistic::compute_g( const double btilde, const uint64_t iter_num, const double delta_u ) const {
    double tmp = (((double) omega) / iter_num + 1./3);
    double err_chern = (log(1./delta_u)) * 1./iter_num * (tmp + sqrt(tmp * tmp + 2 * btilde * omega / (log(1./delta_u))));
    return min(err_chern, 1-btilde);
}

// Outputs the current status.
// INPUT: a Status object describing the current status, and a flag "full".
// If full is true, we output more data.
void Probabilistic::print_status(Status *status, const bool full) const {
    if (full) {
        std::cout << std::setprecision(6) << endl << "Finished after " << status->n_pairs << " iterations." << endl;
    } else {
        std::cout << std::setprecision(6) << endl << "Situation after " << status->n_pairs << " iterations." << endl;
    }
    std::cout << "Edges visited: " << vis_edges << endl;
    std::cout << "Average edges visited: " << vis_edges/status->n_pairs << endl;
    std::cout << "Total time: " << get_time_sec() - start_time << endl;
    std::cout << "Time bfs: " << time_bfs[omp_get_thread_num()] << endl;
    std::cout << "Time critical: " << time_critical[omp_get_thread_num()] << endl;
    std::cout << "Time compute finished: " << time_comp_finished[omp_get_thread_num()] << endl;
    std::cout << "(Printing thread: " << omp_get_thread_num() << ")" << endl;
    // mcera
    //std::cout << "mcera est: " << sup_mcera << endl;
    std::cout << "sup bc est: " << sup_bcest << endl;
    //std::cout << "sup_mcera_node: " << sup_mcera_node << endl;
    //std::cout << "sup_bcest_node: " << sup_bcest_node << endl;
    //std::cout << "mcrade[sup_mcera_node*mctrials]: " << mcrade[sup_mcera_node*mctrials] << endl;
    //std::cout << "mcrade[sup_bcest_node*mctrials]: " << mcrade[sup_bcest_node*mctrials] << endl;
    //std::cout << "approx[sup_mcera_node*mctrials]: " << approx[sup_mcera_node*mctrials] << endl;
    //std::cout << "sup_bcest_supdiscr: " << sup_bcest_supdiscr << endl;
    //std::cout << "approx[sup_bcest_node]: " << approx[sup_bcest_node] << endl;
    std::cout << "num_samples: " << num_samples << endl;


    //#define DEBUG_DELTAS 1

    if (absolute) {
        double max_interval = 0;
        double delta_tot = 0.0;
        for(uint32_t j = 0; j < mctrials; j++){
          max_mcera_exact[j] = -(int64_t)status->n_pairs;
          max_mcera_exact_low[j] = -(int64_t)status->n_pairs;
          max_mcera_exact_high[j] = -(int64_t)status->n_pairs;
        }
        //uint64_t max_mcera_node = 0;
        double sup_bcest_high_freq = 0.0;
        double sup_mcera_high_freq = 0.0;
        double sup_bcest_low_freq = 0.0;
        double sup_mcera_low_freq = 0.0;
        #ifdef DEBUG_DELTAS
        double lower_int = 0.0;
        double upp_int = 0.0;
        #endif
        for (uint32_t i = 0; i < status->k; i++) {
            uint32_t v = status->top_k[i];
            max_interval = max(max_interval, compute_f(status->approx_top_k[i] / status->n_pairs, status->n_pairs, delta_l_guess[v]));
            max_interval = max(max_interval, compute_g(status->approx_top_k[i] / status->n_pairs, status->n_pairs, delta_u_guess[v]));
            delta_tot += delta_l_guess[v]+delta_u_guess[v];
            uint32_t v_rade_idx = v*mctrials;
            for(uint32_t j = 0; j < mctrials; j++){
              max_mcera_exact[j] = max(max_mcera_exact[j] , mcrade[v_rade_idx+j]);
            }
            if(highfreq_set[v] == 1){
              // this node was selected to be high freq
              sup_bcest_high_freq = max(sup_bcest_high_freq , approx[v]);
              //sup_mcera_high_freq = max(sup_mcera_high_freq , mcrade[v]);
              for(uint32_t j = 0; j < mctrials; j++){
                max_mcera_exact_high[j] = max(max_mcera_exact_high[j] , mcrade[v_rade_idx+j]);
              }
            }
            else{
              // this node has low freq
              sup_bcest_low_freq = max(sup_bcest_low_freq , approx[v]);
              //sup_mcera_low_freq = max(sup_mcera_low_freq , mcrade[v]);
              for(uint32_t j = 0; j < mctrials; j++){
                max_mcera_exact_low[j] = max(max_mcera_exact_low[j] , mcrade[v_rade_idx+j]);
              }
            }
            //max_mcera_exact = max(mcrade[v] , max_mcera_exact);

            #ifdef DEBUG_DELTAS
            lower_int = compute_f(status->approx_top_k[i] / status->n_pairs, status->n_pairs, delta_l_guess[v]);
            upp_int = compute_g(status->approx_top_k[i] / status->n_pairs, status->n_pairs, delta_u_guess[v]);
            if( lower_int > err){
              std::cout << " upper interval higher than err " << lower_int << std::endl;
              std::cout << " delta_l_guess[v] " << delta_l_guess[v] << std::endl;
            }
            if( upp_int > err){
              std::cout << " lower interval higher than err " << upp_int << std::endl;
              std::cout << " delta_u_guess[v] " << delta_u_guess[v] << std::endl;
            }
            #endif
        }
        cout << "delta_tot: " << delta_tot << std::endl;
        cout << "max_mcera_exact: " << max_mcera_exact << std::endl;
        //cout << "max_mcera_node: " << max_mcera_node << std::endl;

        // compute averages of supremum discrepancies
        sup_mcera_high_freq = 0.0;
        cout << "  max_mcera_exact_high \n    ";
        for(uint32_t j = 0; j < mctrials; j++){
          sup_mcera_high_freq += max_mcera_exact_high[j];
          cout << " " << max_mcera_exact_high[j] << " ";
        }
        cout << std::endl;
        sup_mcera_high_freq = sup_mcera_high_freq / mctrials;
        sup_mcera_low_freq = 0.0;
        cout << "  max_mcera_exact_low \n    ";
        for(uint32_t j = 0; j < mctrials; j++){
          sup_mcera_low_freq += max_mcera_exact_low[j];
          cout << " " << max_mcera_exact_low[j] << " ";
        }
        cout << std::endl;
        sup_mcera_low_freq = sup_mcera_low_freq / mctrials;
        double sup_mcera_all = 0.0;
        for(uint32_t j = 0; j < mctrials; j++){
          sup_mcera_all += max_mcera_exact[j];
        }
        sup_mcera_all = sup_mcera_all / mctrials;

        cout << "sup_bcest_high_freq: " << sup_bcest_high_freq << std::endl;
        cout << "sup_mcera_high_freq: " << sup_mcera_high_freq << std::endl;
        cout << "sup_bcest_low_freq: " << sup_bcest_low_freq << std::endl;
        cout << "sup_mcera_low_freq: " << sup_mcera_low_freq << std::endl;
        cout << "sup_mcera_all: " << sup_mcera_all << std::endl;

        //cout << "approx[max_mcera_node]: " << approx[max_mcera_node] << std::endl;
        cout << "Maximum confidence interval: " << max_interval;
    }
    else {
        compute_finished(status);
        uint32_t i;
        for (i = 0; i < k; i++) {
            double bet = status->approx_top_k[i] / status->n_pairs;
            if (status->finished[i]) {
                cout << std::setw(8) << to_string(i+1) << ") ";
            } else {
                cout << std::setw(8) << "? " + to_string(i+1) << ") ";
            }
            cout << std::setw(8) << status->top_k[i] << " " << bet-compute_f(bet, status->n_pairs, delta_l_guess[status->top_k[i]]) << " ";
            cout << bet << " " << bet+compute_g(bet, status->n_pairs, delta_u_guess[status->top_k[i]]) << endl;
        }
        if (full) {
            double betk = status->approx_top_k[k-1] / status->n_pairs;
            double lbetk = betk - compute_f(betk, status->n_pairs, delta_l_guess[status->top_k[k-1]]);
            uint32_t pos = k+1;
            for (i = k; i < status->k; i++) {
                double bet = status->approx_top_k[i] / status->n_pairs;
                if (bet+compute_g(bet, status->n_pairs, delta_u_guess[status->top_k[i]]) > lbetk) {
                    cout << std::setw(8) << to_string(pos++) << ") ";
                    cout << std::setw(8) << status->top_k[i] << " " << bet-compute_f(bet, status->n_pairs, delta_l_guess[status->top_k[i]]) << " ";
                    cout << bet << " " << bet+compute_g(bet, status->n_pairs, delta_u_guess[status->top_k[i]]) << endl;
                }
            }
        } else {
            double max_upper = 0;
            for (i = k; i < status->k; i++) {
                double bet = status->approx_top_k[i] / status->n_pairs;
                max_upper = max(max_upper, bet+compute_g(bet, status->n_pairs, delta_u_guess[status->top_k[i]]));
            }
            double bet = status->approx_top_k[status->k-1] / status->n_pairs;
            max_upper = max(max_upper, bet+compute_g(bet, status->n_pairs, delta_u_min_guess));

            cout << std::setw(8) << "Others" << ") <" << max_upper;
        }
    }
    cout << endl;
}

//#define ROUND_DEBUG 1

// Sample one shortest path and updates the ranking of the betweenness approximations.
void Probabilistic::one_round(Sp_sampler &sp_sampler) {
    time_bfs[omp_get_thread_num()] -= get_time_sec();
    vector<uint32_t> path = sp_sampler.random_path();
    time_bfs[omp_get_thread_num()] += get_time_sec();


    time_critical[omp_get_thread_num()] -= get_time_sec();
    #pragma omp critical
    {
        n_pairs++;
        vis_edges += sp_sampler.vis_edges;
        num_samples++;

        // mcrade
        uint64_t maxval_sigmas = 100000000;
        int* sigmas = (int*) calloc( mctrials , sizeof(int));
        for(uint32_t j = 0; j < mctrials; j++){
          sigmas[j] = (mcrade_randgen->get_max(maxval_sigmas) >= (double)maxval_sigmas/2.)*2-1;
          //std::cout << "sigmas[j] " << sigmas[j] << std::endl;
        }

        for(uint32_t u:path){
            approx[u]++;
            top_k->put(u, approx[u]);
            // mcrade
            uint32_t u_idx = u*mctrials;
            for(uint32_t j = 0; j < mctrials; j++){
              mcrade[u_idx+j] += sigmas[j];
            }
            distinct_nodes_top_k += approx[u] == 3;
            //mcrade[u]+=sigma;
            /*if(mcrade[u] > sup_mcera){
              sup_mcera = mcrade[u];
              sup_mcera_node = u;
              //sup_bcest_supdiscr = approx[u];
              #ifdef ROUND_DEBUG
              std::cout << "updating sup_mcera at iter " << num_samples << " ! node " << u << std::endl;
              std::cout << "  mcrade[u] " << mcrade[u] << std::endl;
              std::cout << "  approx[u] " << approx[u] << std::endl;
              std::cout << "    sup_bcest " << sup_bcest << std::endl;
              std::cout << "    sup_mcera " << sup_mcera << std::endl;
              #endif
            }*/
            if(approx[u] > sup_bcest){
              sup_bcest = approx[u];
              sup_bcest_node = u;
              /*std::cout << "sup_bcest " << sup_bcest << std::endl;
              std::cout << "sup_mcera " << sup_mcera << std::endl;
              std::cout << "mcrade[u] " << mcrade[u] << std::endl;
              std::cout << "approx[u] " << approx[u] << std::endl;*/
            }

        }

        free(sigmas);
    }

    time_critical[omp_get_thread_num()] += get_time_sec();
}

// Fills the input variable Status in a synchronized way.
void Probabilistic::get_status (Status *status) const {
    time_critical[omp_get_thread_num()] -= get_time_sec();
    #pragma omp critical
    {
        if (status != NULL) {
            for(uint32_t i=0; i<union_sample; i++) {
                status->top_k[i] = top_k->get(i);
                status->approx_top_k[i] = approx[status->top_k[i]];
            }
            status->n_pairs = n_pairs;
        }
    }
    time_critical[omp_get_thread_num()] += get_time_sec();
}

// Compute the values of err from which the *best* deltas are computed.
// The results are stored in err_l and err_u.
void Probabilistic::compute_bet_err(Status *status, double *bet, double *err_l, double *err_u) const {

    uint32_t i;
    double max_err = sqrt(start_factor) * err / 4;

    for (i = 0; i < status->k; i++) {
        bet[i] = status->approx_top_k[i] / status->n_pairs;
    }
    if (absolute) {
        for (i = 0; i < status->k; i++) {
            err_l[i] = err;
            err_u[i] = err;
        }
    } else {
        err_u[0] = max(err, (bet[0] - bet[1]) / 2.);
        err_l[0] = 10;
        for (i = 1; i < k; i++) {
            err_l[i] = max(err, (bet[i-1]-bet[i]) / 2.);
            err_u[i] = max(err, (bet[i]-bet[i+1]) / 2.);
        }
        for (i = k; i < status->k; i++) {
            err_l[i] = 10;
            err_u[i] = max(err, bet[k-1] + (bet[k-1]-bet[k]) / 2. - bet[i]);
        }
        for (i = 0; i < k-1; i++) {
            if (bet[i] - bet[i+1] < max_err) {
                err_l[i] = err;
                err_u[i] = err;
                err_l[i+1] = err;
                err_u[i+1] = err;
            }
        }
        for (i = k+1; i < status->k; i++) {
            if (bet[k] - bet[i] < max_err) {
                err_l[k] = err;
                err_u[k] = err;
                err_l[i] = err;
                err_u[i] = err;
            }
        }
    }
}

// Compute the *best* deltas for minimizing the stopping time of the algorithm.
// The computation is based on the heuristic in the paper Borassi, Natale (2016).
void Probabilistic::compute_delta_guess() {
    double balancing_factor = 0.001;
    double a = 0, b = 1. / err / err * log(get_nn() * 4 * (1-balancing_factor) / delta), c=(a+b)/2;
    double sum;
    Status status(union_sample);
    get_status(&status);

    double *bet = (double*) malloc(status.k * sizeof(double));
    double *err_l = (double*) malloc(status.k * sizeof(double));
    double *err_u = (double*) malloc(status.k * sizeof(double));

    compute_bet_err(&status, bet, err_l, err_u);
    for (uint32_t i = 0; i < union_sample; i++) {
        uint32_t v = status.top_k[i];
        approx[v] = approx[v] / n_pairs;
    }

    while (b-a > err/10) {
        c = (b+a)/2;
        sum = 0;
        for (uint32_t i = 0; i < union_sample; i++) {
            sum += exp(-c * err_l[i] * err_l[i] / bet[i]);
            sum += exp(-c * err_u[i] * err_u[i] / bet[i]);
        }
        sum += exp(-c * err_l[union_sample-1] * err_l[union_sample-1] / bet[union_sample-1]) * (get_nn() - union_sample);
        sum += exp(-c * err_u[union_sample-1] * err_u[union_sample-1] / bet[union_sample-1]) * (get_nn() - union_sample);

        if (sum >= delta / 2 *(1-balancing_factor)) {
            a = c;
        } else {
            b = c;
        }
    }
    delta_l_min_guess = exp(-b * err_l[union_sample-1] * err_l[union_sample-1] / bet[union_sample-1]) + delta * balancing_factor / 4. / get_nn();;
    delta_u_min_guess = exp(-b * err_u[union_sample-1] * err_u[union_sample-1] / bet[union_sample-1]) + delta * balancing_factor / 4. / get_nn();;

    std::cout << "delta_l_min_guess " << delta_l_min_guess << std::endl;
    std::cout << "delta_u_min_guess " << delta_u_min_guess << std::endl;
    std::cout << "err " << err << std::endl;

    for (uint32_t v = 0; v < get_nn(); v++) {
        delta_l_guess[v] = delta_l_min_guess;
        delta_u_guess[v] = delta_u_min_guess;
    }
    for (uint32_t i = 0; i < union_sample; i++) {
        uint32_t v = status.top_k[i];
        delta_l_guess[v] = exp(-b * err_l[i] * err_l[i] / bet[i]) + delta * balancing_factor / 4. / get_nn();;
        delta_u_guess[v] = exp(-b * err_u[i] * err_u[i] / bet[i]) + delta * balancing_factor / 4. / get_nn();;
    }
    free(bet);
    free(err_l);
    free(err_u);
}


// update the stopping condition schedule
double Probabilistic::get_next_stopping_sample(){
  next_stopping_samples = next_stopping_samples*1.2;
  if(next_stopping_samples*1.2 < (double)omega){
    iteration_index += 1;
  }
  return next_stopping_samples;
}

// function to check that the guarantees on the top-k are satisfied
double Probabilistic::check_topk_guarantees(double* mcrade_clusters , bool verbose , bool last_iter){

  // total delta for this iteration, to split across clusters
  double delta_for_progressive_bound = delta/pow(2.,iteration_index);
  double num_samples_d = (double)num_samples;
  double *eps_clusters = (double*)calloc(num_clusters_top_k , sizeof(double));
  std::cout << "EVALUATING STOPPING CONDITION FOR TOP-K RELATIVE APPROX" << std::endl;

  // compute maximum deviation for each cluster
  for(int i=0; i<num_clusters_top_k; i++){
    // we give higher deltas to cluster with lower values of BC
    double delta_for_this_cluster = delta_for_progressive_bound/pow(2.,i+1);
    // compute the maximum deviation
    double maxestbc_low = (double)sup_bcest_cluster[i]/num_samples_d;
    double mcrade_low = mcrade_clusters[i]/num_samples;
    mcrade_low = max(0. , mcrade_low);
    double log_term_mcrade = log(5./delta_for_this_cluster)/num_samples_d;
    double max_bc_low = maxestbc_low + log_term_mcrade + sqrt( pow(log_term_mcrade,2) + 2.*maxestbc_low*log_term_mcrade  );
    double var_bound_low = (max_bc_low > 0.5) ? 0.25 : max_bc_low*(1.-max_bc_low);
    double log_term_mcrade_n = log_term_mcrade/(double)mctrials;
    double era_low = mcrade_low + 2.*log_term_mcrade_n + sqrt( pow(2.*log_term_mcrade_n,2) + 4.*(mcrade_low+max_bc_low)*log_term_mcrade_n  );
    double rc_low = era_low + log_term_mcrade + sqrt( pow(log_term_mcrade,2) + 2.*era_low*log_term_mcrade );
    double eps_low_mcrade = 2.*rc_low + sqrt( 2.*log_term_mcrade*(var_bound_low+4.*rc_low) ) + log_term_mcrade/3.;
    eps_clusters[i] = eps_low_mcrade;
    std::cout << "  eps_clusters["<<i<<"] " << eps_clusters[i] << std::endl;
    std::cout << "    delta_for_this_cluster " << delta_for_this_cluster << std::endl;
    std::cout << "    mcrade " << mcrade_low << std::endl;
    std::cout << "    maxestbc " << maxestbc_low << std::endl;
    std::cout << "    max_bc " << max_bc_low << std::endl;
    std::cout << "    rc " << rc_low << std::endl;
  }

  // now that we have the sup deviations , we have to check that the guarantees
  // on the top-k approximation are satisfied
  // iterate over nodes in order of estimated BC
  double lb_topk_bc = 0.0;
  double ub_topk_bc = 1.0;
  Ranking_list* sorted_ub_topk = new Ranking_list(k);
  Ranking_list* sorted_lb_topk = new Ranking_list(k);
  uint32_t num_inserted = 0;
  bool is_relative_approx = true;
  for(uint64_t i = 0; i < union_sample && is_relative_approx; i++){
    uint64_t node_id = this->top_k->get(i);
    double approx_node = this->top_k->get_value(i)/num_samples_d;
    int node_cluster_index = highfreq_set[node_id];
    double eps_current_node = eps_clusters[node_cluster_index];
    // we have to check that this eps gives a relative approximation
    // according to the input parameter err
    double lowerbound_bc_ = approx_node-eps_current_node;
    double upperbound_bc_ = approx_node+eps_current_node;
    sorted_ub_topk->put(node_id , upperbound_bc_);
    sorted_lb_topk->put(node_id , lowerbound_bc_);
    num_inserted++;
    if (num_inserted >= k){
      ub_topk_bc = min(ub_topk_bc , sorted_ub_topk->get_value(k-1));
      lb_topk_bc = max(lb_topk_bc , sorted_lb_topk->get_value(k-1));
    }
    // this node should go in the output as a top-k node
    if(upperbound_bc_ > lb_topk_bc){
      bool lower_condition = lowerbound_bc_ >= approx_node/(1.+err);
      bool upper_condition = upperbound_bc_ <= approx_node/(1.-err);
      if(!lower_condition){
        std::cout << "do not stop as current node has to improve lowerbound " << std::endl;
        std::cout << "   lowerbound_bc_ " << lowerbound_bc_ << std::endl;
        std::cout << "   approx_node/(1.+err) " << approx_node/(1.+err) << std::endl;
        std::cout << "   approx_node " << approx_node << std::endl;
        std::cout << "   node_cluster_index " << node_cluster_index << std::endl;
      }
      if(!upper_condition){
        std::cout << "do not stop as current node has to improve upperbound " << std::endl;
        std::cout << "   upperbound_bc_ " << upperbound_bc_ << std::endl;
        std::cout << "   approx_node/(1.-err) " << approx_node/(1.-err) << std::endl;
        std::cout << "   approx_node " << approx_node << std::endl;
        std::cout << "   node_cluster_index " << node_cluster_index << std::endl;
      }
      is_relative_approx = is_relative_approx && lower_condition && upper_condition;
    }

  }

  std::cout << "*** ub_topk_bc " << ub_topk_bc << std::endl;
  std::cout << "*** lb_topk_bc " << lb_topk_bc << std::endl;

  delete(sorted_ub_topk);
  delete(sorted_lb_topk);

  if(is_relative_approx){
    std::cout << "RELATIVE APPROX FOR TOP-K OBTAINED" << std::endl;
    for(uint32_t i=0; i<num_inserted; i++ ){
      uint64_t node_id = this->top_k->get(i);
      double approx_node = this->top_k->get_value(i)/num_samples_d;
      int node_cluster_index = highfreq_set[node_id];
      double eps_current_node = eps_clusters[node_cluster_index];
      double lowerbound_bc_ = approx_node-eps_current_node;
      double upperbound_bc_ = approx_node+eps_current_node;
      if(upperbound_bc_ > lb_topk_bc){
        std::cout << i+1 <<")  ("<<node_id<<") " << lowerbound_bc_ << " " << approx_node << " " << upperbound_bc_ << " " << std::endl;
        eps_final_topk = eps_current_node;
      }
      /*else{
        std::cout << "* " << i+1 <<")  ("<<node_id<<") " << lowerbound_bc_ << " " << approx_node << " " << upperbound_bc_ << " " << std::endl;
      }*/
    }
  }

  free(eps_clusters);

  return is_relative_approx;

}

double Probabilistic::get_max_epsilon_mcrade(double maxestbc_high , double maxestbc_low  , double mcrade_low , bool verbose , bool last_iter) const{

  double delta_for_progressive_bound = delta/pow(2.,iteration_index);
  double delta_for_each_high_freq_f = delta_for_progressive_bound/(2.*num_freq);
  double delta_for_mcrade_low_freq_f = delta_for_progressive_bound/10.;
  double num_samples_d = (double)num_samples;

  // first, compute bound for each function in the high freq clusters
  double eps_high_freq = 0.;
  if(num_freq > 0){
    double log_term_high = log(1./delta_for_each_high_freq_f)/num_samples_d;
    double max_bc_high = maxestbc_high + log_term_high + sqrt( pow(log_term_high,2) + 2.*maxestbc_high*log_term_high  );
    double var_bound_high = (max_bc_high > 0.5) ? 0.25 : max_bc_high*(1.-max_bc_high);
    eps_high_freq = sqrt( 2.*(var_bound_high*log_term_high ) );
  }

  // now, bound on the other functions using mcrade
  double log_term_mcrade = log(1./delta_for_mcrade_low_freq_f)/num_samples_d;
  double max_bc_low = maxestbc_low + log_term_mcrade + sqrt( pow(log_term_mcrade,2) + 2.*maxestbc_low*log_term_mcrade  );
  double var_bound_low = (max_bc_low > 0.5) ? 0.25 : max_bc_low*(1.-max_bc_low);
  double log_term_mcrade_n = log_term_mcrade/(double)mctrials;
  double era_low = mcrade_low + 2.*log_term_mcrade_n + sqrt( pow(2.*log_term_mcrade_n,2) + 4.*(mcrade_low+max_bc_low)*log_term_mcrade_n  );
  double rc_low = era_low + log_term_mcrade + sqrt( pow(log_term_mcrade,2) + 2.*era_low*log_term_mcrade );
  double eps_low_mcrade = 2.*rc_low + sqrt( 2.*log_term_mcrade*(var_bound_low+4.*rc_low) ) + log_term_mcrade/3.;

  double sup_eps = max(eps_low_mcrade , eps_high_freq);

  if(verbose){
    std::cout << "   iteration_index: " << iteration_index << std::endl;
    std::cout << "   delta_for_progressive_bound: " << delta_for_progressive_bound << std::endl;
    std::cout << "   delta_for_each_high_freq_f: " << delta_for_each_high_freq_f << std::endl;
    std::cout << "   delta_for_mcrade_low_freq_f: " << delta_for_mcrade_low_freq_f << std::endl;
    std::cout << "   maxestbc_high: " << maxestbc_high << std::endl;
    std::cout << "   maxestbc_low: " << maxestbc_low << std::endl;
    std::cout << "   max_bc_low: " << max_bc_low << std::endl;
    std::cout << "   era_low: " << era_low << std::endl;
    std::cout << "   rc_low: " << rc_low << std::endl;
    std::cout << "   mcrade_low: " << mcrade_low << std::endl;
    std::cout << "   eps_high_freq: " << eps_high_freq << std::endl;
    std::cout << "   eps_low_mcrade: " << eps_low_mcrade << std::endl;
    std::cout << "***PS sup_eps: " << sup_eps << std::endl;
    std::cout << std::endl;
  }

  return sup_eps;

}


// Runs the algorithm.
// INPUT: k is the number of betweenness that have to be approximated (if k=0 all betweenness
// are approximated with absolute error); delta is the probabilistic guarantee; err is the
// maximum error allowed; union_sample and start_factor are parameters of the algorithm
// that are automatically chosen.
void Probabilistic::run(uint32_t k, double delta, double err, uint32_t union_sample, uint32_t start_factor) {
    this->absolute = (k == 0);
    this->err = err;
    this->delta = delta;
    this->start_factor = start_factor;
    int graph_diameter = estimate_diameter();
    std::cout << "estimated diameter of the graph: " << graph_diameter << std::endl;
    this->omega = 0.5 / err / err * (log2(graph_diameter-1) + 1 + log(2. / delta));
    // here it should be log(2 / delta) ?
    std::cout << "this->omega: " << this->omega << std::endl;
    uint32_t tau = omega / start_factor; // Da sistemare
    std::cout << "tau: " << tau << std::endl;
    tau = max(tau,(uint32_t)1000);
    tau = min(tau,(uint32_t)100000);
    double rel_param = 0.4;
    double test_theta = 0.01;
    double test_rel = 2. / (rel_param*rel_param*test_theta) * (log2(2*(graph_diameter-1))*log(1.0 / test_theta) + log(20.0 / delta) );
    std::cout << "test_rel: " << test_rel << std::endl;
    num_samples = 0;
    sup_bcest = 0;


    if (union_sample == 0) {
        union_sample = min(get_nn(), (uint32_t) max( 2 * sqrt(get_ne()) / omp_get_max_threads(), k+20. ));
    }
    this->union_sample=union_sample;
    this->k=min(k, get_nn());

    last_output = get_time_sec();
    start_time = get_time_sec();
    this->top_k = new Ranking_list(union_sample);
    srand( SEED );
    uint32_t *random_seed = (uint32_t *) malloc( omp_get_max_threads()*sizeof(uint32_t) );
    for( int i=0; i < omp_get_max_threads(); i++ ){
        random_seed[i] = rand();
    }
    distinct_nodes_top_k = 0;
    double high_freq_thr_norm = 1.0;
    double highest_freq = 1.0;

    if(absolute){

      #pragma omp parallel
      {
          Sp_sampler sp_sampler( this, random_seed[omp_get_thread_num()] );
          while (n_pairs <= tau) {
              one_round(sp_sampler);
              double current_time = get_time_sec();

              if (verbose > 0 && current_time - last_output > verbose) {
                  #pragma omp critical
                  {
                      if (current_time - last_output > verbose) {
                          last_output = current_time;
                          cout << "First visits: " << n_pairs << "/" << tau << ".\n";
                      }
                  }
              }
          }
      }

    }
    else{

      int samples_per_step = 10;
      bool stop_first_pass = false;
      #pragma omp parallel
      {
          Sp_sampler sp_sampler( this, random_seed[omp_get_thread_num()] );
          while( !stop_first_pass ) {
              for (int i = 0; i <= samples_per_step; i++) {
                  one_round(sp_sampler);
              }
              if( num_samples > 10000 && !stop_first_pass ){
                #pragma omp critical
                {
                    // check if k distinct nodes have been sampled
                    stop_first_pass = distinct_nodes_top_k >= k;
                }
              }
          }
      }

      theta_top_k = max(1. , top_k->get_value(k-1))/(double)num_samples;
      double theta_top_1 = (double)top_k->get_value(0)/(double)num_samples;
      highest_freq = theta_top_1;
      double log_base = 1.5;
      high_freq_thr_norm = log_base*theta_top_k;
      num_clusters_top_k = log(theta_top_1/theta_top_k)/log(log_base)+2;
      for (uint32_t i = 0; i < get_nn(); i++) {
          highfreq_set[i] = 0;
      }
      std::cout << "stopped after  " << num_samples << " iteration" << std::endl;
      std::cout << "theta_top_k " << theta_top_k << std::endl;
      for(uint64_t i = 0; i < k; i++){
        uint64_t node_index = top_k->get(i);
        double node_freq = top_k->get_value(i)/(double)num_samples;
        int node_cluster_index = log(node_freq/theta_top_k)/log(log_base)+1;
        highfreq_set[node_index] = node_cluster_index;
        std::cout << i+1 << "-th node: \n   node_freq " << node_freq;
        std::cout << "\n   node_cluster_index " << node_cluster_index << std::endl;
      }

      sup_bcest_cluster = (uint64_t*) calloc( num_clusters_top_k , sizeof(uint64_t) );
      max_mcera_cluster = (int64_t*) calloc( mctrials*num_clusters_top_k , sizeof(int64_t) );
      for(int i = 0; i < num_clusters_top_k; i++){
        sup_bcest_cluster[i] = 0;
      }
      std::cout << "num_clusters_top_k  " << num_clusters_top_k << std::endl;
    }

    *time_bfs = 0;
    *time_critical = 0;

    std::cout << "time for first pass " << get_time_sec() - start_time << std::endl;

    if(absolute){
      int k_first_pass = graph_diameter-1;
      Ranking_list* top_k_first_pass = new Ranking_list(k_first_pass);
      for (uint32_t i = 0; i < get_nn(); i++) {
        top_k_first_pass->put(i, approx[i]);
      }
      int high_freq_elements = 0;
      double high_freq_thr = top_k_first_pass->get_value(k_first_pass-1);
      high_freq_thr_norm = high_freq_thr / (double)tau;
      highest_freq = 10.+top_k_first_pass->get_value(0);
      highest_freq = highest_freq / (double)tau;
      high_freq_thr = max(high_freq_thr , 10.);
      for (uint32_t i = 0; i < get_nn(); i++) {
          highfreq_set[i] = approx[i] >= high_freq_thr;
          high_freq_elements += highfreq_set[i];
      }
      cout << "high_freq_thr: " << high_freq_thr << std::endl;
      cout << "high_freq_elements: " << high_freq_elements << std::endl;
      num_freq = high_freq_elements;
      cout << "highest_freq: " << highest_freq << std::endl;
      cout << "high_freq_thr_norm: " << high_freq_thr_norm << std::endl;
      delete(top_k_first_pass);
    }

    double start_time_delta_guess = get_time_sec();

    compute_delta_guess();

    std::cout << "time for compute delta guess " << get_time_sec() - start_time_delta_guess << std::endl;

    n_pairs = 0;
    delete(this->top_k);
    this->top_k = new Ranking_list(union_sample);
    uint32_t v_mc_index;
    for (uint32_t i = 0; i < get_nn(); i++) {
        approx[i] = 0;
        // mcrade
        v_mc_index = i*mctrials;
        for(uint32_t j = 0; j < mctrials; j++){
          mcrade[v_mc_index+j] = 0;
        }
    }
    sup_mcera = 0;
    sup_bcest = 0;
    num_samples = 0;
    if(!absolute){
      omega = pow(10.,15);
    }
    // guess a first stopping condition check according to what we computed in the first phase
    // of the adaptive algorithm
    iteration_index = 1;
    double first_stopping_samples = 0;//2./err/err*( highest_freq*log(2./delta) );
    double eps_guess = 1.;
    // start a binary search to find a good starting sample size
    double first_sample_lower = 1./err*log(2./delta);
    double first_sample_upper = omega;
    while(first_sample_upper - first_sample_lower > 10){
      cout << "--- LOG SEARCH ITER" << std::endl;
      num_samples = (first_sample_upper + first_sample_lower)/2.;
      double guess_mcera = (1+sqrt(high_freq_thr_norm*num_samples))/num_samples;
      guess_mcera = min(guess_mcera , 0.9*high_freq_thr_norm);
      eps_guess = get_max_epsilon_mcrade(highest_freq , high_freq_thr_norm , guess_mcera, true, false);
      if(absolute){
        if(eps_guess > err){
          first_sample_lower = num_samples;
          cout << "lower increased at: " << first_sample_lower << std::endl;
        }
        else{
          first_sample_upper = num_samples;
          cout << "upper decreased at: " << first_sample_upper << std::endl;
        }
      }
      else{
        bool violated_condition_lower = theta_top_k+eps_guess > theta_top_k/(1.-err);
        bool violated_condition_upper = theta_top_k-eps_guess < theta_top_k/(1.+err);
        if(violated_condition_lower || violated_condition_upper){
          first_sample_lower = num_samples;
          cout << "(top-k) lower increased at: " << first_sample_lower << std::endl;
        }
        else{
          first_sample_upper = num_samples;
          cout << "(top-k) upper decreased at: " << first_sample_upper << std::endl;
        }
      }
      cout << "  eps_guess is at: " << eps_guess << std::endl;
    }
    first_stopping_samples = num_samples;
    num_samples = 0;


    cout << "first_stopping_samples: " << first_stopping_samples << std::endl;
    // initialize the next at the first, it will be updated during iterations
    next_stopping_samples = first_stopping_samples;
    // we stop at omega
    double last_stopping_samples = omega;
    cout << "last_stopping_samples: " << last_stopping_samples << std::endl;
    if(first_stopping_samples >= last_stopping_samples){
      first_stopping_samples = last_stopping_samples/8;
      cout << "first_stopping_samples dropped to " << first_stopping_samples << std::endl;
    }
    bool stop_mcrade = false;

    #pragma omp parallel
    {
        Sp_sampler sp_sampler( this, random_seed[omp_get_thread_num()] );
        Status status(union_sample);
        status.n_pairs = 0;
        bool stop = false;

        while( !stop && status.n_pairs < omega ) {
            for (uint32_t i = 0; i <= 10; i++) {
                one_round(sp_sampler);
            }
            //get_status (&status);
            if(absolute){
              get_status (&status);
              time_comp_finished[omp_get_thread_num()] -= get_time_sec();
              stop = compute_finished(&status);
              time_comp_finished[omp_get_thread_num()] += get_time_sec();
            }

            // check stopping condition for mcrade
            if (!stop_mcrade && num_samples < last_stopping_samples && num_samples >= next_stopping_samples) {
                get_status (&status);
                #pragma omp critical(stopcond)
                {
                    if (num_samples < last_stopping_samples && num_samples >= next_stopping_samples) {
                        stop_mcrade = compute_finished_mcrade(&status , false);
                        next_stopping_samples = get_next_stopping_sample();
                        if(stop_mcrade){
                          std::cout << "/* stop_mcrade is true */" << std::endl;
                          if(!absolute){
                            stop = stop_mcrade;
                            omega = status.n_pairs;
                          }
                        }
                    }
                }
            }

            double current_time = get_time_sec();

            if(!absolute && stop_mcrade){
              stop = stop_mcrade;
            }

            if (verbose > 0 && current_time - last_output > verbose) {
                get_status (&status);
                #pragma omp critical(print)
                {
                    if (current_time - last_output > verbose) {
                        last_output = current_time;
                        print_status(&status);
                    }
                }
            }
        }
    }
    cout << "out of second pass " << std::endl;
    if (/*absolute &&*/ verbose > 0) {
        if(!absolute){
          omega = 0.5 / eps_final_topk / eps_final_topk * (log2(graph_diameter-1) + 1 + log(2. / delta));
          cout << "OMEGA WITH PRIOR KNOWLEDGE OF EPS " << omega << std::endl;
        }
        Status status(union_sample);
        get_status(&status);
        print_status(&status, true);
        // compute last bound, using all delta budget and VC dimension bound
        if(!stop_mcrade)
          compute_finished_mcrade(&status , true);
    }
    free(random_seed);
    n_pairs += tau;
    cout << "finished run " << std::endl;
}

// Destructor of the class Probabilistic.
Probabilistic::~Probabilistic() {
    free(approx);
    free(delta_l_guess);
    free(delta_u_guess);
    delete(top_k);
    // mcrade
    free(mcrade);
    delete(mcrade_randgen);
    free(highfreq_set);
    free(max_mcera_exact);
    free(max_mcera_exact_low);
    free(max_mcera_exact_high);
    free(sup_bcest_cluster);
    free(max_mcera_cluster);
}
