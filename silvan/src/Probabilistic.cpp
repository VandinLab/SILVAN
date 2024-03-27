#include <string>
#include <iostream>
#include <math.h>
#include <limits.h>
#include <iomanip>
#include <time.h>
#include <cfloat>
#include <omp.h>
#include <fstream>

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
Probabilistic::Probabilistic( const std::string &filename, const bool directed, const double verb, const double sampling_rate_, const bool alpha_given_, const double empirical_peeling_param_ , const bool enable_m_hat_, const std::string output_file_ ): Graph( filename, directed ), verbose(verb) {
    approx = (double *) calloc( get_nn(), sizeof(double) );
    approx_toadd = (double *) calloc( get_nn(), sizeof(double) );
    time_bfs = (double *) calloc( omp_get_max_threads(), sizeof(double) );
    time_comp_finished = (double *) calloc( omp_get_max_threads(), sizeof(double) );
    time_critical = (double *) calloc( omp_get_max_threads(), sizeof(double) );
    time_critical_round = (double *) calloc( omp_get_max_threads(), sizeof(double) );
    time_mcera = (double *) calloc( omp_get_max_threads(), sizeof(double) );
    emp_wimpy_vars = (double *) calloc( get_nn(), sizeof(double) );
    for (uint32_t i = 0; i < get_nn(); i++) {
        approx[i] = 0;
        approx_toadd[i] = 0;
        emp_wimpy_vars[i] = 0;
    }
    for(int i=0; i<omp_get_max_threads(); i++){
      time_bfs[i] = 0.;
      time_comp_finished[i] = 0.;
      time_critical[i] = 0.;
      time_critical_round[i] = 0.;
      time_mcera[i] = 0.;
    }
    n_pairs = 0;
    vis_edges = 0;
    if (verbose > 0) {
        print_data();
    }
    output_file = output_file_;
    // mcrade
    mcrade = (int64_t *) calloc( get_nn()*mctrials, sizeof(int64_t) );
    partition_index = (int *) calloc( get_nn(), sizeof(int) );
    mcrade_randgen = new Rand_gen( 2021 );
    sup_bcest = 0.0;
    alpha_sp_given = alpha_given_;
    alpha_sp_sampling = sampling_rate_;
    empirical_peeling_a = empirical_peeling_param_;
    enable_m_hat = enable_m_hat_;
}


double numsamplesboundfunction(double x_ , double rho, double delta_ , double eps_){
  double v_x = x_*(1.-x_);
  double arg_h = eps_/v_x;
  double denom = v_x * ( (1+arg_h) * log(1+arg_h) - arg_h  );
  return log( 2*rho / ( x_ * delta_ ) ) / ( denom );
}

double getUpperBoundSamples(double max_bc , double max_var , double avg_spl , double eps , double delta_bound){

  // computation of x_hat
  double x_hat = 0.;
  double x_hat_low = 0.5 - sqrt(eps/3.);
  x_hat_low = max(x_hat_low , 0.);
  double x_hat_high = 0.5;
  bool debug_getUpperBoundSamples = false;
  if(debug_getUpperBoundSamples){
    std::cout << "Call to getUpperBoundSamples " << std::endl;
    std::cout << "   max_bc " << max_bc << std::endl;
    std::cout << "   max_var " << max_var << std::endl;
    std::cout << "   avg_spl " << avg_spl << std::endl;
    std::cout << "   eps " << eps << std::endl;
    std::cout << "   delta_bound " << delta_bound << std::endl;
    std::cout << "   x_hat_low " << x_hat_low << std::endl;
    std::cout << "   x_hat_high " << x_hat_high << std::endl;
  }

  while(x_hat_high - x_hat_low > 0.0001){
    x_hat = (x_hat_high + x_hat_low)/2.;
    if(debug_getUpperBoundSamples)
      std::cout << "   x_hat " << x_hat << std::endl;
    double v_x = x_hat*(1.-x_hat);
    double arg_h = eps/v_x;
    double f_val = v_x * ( (1+arg_h) * log(1+arg_h) - arg_h  );
    if ( f_val <= 2*pow(eps,2) ){
      x_hat_high = x_hat;
    }
    else{
      x_hat_low = x_hat;
    }
  }
  x_hat = x_hat_high;
  if(debug_getUpperBoundSamples)
    std::cout << "final x_hat " << x_hat << std::endl;

  {
    double v_x = x_hat*(1.-x_hat);
    double arg_h = eps/v_x;
    double f_val = v_x * ( (1+arg_h) * log(1+arg_h) - arg_h  );
    bool check = f_val <= 2*pow(eps,2);
    if(debug_getUpperBoundSamples)
      std::cout << "logic check " << check << std::endl;
  }

  // compute the bound on the number of samples
  double x_high = min(x_hat , max_bc);
  double x_hat_var = 1.0;
  if(max_var < 0.25){
    x_hat_var = 0.5 - sqrt(0.25 - max_var);
    x_high = min(x_high , x_hat_var);
  }

  double x_low = x_high;
  double step = x_high/1000.;

  double num_samples_bound_high = numsamplesboundfunction(x_high , avg_spl , delta_bound , eps);
  if(debug_getUpperBoundSamples){
    std::cout << "num_samples_bound_high " << num_samples_bound_high << std::endl;
    std::cout << "x_high " << x_high << std::endl;
  }
  double num_samples_bound_ = num_samples_bound_high+1;
  while(num_samples_bound_ >= num_samples_bound_high){
    x_low = x_high - step;
    if (x_low > 0.){
      num_samples_bound_ = numsamplesboundfunction(x_low , avg_spl , delta_bound , eps);
      if(num_samples_bound_ > num_samples_bound_high){
        x_high = x_low;
        num_samples_bound_high = num_samples_bound_;
      }
    }
    else{
      num_samples_bound_ = num_samples_bound_high-1;
    }
    if(debug_getUpperBoundSamples){
      std::cout << "    num_samples_bound_ " << num_samples_bound_ << std::endl;
      std::cout << "    x_low " << x_low << std::endl;
    }
  }
  if(debug_getUpperBoundSamples){
    std::cout << "  *** num_samples_bound_ " << num_samples_bound_ << std::endl;
    std::cout << "  *** x_low " << x_low << std::endl;
  }

  return num_samples_bound_high;

}

//#define DEBUG_REL_BOUNDS 1


double Probabilistic::computeRelBound(double bc_est , double delta , double rho , double m_samples , bool upper){

  #ifdef DEBUG_REL_BOUNDS
  std::cout << " Call to computeRelBound " << std::endl;
  std::cout << "    bc_est: " << bc_est << std::endl;
  std::cout << "    delta: " << delta << std::endl;
  std::cout << "    rho: " << rho << std::endl;
  std::cout << "    m_samples: " << m_samples << std::endl;
  std::cout << "    upper: " << upper << std::endl;
  #endif

  double quant_condition = 0.;
  double upper_ = 1.0;
  double lower_ = 0.;
  if(upper){
    lower_ = bc_est;
  }
  else{
    upper_ = bc_est;
  }
  double test_ = 1.0;
  int i = 0;
  while(i < 30){
    test_ = (upper_+lower_)/2.;
    double v_test_ = test_*(1.-test_);
    double log_term_ = log( 2./delta * min((double)get_nn() , rho/test_) )/m_samples;
    double sqrt_term = sqrt(2.*log_term_*v_test_);
    if(upper){
      quant_condition = bc_est + log_term_/3. + sqrt_term;
      if(quant_condition >= test_){
        lower_ = test_;
      }
      else{
        upper_ = test_;
      }
    }
    else{
      quant_condition = bc_est - log_term_/3. - sqrt_term;
      if(quant_condition <= test_){
        upper_ = test_;
      }
      else{
        lower_ = test_;
      }
    }
    #ifdef DEBUG_REL_BOUNDS
    std::cout << " iteration i: " << i << std::endl;
    std::cout << "    quant_condition: " << quant_condition << std::endl;
    std::cout << "    upper_: " << upper_ << std::endl;
    std::cout << "    lower_: " << lower_ << std::endl;
    #endif
    i++;
  }
  if(upper){
    #ifdef DEBUG_REL_BOUNDS
    std::cout << " Computed upper bound: " << upper_ << std::endl;
    #endif
    return upper_;
  }
  else{
    #ifdef DEBUG_REL_BOUNDS
    std::cout << " Computed lower bound: " << lower_ << std::endl;
    #endif
    return lower_;
  }

}


// Decides whether the algorithm should terminate
// according to mcera
// INPUT: a Status object describing the current status of the algorithm
bool Probabilistic::compute_finished_mcrade(Status *status) {

  std::cout << std::endl;
  double num_samples_d = (double)num_samples;
  double delta_for_progressive_bound_ = 1.0;
  bool debug_eps_partitions = false;
  if(absolute && enable_m_hat){
    delta_for_progressive_bound_ = delta_for_progressive_bound;
  }
  else{
    delta_for_progressive_bound_ = delta/pow(2.,iteration_index);
  }
  cout << "EVALUATING STOPPING CONDITION at iteration " << iteration_index << " (m="<< num_samples <<", d="<<delta_for_progressive_bound_<<")" << std::endl;
  //std::cout << "(Printing thread: " << omp_get_thread_num() << ")" << endl;

  if(absolute && second_phase_started){
    // first, check updated upper bound on the samples
    // compute upper bound on average shortest path length
    double avg_diam_upperbound = getUpperBoundAvgDiameter(delta_for_progressive_bound_ , false);

    // compute upper bound on top-1 bc
    double top1_est_bc = approx[status->top_k[0]]/ num_samples_d;
    double top1bc_upperbound = getUpperBoundTop1BC(top1_est_bc , delta_for_progressive_bound_);
    double wimpy_var_upperbound = getUpperBoundTop1BC(sup_emp_wimpy_var/num_samples_d , delta_for_progressive_bound_);

    // compute upper limit on number of samples
    if(enable_m_hat){
      uint64_t max_num_samples = getUpperBoundSamples(top1bc_upperbound , wimpy_var_upperbound , avg_diam_upperbound , err , delta_for_progressive_bound_);
      if(last_stopping_samples > (double)max_num_samples){
        last_stopping_samples = (double)max_num_samples;
        omega = last_stopping_samples;
        cout << "******NEW STOPPING CONDITION UPDATE!!! last_stopping_samples: " << last_stopping_samples << std::endl;
        double max_num_iterations = log((double)max_num_samples/(double)first_stopping_samples)/log(1.2)+1;
        delta_for_progressive_bound = delta/(floor(max_num_iterations));
      }
      if(max_num_samples <= num_samples){
        cout << "******NEW STOPPING CONDITION TRUE!!! " << std::endl;
      }
    }
  }

  // first, reset counters for computation of mcera
  std::vector<int> counts_nodes_partition;
  std::vector<double> mcera_full_f_trials;
  counts_nodes_partition.resize(number_of_non_empty_partitions);
  mcera_full_f_trials.resize(mctrials);
  for(int i = 0; i < number_of_non_empty_partitions; i++){
    counts_nodes_partition[i] = 0;
    sup_bcest_partition[i] = 0;
    sup_empwvar_partition[i] = 0;
    epsilon_partition[i] = 1.0;
    int mctrials_idx = i*mctrials;
    for(uint32_t j = 0; j < mctrials; j++){
      max_mcera_partition[mctrials_idx+j] = -(int64_t)status->n_pairs;
    }
  }

  // iterate over nodes to update the mcera
  for (uint32_t i = 0; i < get_nn(); i++) {
      uint32_t v = i;
      uint32_t v_rade_idx = v*mctrials;

      int node_partition_index = partition_index[v];
      int mapped_partition_index = partitions_ids_map[node_partition_index];
      sup_bcest_partition[mapped_partition_index] = max(sup_bcest_partition[mapped_partition_index] , (uint64_t)approx[v]);
      sup_empwvar_partition[mapped_partition_index] = max(sup_empwvar_partition[mapped_partition_index] , emp_wimpy_vars[v]);
      counts_nodes_partition[mapped_partition_index] = counts_nodes_partition[mapped_partition_index]+1;
      int mcera_partition_index = mctrials*mapped_partition_index;
      for(uint32_t j = 0; j < mctrials; j++){
        max_mcera_partition[mcera_partition_index+j] = max(max_mcera_partition[mcera_partition_index+j] , mcrade[v_rade_idx+j]);
      }
  }
  double sup_empwvar_full_f = 0.;
  for(int i = 0; i < number_of_non_empty_partitions; i++){
    sup_empwvar_full_f = max(sup_empwvar_full_f , (double)sup_empwvar_partition[i]);
  }
  sup_empwvar_full_f = sup_empwvar_full_f/num_samples_d;
  double mcera_full_f = 0.;
  for(uint32_t j = 0; j < mctrials; j++){
    mcera_full_f_trials[j] = -(double)status->n_pairs;
    for(int i = 0; i < number_of_non_empty_partitions; i++){
      mcera_full_f_trials[j] = max(mcera_full_f_trials[j] , (double)max_mcera_partition[mctrials*i+j]);
    }
    mcera_full_f += mcera_full_f_trials[j]/(double)mctrials;
  }
  mcera_full_f = mcera_full_f/num_samples_d;
  double eps_full_f = get_epsilon_mcrade(sup_empwvar_full_f , mcera_full_f , delta_for_progressive_bound_ , num_samples_d , -1 , false);

  // compute mcera for all partitions
  // this array containes the mcera for every partition
  std::vector<double> mcera_partition_avg;
  mcera_partition_avg.resize(number_of_non_empty_partitions);
  double mcera_avg_ = 0.0;
  int mcera_partition_index = 0;
  for(int i = 0; i < number_of_non_empty_partitions; i++){
    mcera_avg_ = 0.0;
    mcera_partition_index = mctrials*i;
    for(uint32_t j = 0; j < mctrials; j++){
      mcera_avg_ += max_mcera_partition[mcera_partition_index+j]/(double)mctrials;
    }
    mcera_avg_ = mcera_avg_/num_samples_d;
    mcera_partition_avg[i] = mcera_avg_;
    double sup_emp_wimpy_var_ = sup_empwvar_partition[i]/num_samples_d;
    // compute epsilon for this partition
    double delta_this_partition = delta_for_progressive_bound_/number_of_non_empty_partitions;
    //double delta_this_partition = delta_for_progressive_bound_/(pow(2.,(double)i+1.));
    double current_eps = get_epsilon_mcrade(sup_emp_wimpy_var_ , mcera_avg_ , delta_this_partition , num_samples_d , counts_nodes_partition[i] , false);
    epsilon_partition[i] = current_eps;
  }

  // the mcera and epsilons has been computed , now evaluate the stopping condition
  // stopping condition for absolute (uniform) approximation
  // stop when all supremum deviations are <= than the desired epsilon
  if(absolute){

    double sup_eps = 0;
    for(int i = 0; i < number_of_non_empty_partitions; i++){
      sup_eps = max(sup_eps , epsilon_partition[i]);
    }

    if(sup_eps <= err){
      cout << "MCRADE STOPS WITH eps " << sup_eps << std::endl;
      cout << "MCRADE STOPS at iteration " << iteration_index << std::endl;
      cout << "MCRADE STOPS at sample size " << num_samples << std::endl;
      cout << "MCRADE STOPS after seconds " << get_time_sec() - start_time << std::endl;
    }
    else{
      cout << "   mcrade eps " << sup_eps << std::endl;
      if(debug_eps_partitions){
        cout << "      eps_full_f " << eps_full_f << std::endl;
        cout << "         supvar_full_f " << sup_empwvar_full_f << std::endl;
        cout << "         mcera_full_f " << mcera_full_f << std::endl;
        for(int i = 0; i < number_of_non_empty_partitions; i++){
          cout << "      epsilon_partition["<< i <<"] " << epsilon_partition[i] << std::endl;
          cout << "         supvar["<< i <<"] " << sup_empwvar_partition[i]/num_samples_d << std::endl;
          cout << "         mcera_partition_avg["<< i <<"] " << mcera_partition_avg[i] << std::endl;
          cout << "         counts_nodes_partition["<< i <<"] " << counts_nodes_partition[i] << std::endl;
        }
      }
    }

    return sup_eps <= err;
  }
  // stopping condition for top-k approximation
  else{

    // check if the guarantees on the top-k are satisfied
    bool top_k_done = check_topk_guarantees(true);
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

    return top_k_done;

  }

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
    if (full) {
      std::cout << "Time bfs: " << time_bfs[omp_get_thread_num()] << endl;
      std::cout << "Time critical: " << time_critical[omp_get_thread_num()] << endl;
      std::cout << "Time compute finished: " << time_comp_finished[omp_get_thread_num()] << endl;
      std::cout << "Time mcera: " << time_mcera[omp_get_thread_num()] << endl;
      std::cout << "(Printing thread: " << omp_get_thread_num() << ")" << endl;
    }

    // compute averages time over threads
    double avg_time_mcera = 0.;
    double avg_time_critical = 0.;
    double avg_time_bfs = 0.;
    double avg_time_finished = 0.;
    double avg_time_critical_round = 0.;
    for(int i = 0; i < omp_get_max_threads(); i++){
      avg_time_bfs += time_bfs[i]/(double)omp_get_max_threads();
      avg_time_critical += time_critical[i]/(double)omp_get_max_threads();
      avg_time_critical_round += time_critical_round[i]/(double)omp_get_max_threads();
      avg_time_mcera += time_mcera[i]/(double)omp_get_max_threads();
      avg_time_finished += time_comp_finished[i]/(double)omp_get_max_threads();
    }
    if (full) {
      std::cout << "Avg time bfs: " << avg_time_bfs << endl;
      std::cout << "Avg time critical: " << avg_time_critical << endl;
      std::cout << "Avg time critical (round): " << avg_time_critical_round << endl;
      std::cout << "Avg time compute finished: " << avg_time_finished << endl;
      std::cout << "Avg time mcera: " << avg_time_mcera << endl;

      // mcera
      std::cout << "sup bc est: " << sup_bcest/(double)num_samples << endl;
      std::cout << "sup_emp_wimpy_var: " << sup_emp_wimpy_var/(double)num_samples << std::endl;
      std::cout << "void_samples: " << void_samples << endl;
      std::cout << "num_samples: " << num_samples << endl;
    }

}

//#define ROUND_DEBUG 1

// Sample one shortest path and updates the ranking of the betweenness approximations.
void Probabilistic::one_round(Sp_sampler &sp_sampler) {
    int path_length = 0;
    int num_paths = 0;
    time_bfs[omp_get_thread_num()] -= get_time_sec();
    map<uint32_t, int>/*vector<uint32_t>*/ path = sp_sampler.random_path(path_length , num_paths , alpha_sp_sampling);
    time_bfs[omp_get_thread_num()] += get_time_sec();


    time_critical[omp_get_thread_num()] -= get_time_sec();

    if(path_length > 0){
      #pragma omp critical
      {
          n_pairs++;
          vis_edges += sp_sampler.vis_edges;
          num_samples++;

          // mcrade
          int* sigmas = 0;
          if(path_length > 0 && firstpass == false){
            uint64_t maxval_sigmas = 100000000;
            double maxval_half = (double)maxval_sigmas/2.;
            time_mcera[omp_get_thread_num()] -= get_time_sec();
            sigmas = (int*) calloc( mctrials , sizeof(int));
            for(uint32_t j = 0; j < mctrials; j++){
              sigmas[j] = (mcrade_randgen->get_max(maxval_sigmas) >= maxval_half)*2-1;
              //std::cout << "sigmas[j] " << sigmas[j] << std::endl;
            }
            time_mcera[omp_get_thread_num()] += get_time_sec();
          }
          double one_over_num_paths = 1./(double)num_paths;
          /*for(uint32_t u:path){
            approx_toadd[u] += one_over_num_paths;
          }*/
          for (const auto& elem_path : path) {
          //for(uint32_t u:path){
              //approx[u]++;
              uint32_t u = elem_path.first;
              double appx_to_add_u = elem_path.second*one_over_num_paths;//approx_toadd[u];
              if(appx_to_add_u > 0.){
                //approx_toadd[u] = 0.;
                approx[u] += appx_to_add_u;
                emp_wimpy_vars[u] += pow(appx_to_add_u,2.);
                top_k->put(u, approx[u]);
                if(firstpass == false){
                  // mcrade
                  time_mcera[omp_get_thread_num()] -= get_time_sec();
                  uint32_t u_idx = u*mctrials;
                  for(uint32_t j = 0; j < mctrials; j++){
                    mcrade[u_idx+j] += sigmas[j]*appx_to_add_u;
                  }
                  time_mcera[omp_get_thread_num()] += get_time_sec();
                }
                if(firstpass == true && !absolute){
                  if(approx_toadd[u] == 0 && approx[u] >= 3){
                    distinct_nodes_top_k += 1;
                    approx_toadd[u] = 1;
                  }
                }
              if(approx[u] > sup_bcest){
                sup_bcest = approx[u];
              }
              if(emp_wimpy_vars[u] > sup_emp_wimpy_var){
                sup_emp_wimpy_var = emp_wimpy_vars[u];
              }
            }

          }
          void_samples += path_length == 0;
          sp_lengths[path_length] += 1;

          if(path_length > 0 && firstpass == false){
            time_mcera[omp_get_thread_num()] -= get_time_sec();
            free(sigmas);
            time_mcera[omp_get_thread_num()] += get_time_sec();
          }
      }
    }
    else{
      #pragma omp atomic
      n_pairs++;
      #pragma omp atomic
      vis_edges += sp_sampler.vis_edges;
      #pragma omp atomic
      num_samples++;
      #pragma omp atomic
      void_samples ++;
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



// update the stopping condition schedule
double Probabilistic::get_next_stopping_sample(){
  next_stopping_samples = next_stopping_samples*1.2;
  iteration_index += 1;
  return next_stopping_samples;
}

// function to check that the guarantees on the top-k are satisfied
double Probabilistic::check_topk_guarantees(bool verbose){

  // total delta for this iteration, to split across partitions
  double num_samples_d = (double)num_samples;
  std::cout << "EVALUATING STOPPING CONDITION FOR TOP-K RELATIVE APPROX" << std::endl;
  double delta_for_progressive_bound = delta/pow(2.,iteration_index);

  eps_final_topk = 1.0;

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
    int node_partition_index = partition_index[node_id];
    int map_node_partition_index = partitions_ids_map[node_partition_index];
    double eps_current_node = epsilon_partition[map_node_partition_index];
    // we have to check that this eps gives a relative approximation
    // according to the input parameter err
    double lowerbound_bc_ = approx_node-eps_current_node;
    double upperbound_bc_ = approx_node+eps_current_node;
    double avg_diam_upperbound = getUpperBoundAvgDiameter(delta_for_progressive_bound , false);
    double lowerbound_bc_rel = computeRelBound(approx_node , delta_for_progressive_bound , avg_diam_upperbound , num_samples_d , false);
    double upperbound_bc_rel = computeRelBound(approx_node , delta_for_progressive_bound , avg_diam_upperbound , num_samples_d , true);
    lowerbound_bc_ = max(lowerbound_bc_ , lowerbound_bc_rel);
    upperbound_bc_ = min(upperbound_bc_ , upperbound_bc_rel);
    sorted_ub_topk->put(node_id , upperbound_bc_);
    sorted_lb_topk->put(node_id , lowerbound_bc_);
    num_inserted++;
    if (num_inserted >= k){
      ub_topk_bc = min(ub_topk_bc , sorted_ub_topk->get_value(k-1));
      lb_topk_bc = max(lb_topk_bc , sorted_lb_topk->get_value(k-1));
    }
    // this node should go in the output as a top-k node
    if(upperbound_bc_ > lb_topk_bc){
      eps_final_topk = min(eps_final_topk , eps_current_node);
      bool lower_condition = lowerbound_bc_ >= approx_node/(1.+err);
      bool upper_condition = upperbound_bc_ <= approx_node/(1.-err);
      bool lower_condition_rel = lowerbound_bc_rel >= approx_node/(1.+err);
      bool upper_condition_rel = upperbound_bc_rel <= approx_node/(1.-err);
      if(!lower_condition){
        std::cout << "do not stop as current node has to improve lowerbound " << std::endl;
        std::cout << "   lowerbound_bc_ " << lowerbound_bc_ << std::endl;
        std::cout << "   lowerbound_bc_rel " << lowerbound_bc_rel << std::endl;
        std::cout << "   rel condition " << lower_condition_rel << std::endl;
        std::cout << "   approx_node/(1.+err) " << approx_node/(1.+err) << std::endl;
      }
      if(!upper_condition){
        std::cout << "do not stop as current node has to improve upperbound " << std::endl;
        std::cout << "   upperbound_bc_ " << upperbound_bc_ << std::endl;
        std::cout << "   upperbound_bc_rel " << upperbound_bc_rel << std::endl;
        std::cout << "   rel condition " << upper_condition_rel << std::endl;
        std::cout << "   approx_node/(1.-err) " << approx_node/(1.-err) << std::endl;
      }
      if(!lower_condition || !upper_condition){
        std::cout << "   approx_node " << approx_node << std::endl;
        std::cout << "   node_partition_index " << node_partition_index << std::endl;
      }
      is_relative_approx = is_relative_approx && lower_condition && upper_condition;
    }
    else{
      eps_final_topk = min(eps_final_topk , lb_topk_bc-approx_node);
    }

  }

  std::cout << "*** ub_topk_bc " << ub_topk_bc << std::endl;
  std::cout << "*** lb_topk_bc " << lb_topk_bc << std::endl;

  delete(sorted_ub_topk);
  delete(sorted_lb_topk);

  if(is_relative_approx){
    std::cout << "RELATIVE APPROX FOR TOP-K OBTAINED" << std::endl;
    uint32_t num_output_topk = 0;
    for(uint32_t i=0; i<num_inserted; i++ ){
      uint64_t node_id = this->top_k->get(i);
      double approx_node = this->top_k->get_value(i)/num_samples_d;
      int node_partition_index = partition_index[node_id];
      int map_node_partition_index = partitions_ids_map[node_partition_index];
      double eps_current_node = epsilon_partition[map_node_partition_index];
      double lowerbound_bc_ = approx_node-eps_current_node;
      double upperbound_bc_ = approx_node+eps_current_node;
      if(upperbound_bc_ > lb_topk_bc){
        std::cout << i+1 <<") ("<<node_id<<") " << lowerbound_bc_ << " " << approx_node << " " << upperbound_bc_ << " " << std::endl;
        num_output_topk++;
      }
    }
    this->numresults_topk = num_output_topk;
  }

  return is_relative_approx;

}

double Probabilistic::get_epsilon_mcrade(double sup_emp_wimpy_var_ , double mcera_ , double delta_ , double num_samples_ , double num_nodes_ , bool verbose) const{
  // computes the probabilistic upper bound on supremum deviation given the mcera
  mcera_ = max(mcera_ , 0.);
  double log_term_mcrade = log(5./delta_)/num_samples_;
  double var_ub = sup_emp_wimpy_var_ + log_term_mcrade + sqrt( pow(log_term_mcrade,2.) + 2 * log_term_mcrade * sup_emp_wimpy_var_ );
  var_ub = min(var_ub , 0.25);
  double era_ub = mcera_ + sqrt( 4 * sup_emp_wimpy_var_ * log_term_mcrade / (double)mctrials );
  double ra_ub = era_ub + log_term_mcrade + sqrt( pow(log_term_mcrade,2.) + 2 * log_term_mcrade * era_ub );
  if(num_nodes_ == 1){
    ra_ub = 0;
  }
  double eps_ub = 2 * ra_ub + sqrt( 2 * log_term_mcrade * ( var_ub + 4 * ra_ub ) ) + log_term_mcrade/3;
  double eps_ub_2 = 2 * ra_ub + sqrt( log_term_mcrade / 2. );
  eps_ub = min(eps_ub,eps_ub_2);
  if(num_nodes_ >= 1){
    double eps_unb = sqrt( 2 * var_ub * log(5./delta_*num_nodes_)/num_samples_ ) + log(5./delta_*num_nodes_)/num_samples_/3;
    eps_ub = min(eps_ub,eps_unb);
  }
  return eps_ub;
}



double Probabilistic::getUpperBoundTop1BC(double top1_est_bc , double delta){
  bool debug_getUpperBoundTop1BC = false;
  if(debug_getUpperBoundTop1BC){
    std::cout << "Call to getUpperBoundTop1BC " << std::endl;
    std::cout << "    top1_est_bc: " << top1_est_bc << std::endl;
  }
  double log_term_top1bc = log(1./delta);
  double const_term_top1bc = log_term_top1bc/(double)num_samples;
  double top1bc_upperbound = top1_est_bc + const_term_top1bc + sqrt( 2*const_term_top1bc*top1_est_bc + pow(const_term_top1bc,2.) );
  top1bc_upperbound = min(top1bc_upperbound , 1.0);
  if(debug_getUpperBoundTop1BC){
    std::cout << "    top1bc_upperbound: " << top1bc_upperbound << std::endl;
  }
  return top1bc_upperbound;
}


double Probabilistic::getUpperBoundAvgDiameter(double delta , bool verbose){
  // fist get the avg
  double avg_diam_ = 0.;
  for( int i=0; i <= graph_diameter; i++ ){
    avg_diam_ += (double)i*sp_lengths[i]/(double)(num_samples);
  }

  // upper bound using bernstein
  double log_term_avgspl = log(1./delta);
  double const_term_avgspl = (graph_diameter-2)*log_term_avgspl/(double)num_samples;
  double avg_diam_upperbound_b = avg_diam_ + const_term_avgspl + sqrt( 2*const_term_avgspl*avg_diam_ + pow(const_term_avgspl,2.) );
  bool debug_sp_lengths = false;

  // upper bound using empirical bernstein
  double var_estimate_avg_diam = 0.;
  for( int i=0; i <= graph_diameter; i++ ){
    for( int j=i+1; j <= graph_diameter; j++ ){
      var_estimate_avg_diam += pow(i-j,2.0)*sp_lengths[i]/(double)num_samples*sp_lengths[j]/(double)(num_samples-1);
    }
    if(verbose && debug_sp_lengths && sp_lengths[i] > 0){
      cout << "   sp_lengths["<<i<<"]: " << sp_lengths[i] << std::endl;
    }
  }
  log_term_avgspl = log(2./delta);
  double avg_diam_upperbound_eb = avg_diam_ + 7./3.*(graph_diameter-2)*log_term_avgspl/(double)num_samples + sqrt( 2*var_estimate_avg_diam*log_term_avgspl/(double)num_samples );

  double avg_diam_upperbound = min(avg_diam_upperbound_b , avg_diam_upperbound_eb);
  avg_diam_upperbound = min(avg_diam_upperbound , (double)(graph_diameter-2));

  if(verbose){
    cout << "avg_diam_: " << avg_diam_ << std::endl;
    cout << "avg_diam_upperbound_b: " << avg_diam_upperbound_b << std::endl;
    cout << "avg_diam_upperbound_eb: " << avg_diam_upperbound_eb << std::endl;
    cout << "var_estimate_avg_diam: " << var_estimate_avg_diam << std::endl;
  }
  return avg_diam_upperbound;

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
    start_time = get_time_sec();
    graph_diameter = estimate_diameter();
    //omp_set_num_threads(64);
    std::cout << "estimated diameter of the graph: " << graph_diameter << std::endl;
    std::cout << "time for estimating diameter " << get_time_sec() - start_time << std::endl;
    uint32_t tau = max(1. / err * (log(1. / delta)) , 1000.);
    tau = max((double)tau , 2 * graph_diameter * (log(1. / delta)) );
    std::cout << "Starting first pass. Tentative number of samples: " << tau << std::endl;
    num_samples = 0;
    sup_bcest = 0;
    sup_emp_wimpy_var = 0;
    void_samples = 0;
    second_phase_started = false;
    enable_emp_peel = empirical_peeling_a > 1;
    if(!enable_emp_peel){
      empirical_peeling_a = 100*tau;
    }


    if (union_sample == 0) {
        union_sample = min(get_nn(), (uint32_t) max( 2 * sqrt(get_ne()) / omp_get_max_threads(), k+20. ));
    }
    this->union_sample=union_sample;
    this->k=min(k, get_nn());

    sp_lengths = (uint64_t *) malloc( (graph_diameter+1)*sizeof(uint64_t) );
    for( int i=0; i <= graph_diameter; i++ ){
        sp_lengths[i] = 0;
    }

    last_output = get_time_sec();
    start_time = get_time_sec();
    this->top_k = new Ranking_list(union_sample);
    srand( SEED );
    uint32_t *random_seed = (uint32_t *) malloc( omp_get_max_threads()*sizeof(uint32_t) );
    for( int i=0; i < omp_get_max_threads(); i++ ){
        random_seed[i] = rand();
    }
    distinct_nodes_top_k = 0;

    firstpass = true;
    int samples_per_step = 10;

    // first phase sampling

    bool stop_first_pass = false;
    #pragma omp parallel
    {
        Sp_sampler sp_sampler( this, random_seed[omp_get_thread_num()] );
        while( !stop_first_pass ) {
            for (int i = 0; i <= samples_per_step; i++) {
                one_round(sp_sampler);
            }
            if( !stop_first_pass ){
              #pragma omp critical(stop)
              {
                if(absolute){
                  // check if sampled enough pairs
                  stop_first_pass = n_pairs >= tau;
                }
                else{
                  // check if at least k distinct nodes have been sampled
                  stop_first_pass = distinct_nodes_top_k >= 1.5*k && n_pairs >= tau;
                }
              }
            }
        }
    }

    std::cout << "First pass finished after " << num_samples << " iteration" << std::endl;

    // empirical peeling

    //int max_num_partitions = (int)(log((double)num_samples)/log(empirical_peeling_a)+1);
    number_of_non_empty_partitions = 0;
    map<int, int> non_empty_partitions;
    for (uint32_t i = 0; i < get_nn(); i++) {
      double emp_w_node = emp_wimpy_vars[i]/(double)num_samples;
      double min_inv_w_node = min(1./emp_w_node , (double)num_samples);
      int node_partition_idx = (int)(log(min_inv_w_node)/log(empirical_peeling_a)+1);
      if(!enable_emp_peel){
        node_partition_idx = 0;
      }
      partition_index[i] = node_partition_idx;
      non_empty_partitions[node_partition_idx] = non_empty_partitions[node_partition_idx] + 1;
    }
    number_of_non_empty_partitions = non_empty_partitions.size();
    //map<int, int> partitions_ids_map;
    int part_idx = 0;
    for (const auto &elem : non_empty_partitions){
      partitions_ids_map[elem.first] = part_idx;
      part_idx++;
    }

    sup_bcest_partition = (uint64_t*) calloc( number_of_non_empty_partitions , sizeof(uint64_t) );
    sup_empwvar_partition = (double*) calloc( number_of_non_empty_partitions , sizeof(double) );
    epsilon_partition = (double*) calloc( number_of_non_empty_partitions , sizeof(double) );
    max_mcera_partition = (int64_t*) calloc( mctrials*number_of_non_empty_partitions , sizeof(int64_t) );
    for(int i = 0; i < number_of_non_empty_partitions; i++){
      sup_bcest_partition[i] = 0;
      sup_empwvar_partition[i] = 0;
      epsilon_partition[i] = 1.0;
      for(uint32_t j = 0; j < mctrials; j++){
        max_mcera_partition[i*mctrials+j] = 0;
      }
    }
    std::cout << "number_of_non_empty_partitions  " << number_of_non_empty_partitions << std::endl;
    for (const auto &elem : non_empty_partitions){
      std::cout << "  part. w. index  " << elem.first << " has " << elem.second << " elements, map to " << partitions_ids_map[elem.first] << std::endl;
    }


    *time_bfs = 0;
    *time_critical = 0;
    *time_critical_round = 0;
    *time_mcera = 0;

    std::cout << "time for first pass " << get_time_sec() - start_time << std::endl;
    std::cout << "void_samples first pass " << void_samples << std::endl;
    void_samples = 0;

    firstpass = false;
    double max_num_samples = 0.;

    double time_required_second_pass = get_time_sec();

    // for the absolute approximation, compute upper bound to number of samples
    if(absolute){

      // compute upper bound on average shortest path length
      double avg_diam_upperbound = getUpperBoundAvgDiameter(delta/8. , true);

      // reset counters of shortest path lengths
      for( int i=0; i <= graph_diameter; i++ ){
          sp_lengths[i] = 0;
      }

      // compute upper bound on top-1 bc
      double top1bc_upperbound = getUpperBoundTop1BC(sup_bcest/(double)num_samples , delta/8.);
      double wimpy_var_upperbound = getUpperBoundTop1BC(sup_emp_wimpy_var/(double)num_samples , delta/8.);

      if(alpha_sp_given == false && wimpy_var_upperbound > 0.9 * top1bc_upperbound){
        alpha_sp_sampling = -1.;
      }

      // compute upper limit on number of samples
      if(enable_m_hat){
        max_num_samples = getUpperBoundSamples(top1bc_upperbound , wimpy_var_upperbound , avg_diam_upperbound , err , delta/2.);
      }
      else{
        max_num_samples = 0;
      }

      // number of samples with VC dim bound
      double VC_bound = 0.5 / err / err * (log2(graph_diameter-1) + 1 + log(2. / delta));

      std::cout << "sup_bcest: " << sup_bcest/(double)num_samples << endl;
      std::cout << "sup_emp_wimpy_var: " << sup_emp_wimpy_var/(double)num_samples << std::endl;
      std::cout << "top1bc_upperbound: " << top1bc_upperbound << std::endl;
      std::cout << "wimpy_var_upperbound: " << wimpy_var_upperbound << std::endl;
      std::cout << "max_num_samples: " << max_num_samples << std::endl;
      std::cout << "VC_bound: " << VC_bound << std::endl;
      std::cout << "alpha_sp_sampling: " << alpha_sp_sampling << std::endl;
    }

    // reset variables for second pass
    iteration_index = 1;
    n_pairs = 0;
    double est_kth_bc = 1.0;
    if(!absolute){
      est_kth_bc = this->top_k->get_value(k)/(double)num_samples + 1.0/(double)num_samples;
      cout << "est_kth_bc: " << est_kth_bc << std::endl;
    }
    delete(this->top_k);
    this->top_k = new Ranking_list(union_sample);
    uint32_t v_mc_index;
    for (uint32_t i = 0; i < get_nn(); i++) {
        approx[i] = 0;
        approx_toadd[i] = 0;
        emp_wimpy_vars[i] = 0;
        // mcrade
        v_mc_index = i*mctrials;
        for(uint32_t j = 0; j < mctrials; j++){
          mcrade[v_mc_index+j] = 0;
        }
    }
    sup_bcest = 0;
    omega = pow(10.,15);
    if(max_num_samples > 0.){
      omega = max_num_samples;
    }


    // guess a first sample size according to what we computed in the first phase
    first_stopping_samples = 0;//2./err/err*( highest_freq*log(2./delta) );
    double eps_guess = 1.;
    // start a binary search to find a good starting sample size
    double first_sample_lower = 1./err*log(2./delta);
    double first_sample_upper = omega;
    bool debug_first_sample_size_heuristic = false;
    double sup_emp_wimpy_var_norm = sup_emp_wimpy_var/(double)num_samples + 1/(double)num_samples;
    while(first_sample_upper - first_sample_lower > 10){
      if(debug_first_sample_size_heuristic) cout << "--- LOG SEARCH ITER" << std::endl;
      num_samples = (first_sample_upper + first_sample_lower)/2.;
      if(absolute){
        eps_guess = sqrt(2*sup_emp_wimpy_var_norm*log(2./delta)/(double)num_samples) + log(2./delta)/(double)num_samples/3.;
      }
      else{
        eps_guess = sqrt(2*est_kth_bc*log(2./delta/est_kth_bc)/(double)num_samples) + log(2./delta/est_kth_bc)/(double)num_samples/3.;
      }
      if(debug_first_sample_size_heuristic) {
        double term1 = sqrt(2*sup_emp_wimpy_var_norm*log(2./delta)/(double)num_samples);
        double term2 = log(2./delta)/(double)num_samples/3.;
        cout << "  eps_guess: " << eps_guess << std::endl;
        cout << "  num_samples: " << num_samples << std::endl;
        cout << "  sup_emp_wimpy_var_norm: " << sup_emp_wimpy_var_norm << std::endl;
        cout << "  term1: " << term1 << std::endl;
        cout << "  term2: " << term2 << std::endl;
      }
      if(absolute){
        if(eps_guess > err){
          first_sample_lower = num_samples;
          if(debug_first_sample_size_heuristic) cout << "lower increased at: " << first_sample_lower << std::endl;
        }
        else{
          first_sample_upper = num_samples;
          if(debug_first_sample_size_heuristic) cout << "upper decreased at: " << first_sample_upper << std::endl;
        }
      }
      else{
        bool violated_condition_lower_1 = est_kth_bc+eps_guess > est_kth_bc/(1.-err);
        bool violated_condition_upper_1 = est_kth_bc-eps_guess < est_kth_bc/(1.+err);
        bool violated_condition_lower_2 = eps_guess >= est_kth_bc*err*pow((1-err)/(1+err),2.);
        if(violated_condition_lower_1 || violated_condition_upper_1 || violated_condition_lower_2){
          first_sample_lower = num_samples;
          if(debug_first_sample_size_heuristic) cout << "(top-k) lower increased at: " << first_sample_lower << std::endl;
        }
        else{
          first_sample_upper = num_samples;
          if(debug_first_sample_size_heuristic) cout << "(top-k) upper decreased at: " << first_sample_upper << std::endl;
        }
      }
    }



    first_stopping_samples = num_samples;
    num_samples = 0;
    sup_emp_wimpy_var = 0;

    if(enable_m_hat){
      double max_num_iterations = log((double)max_num_samples/(double)first_stopping_samples)/log(1.2)+1;
      delta_for_progressive_bound = delta/(floor(max_num_iterations));
    }




    cout << "first_stopping_samples: " << first_stopping_samples << std::endl;
    // we stop at omega
    last_stopping_samples = omega;
    cout << "last_stopping_samples: " << last_stopping_samples << std::endl;
    if(first_stopping_samples >= last_stopping_samples/4){
      first_stopping_samples = last_stopping_samples/4;
      cout << "first_stopping_samples dropped to " << first_stopping_samples << std::endl;
    }
    // initialize the next at the first, it will be updated during iterations
    next_stopping_samples = first_stopping_samples;
    bool stop_mcrade = false;
    last_output = get_time_sec();
    second_phase_started = true;

    // start second phase
    cout << "Starting second phase... " << std::endl;

    #pragma omp parallel
    {
        Sp_sampler sp_sampler( this, random_seed[omp_get_thread_num()] );
        Status status(union_sample);
        status.n_pairs = 0;

        while( !stop_mcrade ) {
            for (uint32_t i = 0; i <= 10; i++) {
                one_round(sp_sampler);
            }

            // stop if enough samples have been processed
            if(enable_m_hat && num_samples >= omega){
              #pragma omp critical(stopcond)
              {
                stop_mcrade = true;
              }
            }
            // check stopping condition for mcrade
            if (!stop_mcrade && num_samples < last_stopping_samples && num_samples >= next_stopping_samples) {
                get_status (&status);
                #pragma omp critical(stopcond)
                {
                    if (!stop_mcrade  && num_samples < last_stopping_samples && num_samples >= next_stopping_samples) {
                        stop_mcrade = compute_finished_mcrade(&status);
                        if(stop_mcrade){
                          std::cout << "/* stop_mcrade is true */" << std::endl;
                        }
                        else{
                          next_stopping_samples = get_next_stopping_sample();
                        }
                    }
                }
            }

            double current_time = get_time_sec();

            if (!stop_mcrade && verbose > 0 && current_time - last_output > verbose) {
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
    std::cout << "time for second pass " << get_time_sec() - time_required_second_pass << std::endl;
    std::cout << "void_samples second pass " << void_samples << " (" << (double)void_samples/(double)num_samples << ")" << std::endl;

    if (verbose > 0) {
        if(!absolute){
          cout << "eps_final_topk " << eps_final_topk << std::endl;
        }
        Status status(union_sample);
        get_status(&status);
        print_status(&status, true);
    }

    bool output_results = output_file.length() >= 1;

    if(output_results){
      std::ofstream output_file_str;
      output_file_str.open(output_file);
      if(!absolute){
        for(uint32_t i=0; i < this->numresults_topk; i++){
          uint64_t node_id = this->top_k->get(i);
          double approx_node = this->top_k->get_value(i)/(double)num_samples;
          output_file_str << node_id << ","<< approx_node << "\n";
        }
      }
      else{
          for(uint32_t i=0; i < get_nn(); i++){
            uint64_t node_id = i;
            double approx_node = approx[i]/(double)num_samples;
            output_file_str << node_id << ","<< approx_node << "\n";
          }
      }
      output_file_str.close();
    }

    n_pairs += tau;
}

// Destructor of the class Probabilistic.
Probabilistic::~Probabilistic() {
    free(approx);
    delete(top_k);
    // mcrade
    free(mcrade);
    free(emp_wimpy_vars);
    delete(mcrade_randgen);
    free(partition_index);
    free(epsilon_partition);
    free(sup_bcest_partition);
    free(sup_empwvar_partition);
    free(max_mcera_partition);
    free(sp_lengths);
}
