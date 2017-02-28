#ifndef _IF_SETUP_H
#define _IF_SETUP_H

extern brain pfc;
//const int num_pools = 8;
extern struct poolinfo pool_inh, pool_exc;

void setup(int num_neurons, int num_pools, int num_ext_neurns,
	   double time_step, double inh_percent, double coding_lvl,
	   double g_AMPA_ratio, double ratio);
void set_pool_data();
void set_all_connections(int num_pools, double coding_lvl, double w_plus, double w_I, double w_T, string connectivityFileName);

typedef struct State State;
struct State {
    double nu_inh;
    double nu_1;
    double nu_2;
    double nu_nsl;
};
void setInitState(const State* s);
void getState(State* tmp);
void copyState(State* tmp, const State* s);
void sum_s(State* tmp, const State* s1, const State* s2);
void dif_s(State* tmp, const State* s1, const State* s2);
void set_midpoint(State* tmp, State* s1, State* s2);
double fetch_rate(const State* s, int pool);
void put_rate(State* s, int pool, double rate);
double find_boundary(const State* s1, const State* s2, int frozen_pool);
void find_null_point(const State* s, int frozen_pool);
#endif
