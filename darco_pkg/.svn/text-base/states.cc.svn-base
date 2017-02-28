#include "brain.h"
#include "setup.h"
#include "math.h"

void setInitState(const State* s){
    pfc.setStartPoolrate(0, s->nu_inh);
    pfc.setStartPoolrate(1, s->nu_1);
    pfc.setStartPoolrate(2, s->nu_2);
    pfc.setStartPoolrate(3, s->nu_nsl);
}

// Lleig
void getState(State* tmp){
    tmp->nu_inh = pfc.getPoolrate(0);
    tmp->nu_1 = pfc.getPoolrate(1);
    tmp->nu_2 = pfc.getPoolrate(2);
    tmp->nu_nsl = pfc.getPoolrate(3);
}

void copyState(State* tmp, const State* s){
    // I think this is not needed: tmp=s works fine
    tmp->nu_inh = s->nu_inh;
    tmp->nu_1 = s->nu_1;
    tmp->nu_2 = s->nu_2;
    tmp->nu_nsl = s->nu_nsl;
}

void sum_s(State* tmp, const State* s1, const State* s2){
    tmp->nu_inh = (s1->nu_inh + s2->nu_inh);
    tmp->nu_1 = (s1->nu_1 + s2->nu_1);
    tmp->nu_2 = (s1->nu_2 + s2->nu_2);
    tmp->nu_nsl = (s1->nu_nsl + s2->nu_nsl);
}

void dif_s(State* tmp, const State* s1, const State* s2){
    tmp->nu_inh = (s1->nu_inh - s2->nu_inh);
    tmp->nu_1 = (s1->nu_1 - s2->nu_1);
    tmp->nu_2 = (s1->nu_2 - s2->nu_2);
    tmp->nu_nsl = (s1->nu_nsl - s2->nu_nsl);
}

void set_midpoint(State* tmp, State* s1, State* s2){
    // Given two states s1 and s2, set tmp to the midpoint (in
    // the 4D space (nu0,nu1,nu2,nu3)) between the two.
    tmp->nu_inh = (s1->nu_inh + s2->nu_inh)/2.0;
    tmp->nu_1 = (s1->nu_1 + s2->nu_1)/2.0;
    tmp->nu_2 = (s1->nu_2 + s2->nu_2)/2.0;
    tmp->nu_nsl = (s1->nu_nsl + s2->nu_nsl)/2.0;
}

double fetch_rate(const State* s, int pool){
    if(pool==1){
	return s->nu_1;
    } else if(pool==2){
	return s->nu_2;
    } else
	exit(2);
}

void put_rate(State* s, int pool, double rate){
    if(pool==1){
	s->nu_1 = rate;;
    } else if(pool==2){
	s->nu_2 = rate;
    } else
	exit(2);
}

void find_null_point(const State* s, int frozen_pool)
{
    // Find the fixed states starting from different
    // initial conditions.
    setInitState(s);
    pfc.resetPools();

    // Integrate the diff equation quenching pool "frozen_pool"
    pfc.get_fixed_rates(0.1, frozen_pool);
}

double find_boundary(const State* s1, const State* s2, int frozen_pool)
{
    int not_frozen = (frozen_pool == 1) ? 2 : 1;
    double difference, tmp_rate;

    double v1, v2;
    v1 = fetch_rate(s1, not_frozen);
    v2 = fetch_rate(s2, not_frozen);

    State midpoint, up, dn;
    copyState(&up, s1);
    copyState(&dn, s2);
  
    // Set initial conditions in the (4D) midpoint of the 2 states
    difference = fabs(v1 - v2);

    while (difference > 1e-5) {
	set_midpoint(&midpoint, &up, &dn);
	setInitState(&midpoint);
	pfc.resetPools();
	pfc.get_fixed_rates(0.1, frozen_pool);
	tmp_rate = pfc.getPoolrate(not_frozen);

	if (fabs(tmp_rate - v1) < 1.0) {
	    up = midpoint;
	    //copyState(&up, &midpoint);
	}
	else if (fabs(tmp_rate - v2) < 1.0) {
	    // copyState(&dn, &midpoint);
	    dn = midpoint;
	}
	else {
	    cerr << "Not converged to any of the stable points\n";
	    cerr << tmp_rate << "\t" << v1 << "\t" << v2 << "\n";
	    //exit(2);
	    break;
	}
	difference = fabs(fetch_rate(&up,not_frozen) 
		- fetch_rate(&dn,not_frozen))/2.0;
    }
    set_midpoint(&midpoint, &up, &dn);
    return fetch_rate(&midpoint, not_frozen);
}

