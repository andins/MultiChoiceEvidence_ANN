#ifndef _IF_STATES_H
#define _IF_STATES_H

extern brain pfc;

typedef struct State State;
struct State {
    double nu_inh;
    double nu_1;
    double nu_2;
    double nu_nsl;
};

#endif
