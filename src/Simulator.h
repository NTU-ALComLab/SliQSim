#ifndef _SIMULATOR_H_
#define _SIMULATOR_H_

#include <iostream>
#include <stdio.h> // FILE
#include <unordered_map>
#include <sys/time.h> //estimate time
#include <fstream> //fstream
#include <sstream> // int to string
#include <cstdlib> //atoi
#include <string> //string
#include <sstream>
#include <random>
#include <cmath>
#include <vector>
#include "../cudd/cudd/cudd.h"
#include "../cudd/cudd/cuddInt.h"
#include "../cudd/util/util.h"


#define PI 3.14159265358979323846264338327950288419716939937510582097494459230781640628620899

class Simulator
{
public:
    // constructor and destructor
    Simulator(int type, int nshots, int seed) : 
    n(0), r(sizeof(int) * 8), w(4), k(0), inc(3), error(0), 
    normalize_factor(1), gatecount(0), NodeCount(0), isMeasure(0), measured_qubits(0), shots(nshots)
    , sim_type(type), gen(std::default_random_engine(seed)){
    }
    ~Simulator()  { 
        clear();
    }

    /* gates */
    void Toffoli(int targ, std::vector<int> &cont, std::vector<int> &ncont);
    void Fredkin(int swapA , int swapB, std::vector<int> &cont);
    void Peres(int a, int b, int c);
    void Peres_i(int a, int b, int c);
    void Hadamard(int iqubit);
    // void X_1_2(int iqubit); // have not been modified for modified algebraic
    // void Y_1_2(int iqubit); // have not been modified for modified algebraic
    void rx_pi_2(int iqubit); 
    void ry_pi_2(int iqubit);
    void Phase_shift(int phase, int *iqubit, int napplied); // phase can only be the power of 2
    void Phase_shift_dagger(int phase, int *iqubit, int napplied);
    void PauliX(int iqubit);
    void PauliY(int iqubit);
    void PauliZ(int *iqubit, int napplied); // Z or CZ

    /* measurement */
    double measure_probability(DdNode *node, int kd2, int nVar, int nAnci_fourInt, int edge);
    void measure_one(int position, int kd2, double H_factor, int nVar, int nAnci_fourInt, std::string *outcome);
    void measurement();
    void getStatevector();

    /* simulation */
    void sim_qasm_file(std::string qasm);
    void sim_qasm(std::string qasm);

    /* misc */
    void init_state(int *constants);
    void decode_entries();
    void reorder();
    void alloc_BDD(bool extend);
    int overflow(DdNode *g, DdNode *h, DdNode *crin);
    void nodecount();
    void print_info(double runtime, size_t memPeak);

private:
    DdManager *manager;
    DdNode ***All_Bdd;
    int n; // # of qubits
    int r; // resolution of integers
    int w; // # of integers
    int k; // k in algebraic representation
    int inc; // add inc BDDs when overflow occurs, used in alloc_BDD
    int shots;
    int sim_type; // 0: statevector, 1: measure
    bool isMeasure;
    std::vector<int> measured_qubits; // i-th element is the i-th qubit measured
    int *measured_qubits_to_clbits; // -1 if not measured
    std::string measure_outcome;
    double normalize_factor; // normalization factor used in measurement
    DdNode *bigBDD; // big BDD used if measurement
    std::default_random_engine gen; // random generator
    std::unordered_map<DdNode *, double> Node_Table; // key: node, value: summed prob
    std::unordered_map<std::string, int> state_count;
    std::string statevector;
    std::string run_output; // output string for Qiskit

    unsigned long gatecount;
    unsigned long NodeCount;
    double error;
    struct timeval t1, t2;

    // Clean up Simulator
    void clear() {
        for (int i = 0; i < w; i++)
            for (int j = 0; j < r; j++)
                Cudd_RecursiveDeref(manager, All_Bdd[i][j]);
        for (int i = 0; i < w; i++)
            delete[] All_Bdd[i];
        delete [] All_Bdd;
        delete [] measured_qubits_to_clbits;
        measured_qubits.clear();
        measure_outcome.clear();
        Node_Table.clear();
        state_count.clear();
        Cudd_Quit(manager);
    };
};

#endif
