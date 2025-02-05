#ifndef _SIMULATOR_H_
#define _SIMULATOR_H_

#include <iostream>
#include <stdio.h>    // FILE
#include <unordered_map>
#include <unordered_set>
#include <tuple>
#include <sys/time.h> // estimate time
#include <fstream>    // fstream
#include <sstream>    // int to string
#include <cstdlib>    // atoi
#include <string>     // string
#include <regex>      // replace
#include <iomanip>    // setw
#include <sstream>
#include <random>
#include <cmath>
#include <vector>
#include <gmpxx.h>
#include "../cudd/cudd/cudd.h"
#include "../cudd/cudd/cuddInt.h"
#include "../cudd/util/util.h"

#define KW_ASSIGN     "assign"
#define KW_DIST       "dist"
#define KW_AMP        "amp"
#define KW_BOOL       "bf"
#define KW_EXPT       "expt"
#define KW_INTEQ      "inteq"
#define KW_INTNEQ     "intneq"
#define KW_INTGT      "intgt"
#define KW_INTLT      "intlt"
#define KW_HWEQ       "hweq"
#define KW_HWNEQ      "hwneq"
#define KW_HWGT       "hwgt"
#define KW_HWLT       "hwlt"
#define KW_WS         "weightedsum"
#define KW_EWS        "endweightedsum"
#define KW_BTN        "between"
#define KW_OOF        "outof"
#define KW_GEQ        "geq"
#define KW_LEQ        "leq"


#define PI 3.14159265358979323846264338327950288419716939937510582097494459230781640628620899

class Simulator
{
public:
    // constructor and destructor
    Simulator(int type, int nshots, int seed, int bitSize, bool reorder, bool isQuery, bool alloc) :
        n(0), r(bitSize), w(4), k(0), inc(3), shift(0), error(0),
        normalize_factor(1), gatecount(0), NodeCount(0), isMeasure(0), isQuery(isQuery), shots(nshots), isReorder(reorder), isAlloc(alloc)
        , sim_type(type), statevector("null"), gen(std::default_random_engine(seed)) {}
    Simulator(int nshots, int seed, int bitSize, bool reorder, bool isQuery, bool alloc) :
        n(0), r(bitSize), w(4), k(0), inc(3), shift(0), error(0),
        normalize_factor(1), gatecount(0), NodeCount(0), isMeasure(0), isQuery(isQuery), shots(nshots), isReorder(reorder), isAlloc(alloc)
        , sim_type(0), statevector("null"), gen(std::default_random_engine(seed)) {}
    ~Simulator()  {
        clear();
    }

    /* gates */
    void Toffoli(int targ, std::vector<int> cont, std::vector<int> ncont);
    void Fredkin(int swapA , int swapB, std::vector<int> cont);
    void Peres(int a, int b, int c);
    void Peres_i(int a, int b, int c);
    void Hadamard(int iqubit);
    void rx_pi_2(int iqubit);
    void ry_pi_2(int iqubit);
    void Phase_shift(int phase, int iqubit); // phase can only be 2 to the power of an integer
    void Phase_shift_dagger(int phase, int iqubit);
    void PauliX(int iqubit);
    void PauliY(int iqubit);
    void PauliZ(std::vector<int> iqubit); // Z or CZ
    void measure(int qreg, int creg);

    /* measurement */
    void measurement_obs(std::string obsfile);
    void measurement();
    void getStatevector();

    /* simulation */
    void init_simulator(int n);
    void sim_qasm_file(std::string qasm);
    void sim_qasm(std::string qasm);
    void print_results();

    /* misc */
    void reorder();
    void decode_entries();
    void print_info(double runtime, size_t memPeak);

private:
    DdManager *manager;
    DdNode ***All_Bdd;
    int n; // # of qubits
    int r; // resolution of integers
    int w; // # of integers
    int k; // k in algebraic representation
    int inc; // add inc BDDs when overflow occurs, used in alloc_BDD
    int shift; // # of right shifts
    int shots;
    int sim_type; // 0: statevector, 1: measure
    bool isMeasure;
    bool isQuery;
    bool isReorder;
    bool isAlloc;
    int nClbits;
    std::vector<std::vector<int>> measured_qubits_to_clbits; // empty if not measured
    std::string measure_outcome;
    double normalize_factor; // normalization factor used in measurement
    DdNode *bigBDD; // big BDD used if measurement
    std::default_random_engine gen; // random generator
    std::unordered_map<DdNode *, double> Node_Table; // key: node, value: summed prob
    std::unordered_map<std::string, int> state_count;
    std::string statevector;
    std::vector<std::tuple<std::string, std::string>> property;
    std::string run_output; // output string for Qiskit
    std::unordered_map<std::string, DdNode*> defined_var;
    const std::unordered_set<std::string> keyword_truth = {KW_BOOL, KW_INTEQ, KW_INTNEQ, KW_INTGT, KW_INTLT, KW_HWEQ, KW_HWNEQ, KW_HWGT, KW_HWLT};
    const std::unordered_set<std::string> keyword_value = {KW_BOOL, KW_INTEQ, KW_INTNEQ, KW_INTGT, KW_INTLT, KW_HWEQ, KW_HWNEQ, KW_HWGT, KW_HWLT, KW_EXPT};
    const std::unordered_set<std::string> keyword_condition = {KW_BTN, KW_OOF, KW_GEQ, KW_LEQ};
    std::string condition_stack = "";
    std::vector<double> range = {0, 0};

    unsigned long gatecount;
    unsigned long NodeCount;
    double error;

    /* measurement */
    void create_bigBDD();
    //void record_property_value(std::string& op_name, std::vector<std::string>& followings);
    int handle_property(std::string& raw_line, std::vector<std::string>& words, bool is_in_ws);
    DdNode* get_property(std::string& op_name, std::vector<std::string>& followings);
    std::string getCondStr(std::string& line);
    std::string getCondRslt(double value);
    
    double get_prob(DdNode* function);
    std::string get_amplitude_string(int *assign);
    double measure_probability(DdNode *node, int kd2, int nVar, int nAnci_fourInt, int edge);
    void measure_one(int position, int kd2, double H_factor, int nVar, int nAnci_fourInt, std::string *outcome);

    /* misc */
    void init_state(int *constants);
    void alloc_BDD(DdNode ***Bdd, bool extend);
    void dropLSB(DdNode ***Bdd);
    int overflow3(DdNode *g, DdNode *h, DdNode *crin);
    int overflow2(DdNode *g, DdNode *crin);
    void nodecount();
    std::vector<std::string> boolean_parser(std::string& inStr);
    DdNode* node_equiv(std::vector<DdNode*>& int_1, std::vector<DdNode*>& int_2);
    DdNode* node_larger(std::vector<DdNode*>& int_1, std::vector<DdNode*>& int_2);
    DdNode* func2node(std::vector<std::string>& func);

    // Clean up Simulator
    void clear() {
        if (!manager) return;    // hasn't initialized yet; maybe due to lack of qasm file
        
        for (int i = 0; i < w; i++)
            for (int j = 0; j < r; j++)
                Cudd_RecursiveDeref(manager, All_Bdd[i][j]);
        for (int i = 0; i < w; i++)
            delete[] All_Bdd[i];
        delete [] All_Bdd;
        if (isMeasure == 1)
            Cudd_RecursiveDeref(manager, bigBDD);
        for (auto& it: defined_var) {
            Cudd_RecursiveDeref(manager, it.second);
        }
        defined_var.clear();
        measured_qubits_to_clbits.clear();
        measure_outcome.clear();
        Node_Table.clear();
        state_count.clear();
        Cudd_Quit(manager);
    };
};

#endif
