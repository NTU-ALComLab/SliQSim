#include "Simulator.h"
#include "util_sim.h"

/**Function*************************************************************

  Synopsis    [calculate probabilities for measurement]

  Description []
               
  SideEffects []

  SeeAlso     []

***********************************************************************/
double Simulator::measure_probability(DdNode *node, int kd2, int nVar, int nAnci_fourInt, int edge)
{
    double then_edge, else_edge, probability;
    double re, im;
    DdNode *child = Cudd_Child(manager, node, edge);
    Cudd_Ref(child);
    int position_node = Cudd_ReadPerm(manager, Cudd_NodeReadIndex(node)), position_child = Cudd_ReadPerm(manager, Cudd_NodeReadIndex(child));
    int skip_level, i, j, k;
    std::unordered_map<DdNode *, double>::iterator it;
    it = Node_Table.find(child);

    /* deal with the nodes which are reduced */
    if (position_child > n)
        skip_level = n - position_node - 1;
    else
        skip_level = position_child - position_node - 1;

    /* compute probability of the node */
    if (it != Node_Table.end()) // node has been recorded
    {
        Cudd_RecursiveDeref(manager, child);
        return it->second * pow(2, skip_level);
    }
    else if (Cudd_IsConstant(child)) // constant
    {
        if (!(Cudd_IsComplement(child))) // constant 1
        {
            re = 0;
            im = 0;
            for (i = 0; i < w; i++)
            {
                re -= cos(i/w * PI);
                im -= sin(i/w * PI);
            }
            probability = pow(re, 2) + pow(im, 2);
            Cudd_RecursiveDeref(manager, child);
            return probability * pow(2, n - position_node - 1);
        }
        else // constant 0
        {
            Cudd_RecursiveDeref(manager, child);
            return 0;
        }
    }
    else
    {
        if (position_child >= n) // compute entry
        {
            double int_value;
            DdNode *tmp;
            int oneEntry;
            re = 0;
            im = 0;
            int *assign = new int[nVar];
            for (i = 0; i < nVar; i++)
                assign[i] = 0;
            /* TODO: BDD to truth table */
            for (i = 0; i < w; i++) //compute each complex value
            {
                int_value = 0;
                for (j = 0; j < r; j++) //compute each integer
                {
                    tmp = Cudd_Eval(manager, child, assign);
                    Cudd_Ref(tmp);
                    oneEntry = !(Cudd_IsComplement(tmp));
                    Cudd_RecursiveDeref(manager, tmp);
                    if (j == r - 1)
                        int_value -= oneEntry * pow(2, j - kd2);
                    else
                        int_value += oneEntry * pow(2, j - kd2);
                    full_adder_plus_1_start(nVar, assign, n + nAnci_fourInt);
                }
                /* translate to re and im */
                re += int_value * cos((double) (w - i - 1)/w * PI);
                im += int_value * sin((double) (w - i - 1)/w * PI);
                full_adder_plus_1_start(nVar, assign, n);
                for (j = n + nAnci_fourInt; j < nVar; j++) // reset array: not necessary but straightforward
                    assign[j] = 0;
            }
            probability = pow(re, 2) + pow(im, 2);
            delete[] assign;
        }
        else // trace edges
        {
            then_edge = measure_probability(child, kd2, nVar, nAnci_fourInt, 1);
            else_edge = measure_probability(child, kd2, nVar, nAnci_fourInt, 0);
            probability = then_edge + else_edge;
        }
    }

    Node_Table[child] = probability;
    Cudd_RecursiveDeref(manager, child);

    return probability * pow(2, skip_level);
}

/**Function*************************************************************

  Synopsis    [measure one qubit]

  Description []
               
  SideEffects []

  SeeAlso     []

***********************************************************************/
void Simulator::measure_one(int position, int kd2, double H_factor, int nVar, int nAnci_fourInt, std::string *outcome)
{
    int index = Cudd_ReadInvPerm(manager, position);
    int comple;
    int noNode_f = 0; //flag: 1 if the node we want to measure is reduced
    std::uniform_real_distribution<double> dis(0.0, 1.0);
	double p0, p1, p;
    double rand_num;
    double epsilon = 0.001;
    DdNode *node, *tmp;
    Cudd_Ref(bigBDD); // if this line is not added, there might be some error in CUDD.

    /* traverse until the qubit measured */
    comple = Cudd_IsComplement(bigBDD);
    tmp = Cudd_Regular(bigBDD);
    while (!(tmp->index == index))
    {
        if (cuddIsConstant(tmp))
        {
            noNode_f = 1;
            break;
        }
        if ((tmp->index < n) && (outcome->substr(n - 1 - tmp->index, 1) == "1"))
        {
            tmp = cuddT(tmp);
        }
        else
        {
            comple ^= Cudd_IsComplement(cuddE(tmp));
            tmp = Cudd_Regular(cuddE(tmp));
        }
    }

    if (noNode_f)
    {
        p0 = 0.5;
        p1 = 0.5;
    }
    else
    {
        node = (Cudd_NotCond(tmp, comple));
        Cudd_Ref(node);
        p0 = measure_probability(node, kd2, nVar, nAnci_fourInt, 0);
        p0 *= H_factor * H_factor * normalize_factor * normalize_factor;
        p1 = measure_probability(node, kd2, nVar, nAnci_fourInt, 1);
        p1 *= H_factor * H_factor * normalize_factor * normalize_factor;
        Cudd_RecursiveDeref(manager, node);
        p = p0 + p1;
        if (abs(p - 1) > epsilon)
		{
            std::cerr << "[error]Numerical error: p0 + p1 = " << p << ", not 1" << std::endl;
			std::exit(1);
		}
    }
	
    double error_tmp = abs(p0 + p1 - 1);
    if (error_tmp > error)
	    error = error_tmp;

    /* sample */
    rand_num = dis(gen);
    if (rand_num > p0)
    {
        (*outcome)[n - 1 - index] = '1'; // LSB: q0
        // (*outcome)[index] = '1'; // LSB: qn-1
        normalize_factor /= sqrt(p1);
    }
    else
    {
        // (*outcome)[n - 1 - index] = '0'; // LSB: q0
        // (*outcome)[index] = '0'; // LSB: qn-1
        normalize_factor /= sqrt(p0);
    }
}

/**Function*************************************************************

  Synopsis    [measurement]

  Description []
               
  SideEffects []

  SeeAlso     []

***********************************************************************/
void Simulator::measurement()
{
    double oneroot2 = 1 / sqrt(2);
    double H_factor = pow(oneroot2, k%2);
    int nAnci_oneInt = ceil(log(r) / log(2)), nAnci_fourInt = ceil(log(w) / log(2)), nAnci = nAnci_oneInt + nAnci_fourInt, nnAnci_fourInt = n + nAnci_fourInt, nVar = n + nAnci;
    DdNode *tmp1, *tmp2, *tmp3;
    int i, j, h;
    
    int *arrAnci_fourInt = new int[nAnci_fourInt];
    for (i = 0; i < nAnci_fourInt; i++)
        arrAnci_fourInt[i] = 0;
    int *arrAnci_oneInt = new int[nAnci_oneInt];
    for (i = 0; i < nAnci_oneInt; i++)
        arrAnci_oneInt[i] = 0;
    
    bigBDD = Cudd_Not(Cudd_ReadOne(manager));
    Cudd_Ref(bigBDD);
    for (i = 0; i < w; i++)
    {
        tmp3 = Cudd_Not(Cudd_ReadOne(manager));
        Cudd_Ref(tmp3);
        for (j = 0; j < r; j++)
        {
            tmp1 = Cudd_ReadOne(manager);
            Cudd_Ref(tmp1);
            for (h = n + nAnci - 1; h >= nnAnci_fourInt; h--)
            {
                if (arrAnci_oneInt[h - nnAnci_fourInt])
                    tmp2 = Cudd_bddAnd(manager, Cudd_bddIthVar(manager, h), tmp1);
                else
                    tmp2 = Cudd_bddAnd(manager, Cudd_Not(Cudd_bddIthVar(manager, h)), tmp1);
                Cudd_Ref(tmp2);
                Cudd_RecursiveDeref(manager, tmp1);
                tmp1 = tmp2;
            }
            tmp2 = Cudd_bddAnd(manager, All_Bdd[i][j], tmp1);
            Cudd_Ref(tmp2);
            Cudd_RecursiveDeref(manager, tmp1);
            // Cudd_RecursiveDeref(manager, All_Bdd[i][j]); // keep these BDDs for statevector after measurement
            tmp1 = tmp2;
            tmp2 = Cudd_bddOr(manager, tmp3, tmp1);
            Cudd_Ref(tmp2);
            Cudd_RecursiveDeref(manager, tmp1);
            Cudd_RecursiveDeref(manager, tmp3);
            tmp3 = tmp2;
            full_adder_plus_1(nAnci_oneInt, arrAnci_oneInt);
        }
        tmp1 = Cudd_ReadOne(manager);
        Cudd_Ref(tmp1);
        for (j = nnAnci_fourInt - 1; j >= n; j--)
        {
            if (arrAnci_fourInt[j - n])
                tmp2 = Cudd_bddAnd(manager, Cudd_bddIthVar(manager, j), tmp1);
            else
                tmp2 = Cudd_bddAnd(manager, Cudd_Not(Cudd_bddIthVar(manager, j)), tmp1);
            Cudd_Ref(tmp2);
            Cudd_RecursiveDeref(manager, tmp1);
            tmp1 = tmp2;
        }
        tmp2 = Cudd_bddAnd(manager, tmp3, tmp1);
        Cudd_Ref(tmp2);
        Cudd_RecursiveDeref(manager, tmp1);
        Cudd_RecursiveDeref(manager, tmp3);
        tmp1 = tmp2;
        tmp2 = Cudd_bddOr(manager, bigBDD, tmp1);
        Cudd_Ref(tmp2);
        Cudd_RecursiveDeref(manager, tmp1);
        Cudd_RecursiveDeref(manager, bigBDD);
        bigBDD = tmp2;
        full_adder_plus_1(nAnci_fourInt, arrAnci_fourInt);
        for (j = 0; j < nAnci_oneInt; j++) // reset array: not necessary but straightforward
            arrAnci_oneInt[j] = 0;
    }
    nodecount();

    // move measured qubits to the top
    int *permutation = new int[nVar];
    int indCount1 = 0;
    int indCount0 = measured_qubits.size();
    for (i = 0; i < n; i++)
    {
        int index = Cudd_ReadInvPerm(manager, i);
        if (measured_qubits_to_clbits[index] != -1)
        {
            permutation[indCount1] = index;
            indCount1++;
        }
        else
        {
            permutation[indCount0] = index;
            indCount0++;
        }
    }
    for (i = n; i < nVar; i++)
        permutation[i] = Cudd_ReadInvPerm(manager, i);
    int dum = Cudd_ShuffleHeap(manager, permutation);
    nodecount();

    std::unordered_map<std::string, int>::iterator it;
    std::string measure_outcome_qubits;
    std::string measure_outcome_clbits;
    for (i = 0; i < shots; i++)
    {
        for (j = 0; j < n; j++)
        {
            measure_outcome_qubits += '0';
            measure_outcome_clbits += '0';
        }
        normalize_factor = 1;

        for (j = 0; j < measured_qubits.size(); j++)
            measure_one(j, k/2, H_factor, nVar, nAnci_fourInt, &measure_outcome_qubits);

        // convert measurement outcome of qubits to clbits
        for (j = 0; j < measured_qubits.size(); j++)
        {
            int qIndex = measured_qubits[j];
            int cIndex = measured_qubits_to_clbits[qIndex];
            measure_outcome_clbits[n - 1 - cIndex] = measure_outcome_qubits[n - 1 - qIndex];
        }
        it = state_count.find(measure_outcome_clbits);
        if (it != state_count.end())
            state_count[measure_outcome_clbits] = it->second + 1;
        else
            state_count[measure_outcome_clbits] = 1;
        measure_outcome = measure_outcome_qubits;
        measure_outcome_qubits.clear();
        measure_outcome_clbits.clear();
    }

    Cudd_RecursiveDeref(manager, bigBDD);
    delete[] arrAnci_fourInt;
    delete[] arrAnci_oneInt;
    delete[] permutation;
}

/**Function*************************************************************

  Synopsis    [get statevector string based on BDDs]

  Description []
               
  SideEffects []

  SeeAlso     []

***********************************************************************/
void Simulator::getStatevector()
{
    double oneroot2 = 1 / sqrt(2);
    double H_factor = pow(oneroot2, k);
    double re = 0, im = 0;
    int *assign = new int[n];
    unsigned long long nEntries = pow(2, n);
    int oneEntry;
    long long int_value = 0;
    unsigned long long i, j, h;
    DdNode *tmp;

    for (i = 0; i < n; i++) //initialize assignment
        assign[i] = 0;

    statevector = "[";

    for (i = 0; i < nEntries; i++) // compute every entry
    {
        re = 0;
        im = 0;
        bool isZero = 0;
        for (j = 0; j < n; j++)
        {
            if (measured_qubits_to_clbits[j] != -1)
                if (assign[j] != stoi(measure_outcome.substr(n - 1 - j, 1)))
                {
                    isZero = 1;
                    break;
                }
        }
        
        if (isZero == 0)
        {
            for (j = 0; j < w; j++) // compute every complex value
            {
                int_value = 0;
                for (h = 0; h < r; h++) // compute every integer
                {
                    tmp = Cudd_Eval(manager, All_Bdd[j][h], assign);
                    Cudd_Ref(tmp);
                    oneEntry = !(Cudd_IsComplement(tmp));
                    Cudd_RecursiveDeref(manager, tmp);
                    if (h == r - 1)
                        int_value -= oneEntry * pow(2, h);
                    else
                        int_value += oneEntry * pow(2, h);
                }
                /* translate to re and im */
                re += int_value * cos((double) (w - j - 1)/w * PI);
                im += int_value * sin((double) (w - j - 1)/w * PI);
            }
            re *= H_factor*normalize_factor;
            im *= H_factor*normalize_factor;
        }

        if ((re == 0)&&(im == 0))
            statevector = statevector + "\"0\"";
        else if (re == 0)
            statevector = statevector + "\"" + std::to_string(im) + "i\"";
        else if (im == 0)
            statevector = statevector + "\"" + std::to_string(re) + "\"";
        else
        {
            if (im < 0)
                statevector = statevector + "\"" + std::to_string(re) + std::to_string(im) + "i\"";
            else
                statevector = statevector + "\"" + std::to_string(re) + "+" + std::to_string(im) + "i\"";
        }
        if (i != nEntries - 1)
            statevector = statevector + ", ";
        full_adder_plus_1(n, assign);
    }
    statevector += "]";

    delete[] assign;
}