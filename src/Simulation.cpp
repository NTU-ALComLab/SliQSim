#include "Simulator.h"
#include "util_sim.h"


/**Function*************************************************************

  Synopsis    [parse and simulate the qasm file]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
void Simulator::sim_qasm_file(std::string qasm)
{
    int i;
    std::string inStr;
    int targ, swapA, swapB;
    std::stringstream inFile_ss(qasm);
    while (getline(inFile_ss, inStr))
    {
        inStr = inStr.substr(0, inStr.find("//"));
        if (inStr.find_first_not_of("\t\n ") != std::string::npos)
        {
            std::stringstream inStr_ss(inStr);
            getline(inStr_ss, inStr, ' ');
            if (inStr == "qreg")
            {
                getline(inStr_ss, inStr, '[');
                getline(inStr_ss, inStr, ']');
                n = stoi(inStr); // get the number n here
                manager = Cudd_Init(n, n, CUDD_UNIQUE_SLOTS, CUDD_CACHE_SLOTS, 0);
                int *constants = new int[n];
                for (i = 0; i < n; i++)
                    constants[i] = 0; // TODO: costom initial state
                measured_qubits_to_clbits = new int[n];
                for (i = 0; i < n; i++)
                    measured_qubits_to_clbits[i] = -1;
                init_state(constants);
                delete[] constants;
                Cudd_AutodynEnable(manager, CUDD_REORDER_SYMM_SIFT);
                getline(inStr_ss, inStr, '\n');
            }
            else if (inStr == "creg"){;}
            else if (inStr == "OPENQASM"){;}
            else if (inStr == "include"){;}
            else if (inStr == "measure")
            {
                isMeasure = 1;
                getline(inStr_ss, inStr, '[');
                getline(inStr_ss, inStr, ']');
                int qIndex = stoi(inStr);
                getline(inStr_ss, inStr, '[');
                getline(inStr_ss, inStr, ']');
                int cIndex = stoi(inStr);
                measured_qubits.push_back(qIndex);
                measured_qubits_to_clbits[qIndex] = cIndex;
            }
            else
            {
                if (inStr == "x")
                {
                    getline(inStr_ss, inStr, '[');
                    getline(inStr_ss, inStr, ']');
                    PauliX(stoi(inStr));
                }
                else if (inStr == "y")
                {
                    getline(inStr_ss, inStr, '[');
                    getline(inStr_ss, inStr, ']');
                    PauliY(stoi(inStr));
                }
                else if (inStr == "z")
                {
                    int *iqubit = new int[1];
                    getline(inStr_ss, inStr, '[');
                    getline(inStr_ss, inStr, ']');
                    iqubit[0] = stoi(inStr);
                    PauliZ(iqubit, 1);
                    delete[] iqubit;
                }
                else if (inStr == "h")
                {
                    getline(inStr_ss, inStr, '[');
                    getline(inStr_ss, inStr, ']');
                    Hadamard(stoi(inStr));
                }
                else if (inStr == "s")
                {
                    getline(inStr_ss, inStr, '[');
                    getline(inStr_ss, inStr, ']');
                    int iqubit = stoi(inStr);
                    Phase_shift(2, &iqubit, 1);
                }
                else if (inStr == "sdg")
                {
                    getline(inStr_ss, inStr, '[');
                    getline(inStr_ss, inStr, ']');
                    int iqubit = stoi(inStr);
                    Phase_shift(-2, &iqubit, 1);
                }
                else if (inStr == "t")
                {
                    getline(inStr_ss, inStr, '[');
                    getline(inStr_ss, inStr, ']');
                    int iqubit = stoi(inStr);
                    Phase_shift(4, &iqubit, 1);
                }
                else if (inStr == "tdg")
                {
                    getline(inStr_ss, inStr, '[');
                    getline(inStr_ss, inStr, ']');
                    int iqubit = stoi(inStr);
                    Phase_shift(-4, &iqubit, 1);
                }
                else if (inStr == "rx(pi/2)")
                {
                    getline(inStr_ss, inStr, '[');
                    getline(inStr_ss, inStr, ']');
                    rx_pi_2(stoi(inStr));
                }
                else if (inStr == "ry(pi/2)")
                {
                    getline(inStr_ss, inStr, '[');
                    getline(inStr_ss, inStr, ']');
                    ry_pi_2(stoi(inStr));
                }
                else if (inStr == "cx")
                {
                    std::vector<int> cont(1);
                    std::vector<int> ncont(0);
                    int targ;
                    getline(inStr_ss, inStr, '[');
                    getline(inStr_ss, inStr, ']');
                    cont[0] = stoi(inStr);
                    getline(inStr_ss, inStr, '[');
                    getline(inStr_ss, inStr, ']');
                    targ = stoi(inStr);
                    Toffoli(targ, cont, ncont);
                    cont.clear();
                    ncont.clear();
                }
                else if (inStr == "cz")
                {
                    int *iqubit = new int[2];
                    for (i = 0; i < 2; i++)
                    {
                        getline(inStr_ss, inStr, '[');
                        getline(inStr_ss, inStr, ']');
                        iqubit[i] = stoi(inStr);
                    }
                    PauliZ(iqubit, 2);
                    delete[] iqubit;
                }
                else if (inStr == "swap")
                {
                    std::vector<int> cont(0);
                    for (i = 0; i < 2; i++)
                    {
                        getline(inStr_ss, inStr, '[');
                        getline(inStr_ss, inStr, ']');
                        if (i == 0)
                            swapA = stoi(inStr);
                        else
                            swapB = stoi(inStr);
                    }
                    Fredkin(swapA, swapB, cont);
                    cont.clear();
                }
                else if (inStr == "cswap")
                {
                    std::vector<int> cont(1);
                    for (i = 0; i < 3; i++)
                    {
                        getline(inStr_ss, inStr, '[');
                        getline(inStr_ss, inStr, ']');
                        if (i == 0)
                            cont[i] = stoi(inStr);
                        else if (i == 1)
                            swapA = stoi(inStr);
                        else
                            swapB = stoi(inStr);
                    }
                    Fredkin(swapA, swapB, cont);
                    cont.clear();
                }
                else if (inStr == "ccx" || inStr == "mcx")
                {
                    std::vector<int> cont(0);
                    std::vector<int> ncont(0);
                    int targ;
                    getline(inStr_ss, inStr, '[');
                    while(getline(inStr_ss, inStr, ']'))
                    {
                        cont.push_back(stoi(inStr));
                        getline(inStr_ss, inStr, '[');
                    }
                    targ = cont.back();
                    cont.pop_back();
                    Toffoli(targ, cont, ncont);
                    cont.clear();
                    ncont.clear();
                }
                else
                {
                    gatecount--;
                    std::cerr << std::endl
                            // << "[warning]: Gate \'" << inStr << "\' is not supported in this simulator. The gate is ignored ..." << std::endl;
                            << "[warning]: Syntax \'" << inStr << "\' is not supported in this simulator. The line is ignored ..." << std::endl;
                }
                gatecount++;
                if (manager != NULL)
                    nodecount();
            }
        }
    }
    Cudd_AutodynDisable(manager);
}

/**Function*************************************************************

  Synopsis    [simulate the circuit described by a qasm file]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
void Simulator::sim_qasm(std::string qasm)
{
    int i;
    std::unordered_map<std::string, int>::iterator it;
    sim_qasm_file(qasm); // simulate

    // measure based on simulator type
    if (sim_type == 0) // weak
    {
        if (isMeasure == 0)
        {
            std::string noMeasureBasis;
            for (i = 0; i < n; i++)
                noMeasureBasis += '0';
            state_count[noMeasureBasis] = shots;
        }
        else
            measurement();
        statevector = "null";
    }
    else if (sim_type == 1) // strong
    {
        if (isMeasure == 0)
        {
            std::string noMeasureBasis;
            for (i = 0; i < n; i++)
                noMeasureBasis += '0';
            state_count[noMeasureBasis] = 1;
        }
        else
            measurement();
        getStatevector();
    }

    // write output string based on state_count and statevector
    run_output = "{ \"counts\": { ";
    for (it = state_count.begin(); it != state_count.end(); it++)
    {
        if (std::next(it) == state_count.end())
            run_output = run_output + "\"" + it->first + "\": " + std::to_string(it->second) + " }";
        else
            run_output = run_output + "\"" + it->first + "\": " + std::to_string(it->second) + ", ";
    }
    run_output += (statevector != "null") ? ", \"statevector\": " + statevector + " }" : " }";

    std::cout << run_output << std::endl;
}
