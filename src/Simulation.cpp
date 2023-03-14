#include "Simulator.h"
#include "util_sim.h"


/**Function*************************************************************

  Synopsis    [Initailize simulator]

  Description [This function will set #qubits n, construct initial state, and enable dynamic reordering]

  SideEffects []

  SeeAlso     []

***********************************************************************/
void Simulator::init_simulator(int nQubits)
{
    n = nQubits; // set the number n here
    manager = Cudd_Init(n, n, CUDD_UNIQUE_SLOTS, CUDD_CACHE_SLOTS, 0);
    int *constants = new int[n];
    for (int i = 0; i < n; i++)
        constants[i] = 0; // TODO: costom initial state
    measured_qubits_to_clbits = std::vector<std::vector<int>>(n, std::vector<int>(0));
    init_state(constants);
    delete[] constants;
    if (isReorder) Cudd_AutodynEnable(manager, CUDD_REORDER_SYMM_SIFT);
}


/**Function*************************************************************

  Synopsis    [parse and simulate the qasm file]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
void Simulator::sim_qasm_file(std::string qasm)
{
    std::string inStr;
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
                init_simulator(stoi(inStr));
            }
            else if (inStr == "creg")
            {
                getline(inStr_ss, inStr, '[');
                getline(inStr_ss, inStr, ']');
                nClbits = stoi(inStr);                
            }
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
                measure(qIndex, cIndex);
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
                    std::vector<int> iqubit(1);
                    getline(inStr_ss, inStr, '[');
                    getline(inStr_ss, inStr, ']');
                    iqubit[0] = stoi(inStr);
                    PauliZ(iqubit);
                    iqubit.clear();
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
                    Phase_shift(2, iqubit);
                }
                else if (inStr == "sdg")
                {
                    getline(inStr_ss, inStr, '[');
                    getline(inStr_ss, inStr, ']');
                    int iqubit = stoi(inStr);
                    Phase_shift_dagger(-2, iqubit);
                }
                else if (inStr == "t")
                {
                    getline(inStr_ss, inStr, '[');
                    getline(inStr_ss, inStr, ']');
                    int iqubit = stoi(inStr);
                    Phase_shift(4, iqubit);
                }
                else if (inStr == "tdg")
                {
                    getline(inStr_ss, inStr, '[');
                    getline(inStr_ss, inStr, ']');
                    int iqubit = stoi(inStr);
                    Phase_shift_dagger(-4, iqubit);
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
                    std::vector<int> iqubit(2);
                    for (int i = 0; i < 2; i++)
                    {
                        getline(inStr_ss, inStr, '[');
                        getline(inStr_ss, inStr, ']');
                        iqubit[i] = stoi(inStr);
                    }
                    PauliZ(iqubit);
                    iqubit.clear();
                }
                else if (inStr == "swap")
                {
                    int swapA, swapB;
                    std::vector<int> cont(0);
                    for (int i = 0; i < 2; i++)
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
                    int swapA, swapB;
                    std::vector<int> cont(1);
                    for (int i = 0; i < 3; i++)
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
                    std::cerr << std::endl
                            // << "[warning]: Gate \'" << inStr << "\' is not supported in this simulator. The gate is ignored ..." << std::endl;
                            << "[warning]: Syntax \'" << inStr << "\' is not supported in this simulator. The line is ignored ..." << std::endl;
                }
            }
        }
    }
    if (isReorder) Cudd_AutodynDisable(manager);
}

/**Function*************************************************************

  Synopsis    [simulate the circuit described by a qasm file]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
void Simulator::sim_qasm(std::string qasm)
{
    sim_qasm_file(qasm); // simulate
    
    if (sim_type == 0 && isMeasure == 0)
    {
        std::cout << "Error: no measurement detected. Cannot do sampling.\n" << std::flush;
        assert(sim_type != 0 || isMeasure != 0);
    }
    if (sim_type == 1)
    {
        if (isMeasure == 1)
        {
            std::cout << "Warning: measurement detected. The final statevector will collapse based on the measurement outcome.\n" << std::flush;
            if (shots != 1)
            {
                shots = 1;
                std::cout << "Warning: shot number is limited to 1 in all_amplitude mode.\n" << std::flush;
            }
        }    
        else 
        {
            if (shots != 1)
            {
                std::cout << "Warning: no measurement detected. The --shots argument is ignored.\n" << std::flush;
            }
        }  
    }

    // measure based on simulator type
    if (sim_type == 0) // sampling mode
    {
        measurement();
    }
    else if (sim_type == 1) // all_amplitude mode
    {
        if (isMeasure == 1)
        {
            measurement();
        }
        getStatevector();
    }
    print_results();
}



/**Function*************************************************************

  Synopsis    [print state vector and distribution of sampled outcomes]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
void Simulator::print_results()
{
    // write output string based on state_count and statevector
    std::unordered_map<std::string, int>::iterator it;
    
    run_output = "{";
    if (state_count.begin() != state_count.end()){
        run_output += "\"counts\": { ";
        for (it = state_count.begin(); it != state_count.end(); it++)
        {
            if (std::next(it) == state_count.end())
                run_output = run_output + "\"" + it->first + "\": " + std::to_string(it->second);
            else
                run_output = run_output + "\"" + it->first + "\": " + std::to_string(it->second) + ", ";
        }
        run_output += " }";
        run_output += (statevector != "null") ? ", " : ""; 
    }    

    run_output += (statevector != "null") ? "\"statevector\": " + statevector + " }" : " }";
    //return;
    std::cout << run_output << std::endl;
}