#include "Simulator.h"
#include "util_sim.h"


/**Function*************************************************************

  Synopsis    [initialize state vector by a basis state]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
void Simulator::init_state(int *constants)
{
    DdNode *var, *tmp;
    All_Bdd = new DdNode **[w];
    for (int i = 0; i < w; i++)
        All_Bdd[i] = new DdNode *[r];

    for (int i = 0; i < r; i++)
    {
        if (i == 0)
        {
            for (int j = 0; j < w - 1; j++)
            {
                All_Bdd[j][i] = Cudd_Not(Cudd_ReadOne(manager));
                Cudd_Ref(All_Bdd[j][i]);
            }
            All_Bdd[w - 1][i] = Cudd_ReadOne(manager);
            Cudd_Ref(All_Bdd[w - 1][i]);
            for (int j = n - 1; j >= 0; j--)
            {
                var = Cudd_bddIthVar(manager, j);
                if (constants[j] == 0)
                    tmp = Cudd_bddAnd(manager, Cudd_Not(var), All_Bdd[w - 1][i]);
                else
                    tmp = Cudd_bddAnd(manager, var, All_Bdd[w - 1][i]);
                Cudd_Ref(tmp);
                Cudd_RecursiveDeref(manager, All_Bdd[w - 1][i]);
                All_Bdd[w - 1][i] = tmp;
            }
        }
        else
        {
            for (int j = 0; j < w; j++)
            {
                All_Bdd[j][i] = Cudd_Not(Cudd_ReadOne(manager));
                Cudd_Ref(All_Bdd[j][i]);
            }
        }
    }
}

/**Function*************************************************************

  Synopsis    [allocate new BDDs for each integer vector]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
void Simulator::alloc_BDD(DdNode ***Bdd, bool extend)
{
    DdNode *tmp;

    DdNode ***W = new DdNode **[w];
    for (int i = 0; i < w; i++)
        W[i] = new DdNode *[r];

    for (int i = 0; i < r - inc; i++)
        for (int j = 0; j < w; j++)
            W[j][i] = Bdd[j][i];

    for (int i = 0; i < w; i++)
        delete[] Bdd[i];

    for (int i = 0; i < w; i++)
        Bdd[i] = W[i];

    if (extend)
    {
        for (int i = r - inc; i < r; i++)
        {
            for (int j = 0; j < w; j++)
            {
                Bdd[j][i] = Cudd_ReadOne(manager);
                Cudd_Ref(Bdd[j][i]);
                tmp = Cudd_bddAnd(manager, Bdd[j][r - inc - 1], Bdd[j][i]);
                Cudd_Ref(tmp);
                Cudd_RecursiveDeref(manager, Bdd[j][i]);
                Bdd[j][i] = tmp;
            }
        }
    }
}

/**Function*************************************************************

  Synopsis    [Drop LSB and shift right by 1 bit]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
void Simulator::dropLSB(DdNode ***Bdd)
{
    DdNode *tmp;

    for (int i = 0; i < w; i++)
    {
        Cudd_RecursiveDeref(manager, Bdd[i][0]); // drop LSB
        // right shift
        for (int j = 1; j < r; j++)
        {
            Bdd[i][j - 1] = Bdd[i][j];
        }
        // sign extension
        Bdd[i][r - 1] = Cudd_ReadOne(manager);
        Cudd_Ref(Bdd[i][r - 1]);
        tmp = Cudd_bddAnd(manager, Bdd[i][r - 2], Bdd[i][r - 1]);
        Cudd_Ref(tmp);
        Cudd_RecursiveDeref(manager, Bdd[i][r - 1]);
        Bdd[i][r - 1] = tmp;
    }
}

/**Function*************************************************************

  Synopsis    [detect overflow in integer vectors]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
int Simulator::overflow3(DdNode *g, DdNode *h, DdNode *crin)
{
    DdNode *tmp, *dd1, *dd2;
    int overflow;

    dd1 = Cudd_bddXor(manager, g, crin);
    Cudd_Ref(dd1);

    dd2 = Cudd_bddXnor(manager, g, h);
    Cudd_Ref(dd2);

    tmp = Cudd_bddAnd(manager, dd1, dd2);
    Cudd_Ref(tmp);
    Cudd_RecursiveDeref(manager, dd1);
    Cudd_RecursiveDeref(manager, dd2);

    if (Cudd_CountPathsToNonZero(tmp))
        overflow = 1;
    else
        overflow = 0;
    Cudd_RecursiveDeref(manager, tmp);

    return overflow;
}

/**Function*************************************************************

  Synopsis    [detect overflow in integer vectors -- for the case that h is 0]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
int Simulator::overflow2(DdNode *g, DdNode *crin){
    DdNode *tmp;
    int overflow;

    tmp = Cudd_bddAnd(manager, Cudd_Not(g), crin);
    Cudd_Ref(tmp);

    if (Cudd_CountPathsToNonZero(tmp))
        overflow = 1;
    else
        overflow = 0;
    Cudd_RecursiveDeref(manager, tmp);

    return overflow;
}

/**Function*************************************************************

  Synopsis    [decode and print each entry of the state vector]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
void Simulator::decode_entries()
{
    double oneroot2 = 1 / sqrt(2);
    double H_factor = pow(oneroot2, k);
    double re = 0, im = 0;
    int *assign = new int[n];
    unsigned long long nEntries = pow(2, n);
    int oneEntry;
    long long int_value = 0;
    DdNode *tmp;

    for (int i = 0; i < n; i++) //initialize assignment
        assign[i] = 0;

    std::cout << "Amplitudes of the Computational Basis States:" << std::endl;

    for (unsigned long long i = 0; i < nEntries; i++) // compute every entry
    {
        re = 0;
        im = 0;
        for (int j = 0; j < w; j++) // compute every complex value
        {
            int_value = 0;
            for (int h = 0; h < r; h++) // compute every integer
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
        re *= H_factor;
        im *= H_factor;
        std::cout << i << ": " << re << " + " << im << "i" << std::endl;
        full_adder_plus_1(n, assign);
    }

    delete[] assign;
}

/**Function*************************************************************

  Synopsis    [reorder BDDs]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
void Simulator::reorder()
{
    int reorder_signal = Cudd_ReduceHeap(manager, CUDD_REORDER_SYMM_SIFT, 0);
    if (!reorder_signal)
        std::cout << "reorder fails" << std::endl;
}

/**Function*************************************************************

  Synopsis    [update max #nodes]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
void Simulator::nodecount()
{
    unsigned long NodeCount_new = Cudd_ReadNodeCount(manager);
    if (NodeCount_new > NodeCount)
         NodeCount = NodeCount_new;
}

/**Function*************************************************************

  Synopsis    [print statistics]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
void Simulator::print_info(double runtime, size_t memPeak)
{
    std::cout << "  Runtime: " << runtime << " seconds" << std::endl;
    std::cout << "  Peak memory usage: " << memPeak << " bytes" << std::endl; //unit in bytes
    std::cout << "  #Applied gates: " << gatecount << std::endl;
    std::cout << "  Max #nodes: " << NodeCount << std::endl;
    std::cout << "  Integer bit size: " << r << std::endl;
    std::cout << "  Accuracy loss: " << error << std::endl;
    // std::cout << "  #Integers: " << w << std::endl;

    // std::unordered_map<std::string, int>::iterator it;
    // std::cout << "  Measurement: " << std::endl;
    // for(it = state_count.begin(); it != state_count.end(); it++)
    //     std::cout << "      " << it->first << ": " << it->second << std::endl;
}

/**Function*************************************************************

  Synopsis    [parse an infix Boolean function to postfix]

  Description [use shunting yard algorithm]

  SideEffects []

  SeeAlso     []

***********************************************************************/

std::vector<std::string> Simulator::boolean_parser(std::string& inStr)
{
    inStr = std::regex_replace(inStr, std::regex("!"), " not ");
    inStr = std::regex_replace(inStr, std::regex("~"), " not ");
    inStr = std::regex_replace(inStr, std::regex("\\*"), " and ");
    inStr = std::regex_replace(inStr, std::regex("&"), " and ");
    inStr = std::regex_replace(inStr, std::regex("\\+"), " or ");
    inStr = std::regex_replace(inStr, std::regex("\\|"), " or ");
    inStr = std::regex_replace(inStr, std::regex("\\^"), " xor ");
    inStr = std::regex_replace(inStr, std::regex("\\("), " ( ");
    inStr = std::regex_replace(inStr, std::regex("\\)"), " ) ");
    
    std::vector<std::string> result;
    std::map<std::string, int> operator_priority = { {"not", 4},
                                                     {"xor", 3},
                                                     {"and", 2},
                                                     {"or",  1}, };
        
    std::stringstream inStr_ss(inStr);
    std::stack<std::string> waiting_operators;
    int nLeftTuple = 0;
    while(getline(inStr_ss, inStr, ' '))
    {
        if (inStr == "") continue;
                
        if (inStr == "(")
        {
            waiting_operators.push(inStr);
            nLeftTuple++;
        }
        else if (inStr == ")")
        {
            assert(nLeftTuple > 0);
            while (waiting_operators.top() != "(") 
            {
                result.emplace_back(waiting_operators.top());
                waiting_operators.pop();
            }
            waiting_operators.pop();
            nLeftTuple--;
        }
        else if (operator_priority.count(inStr) != 0)
        {
            while (!waiting_operators.empty() && waiting_operators.top() != "(" && operator_priority[inStr] <= operator_priority[waiting_operators.top()])
            {
                result.emplace_back(waiting_operators.top());
                waiting_operators.pop();
            }
            waiting_operators.push(inStr);
        }
        else if (inStr.find( "[" ) != std::string::npos) 
        {
            std::stringstream inStr_sss(inStr);
            getline(inStr_sss, inStr, '[');
            getline(inStr_sss, inStr, ']');
            assert(inStr != "");
            result.emplace_back(inStr);
        }
        else  // self-defined variables
        {
            result.emplace_back(inStr);
        }
    }
    assert(nLeftTuple == 0);
    
    while (!waiting_operators.empty())
    {
        result.emplace_back(waiting_operators.top());
        waiting_operators.pop();
    }
    return result;
}

/**Function*************************************************************

  Synopsis    [compare two node lists and return a node representing 
               whether two node lists represents the same integers ]

  Description []
               
  SideEffects [returned node already referenced]

  SeeAlso     []

***********************************************************************/
DdNode* Simulator::node_equiv(std::vector<DdNode*>& int_1, std::vector<DdNode*>& int_2)
{
    DdNode* result = Cudd_ReadOne(manager);
    Cudd_Ref(result);
    
    for (int i = 0; i < int_1.size(); ++i)
    {
        DdNode* tmp = Cudd_bddXnor(manager, int_1[i], int_2[i]);
        Cudd_Ref(tmp);
        DdNode* new_result = Cudd_bddAnd(manager, tmp, result);
        Cudd_Ref(new_result);
        Cudd_RecursiveDeref(manager, tmp);
        Cudd_RecursiveDeref(manager, result);
        result = new_result;
    }
    
    return result;
}

/**Function*************************************************************

  Synopsis    [compare two node lists and return a node representing 
               whether the first node list represents a larger integer ]

  Description []
               
  SideEffects [returned node already referenced]

  SeeAlso     []

***********************************************************************/
DdNode* Simulator::node_larger(std::vector<DdNode*>& int_1, std::vector<DdNode*>& int_2)
{
    DdNode* result = Cudd_Not(Cudd_ReadOne(manager));
    Cudd_Ref(result);
    
    DdNode* even = Cudd_ReadOne(manager);
    Cudd_Ref(even);
    
    for (int i = 0; i < int_1.size(); ++i)
    {
        DdNode* tmp1 = Cudd_bddAnd(manager, int_1[i], Cudd_Not(int_2[i]));
        Cudd_Ref(tmp1);
        DdNode* new_larger = Cudd_bddAnd(manager, tmp1, even);
        Cudd_Ref(new_larger);
        DdNode* new_result = Cudd_bddOr(manager, result, new_larger);
        Cudd_Ref(new_result);
        Cudd_RecursiveDeref(manager, tmp1);
        Cudd_RecursiveDeref(manager, new_larger);
        Cudd_RecursiveDeref(manager, result);
        result = new_result;
        
        DdNode* tmp2 = Cudd_bddXnor(manager, int_1[i], int_2[i]);
        Cudd_Ref(tmp2);
        DdNode* new_even = Cudd_bddAnd(manager, tmp2, even);
        Cudd_Ref(new_even);
        Cudd_RecursiveDeref(manager, tmp2);
        Cudd_RecursiveDeref(manager, even);
        even = new_even;
    }
    Cudd_RecursiveDeref(manager, even);
    
    return result;
}

/**Function*************************************************************

  Synopsis    [build a DdNode* type from a function    ]

  Description []

  SideEffects [returned node already referenced]

  SeeAlso     []

***********************************************************************/


DdNode* Simulator::func2node(std::vector<std::string>& func)
{
    assert(!func.empty());
    
    std::stack<DdNode*> waiting;
    for (int i = 0; i < func.size(); ++i)
    {
        if (func[i] == "not")
        {
            assert(waiting.size() >= 1);
            DdNode* temp = waiting.top();
            waiting.pop();
            
            DdNode* new_item = Cudd_Not(temp);
            Cudd_Ref(new_item);
            waiting.push(new_item);
            
            Cudd_RecursiveDeref(manager, temp);
        }
        else if (func[i] == "and" || func[i] == "or" || func[i] == "xor")
        {
            assert(waiting.size() >= 2);
            DdNode* temp_1 = waiting.top();
            waiting.pop();
            DdNode* temp_2 = waiting.top();
            waiting.pop();
             
            DdNode* new_item;              
                 if (func[i] == "and") new_item = Cudd_bddAnd(manager, temp_1, temp_2);
            else if (func[i] == "or")  new_item = Cudd_bddOr (manager, temp_1, temp_2);
            else if (func[i] == "xor") new_item = Cudd_bddXor(manager, temp_1, temp_2);
            Cudd_Ref(new_item);
            waiting.push(new_item);
            
            Cudd_RecursiveDeref(manager, temp_1);
            Cudd_RecursiveDeref(manager, temp_2);
        }
        else
        {
            if (func[i].find_first_not_of( "0123456789" ) == std::string::npos) 
            {
                int ith_var = stoi(func[i]);
                waiting.push(Cudd_bddIthVar(manager, ith_var));
                Cudd_Ref(waiting.top());
            }
            else 
            {
                std::unordered_map<std::string, DdNode*>::iterator iter = defined_var.find(func[i]);
                if (iter == defined_var.end()) 
                {
                    std::cerr << std::endl
                            << "[warning]: Variable \'" << func[i] << "\' is not defined. The variable is taken as constant 0 ..." << std::endl;
                    waiting.push(Cudd_Not(Cudd_ReadOne(manager)));
                    Cudd_Ref(waiting.top());
                }
                else
                {
                    waiting.push(iter->second);
                    Cudd_Ref(waiting.top());
                }
            }
        }
    }
    assert(waiting.size() == 1);
    DdNode *result = waiting.top();
    waiting.pop();
    return result;
}