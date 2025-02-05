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
    int skip_level;
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
            for (int i = 0; i < w; i++)
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
            for (int i = 0; i < nVar; i++)
                assign[i] = 0;
            /* TODO: BDD to truth table */
            for (int i = 0; i < w; i++) //compute each complex value
            {
                int_value = 0;
                for (int j = 0; j < r; j++) //compute each integer
                {
                    tmp = Cudd_Eval(manager, child, assign);
                    Cudd_Ref(tmp);
                    oneEntry = !(Cudd_IsComplement(tmp));
                    Cudd_RecursiveDeref(manager, tmp);
                    if (j == r - 1)
                        int_value -= oneEntry * pow(2, j + shift - kd2);
                    else
                        int_value += oneEntry * pow(2, j + shift - kd2);
                    full_adder_plus_1_start(nVar, assign, n + nAnci_fourInt);
                }
                /* translate to re and im */
                re += int_value * cos((double) (w - i - 1)/w * PI);
                im += int_value * sin((double) (w - i - 1)/w * PI);
                full_adder_plus_1_start(nVar, assign, n);
                for (int j = n + nAnci_fourInt; j < nVar; j++) // reset array: not necessary but straightforward
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

  Synopsis    [create bigBDD for measurement]

  Description []
               
  SideEffects []

  SeeAlso     []

***********************************************************************/
void Simulator::create_bigBDD()
{
    double oneroot2 = 1 / sqrt(2);
    double H_factor = pow(oneroot2, k%2);
    int nAnci_oneInt = ceil(log(r) / log(2)), nAnci_fourInt = ceil(log(w) / log(2)), nAnci = nAnci_oneInt + nAnci_fourInt, nnAnci_fourInt = n + nAnci_fourInt, nVar = n + nAnci;
    DdNode *tmp1, *tmp2, *tmp3;
    
    int *arrAnci_fourInt = new int[nAnci_fourInt];
    for (int i = 0; i < nAnci_fourInt; i++)
        arrAnci_fourInt[i] = 0;
    int *arrAnci_oneInt = new int[nAnci_oneInt];
    for (int i = 0; i < nAnci_oneInt; i++)
        arrAnci_oneInt[i] = 0;
    
    bigBDD = Cudd_Not(Cudd_ReadOne(manager));
    Cudd_Ref(bigBDD);
    for (int i = 0; i < w; i++)
    {
        tmp3 = Cudd_Not(Cudd_ReadOne(manager));
        Cudd_Ref(tmp3);
        for (int j = 0; j < r; j++)
        {
            tmp1 = Cudd_ReadOne(manager);
            Cudd_Ref(tmp1);
            for (int h = n + nAnci - 1; h >= nnAnci_fourInt; h--)
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
        for (int j = nnAnci_fourInt - 1; j >= n; j--)
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
        for (int j = 0; j < nAnci_oneInt; j++) // reset array: not necessary but straightforward
            arrAnci_oneInt[j] = 0;
    }
    nodecount();
    
    delete[] arrAnci_fourInt;
    delete[] arrAnci_oneInt;
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
    
    // move measured qubits to the top
    int *permutation = new int[nVar];
    int indCount1 = 0;
    int indCount0 = 0;  // number of measured qubits
    for (int i = 0; i < n; i++){
        if (!measured_qubits_to_clbits[i].empty())
        {
            indCount0++;
        }
    }
    for (int i = 0; i < n; i++)
    {
        int index = Cudd_ReadInvPerm(manager, i);
        if (!measured_qubits_to_clbits[index].empty())
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
    for (int i = n; i < nVar; i++)
        permutation[i] = Cudd_ReadInvPerm(manager, i);      
    int dum = Cudd_ShuffleHeap(manager, permutation);
    Node_Table.clear();    // need to re-calculate
    nodecount();
    
    std::unordered_map<std::string, int>::iterator it;
    std::string measure_outcome_qubits;
    std::string measure_outcome_clbits;
    for (int i = 0; i < shots; i++)
    {
        for (int j = 0; j < n; j++)
        {
            measure_outcome_qubits += '0';
        }
        for (int j = 0; j < nClbits; j++) 
        {
            measure_outcome_clbits += '0';
        }
        normalize_factor = 1;

        for (int j = 0; j < indCount1; j++)  // measure for the (top) j^th level variable
        {
            measure_one(j, k/2, H_factor, nVar, nAnci_fourInt, &measure_outcome_qubits);            
        }
        
        // convert measurement outcome of qubits to clbits
        for (int qIndex = 0; qIndex < n; qIndex++)
        {
            for (int cIndex : measured_qubits_to_clbits[qIndex])
            {
                measure_outcome_clbits[nClbits - 1 - cIndex] = measure_outcome_qubits[n - 1 - qIndex];
                // the order is reversed
            }
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

    delete[] permutation;
}

/**Function*************************************************************

  Synopsis    [measurement with self-defined functions]

  Description []
               
  SideEffects []

  SeeAlso     []
  
***********************************************************************/
void Simulator::measurement_obs(std::string obsfile)
{
    normalize_factor = 1;
    
    obsfile = std::regex_replace(obsfile, std::regex("\r\n"), "\n");
    defined_var["0"]     = Cudd_Not(Cudd_ReadOne(manager));        Cudd_Ref(defined_var["0"]);
    defined_var["false"] = Cudd_Not(Cudd_ReadOne(manager));        Cudd_Ref(defined_var["false"]);
    defined_var["1"]     = Cudd_ReadOne(manager);                  Cudd_Ref(defined_var["1"]);
    defined_var["true"]  = Cudd_ReadOne(manager);                  Cudd_Ref(defined_var["true"]);
    
    std::string line;
    std::stringstream inFile_ss(obsfile);
           
    while (getline(inFile_ss, line))
    {   
        std::string raw_line = line;
        
        line = line.substr(0, line.find("//"));
        line = line.substr(0, line.find(";"));  
        std::transform(line.begin(), line.end(), line.begin(),
                       [](unsigned char c){ return std::tolower(c); }); 
        if (line.find_first_not_of("\t\n\r ") == std::string::npos)
            continue;
            
        std::string word;
        std::stringstream line_ss(line);
        std::vector<std::string> words;
        while (getline(line_ss, word, ' '))
            words.push_back(word);
        
        if (words.size() == 0)
            continue;
        if (words[0] == KW_WS)
        {
            std::string all_lines = (std::string)KW_WS + "\n";
            double value = 0;
            while(getline(inFile_ss, line))
            {
                raw_line = line;
                line = line.substr(0, line.find("//"));
                line = line.substr(0, line.find(";"));  
                std::transform(line.begin(), line.end(), line.begin(),
                               [](unsigned char c){ return std::tolower(c); }); 
                if (line.find_first_not_of("\t\n\r ") == std::string::npos)
                    continue;
                    
                std::stringstream line_ss(line);
                std::vector<std::string> words;
                while (getline(line_ss, word, ' '))
                    if (word != "")
                        words.push_back(word);
                    
                if (words.size() == 0)
                    continue;
                if (words[0] == KW_EWS)
                    break;
                    
                char *p;       
                double weight = strtof(words[0].c_str(), &p);
                if (*p != 0)
                {
                    std::cerr << std::endl << "[warning]: In line \"" << raw_line << "\":" << std::endl;
                    std::cerr << "Syntax \'" << words[0] << "\' is not supported inside " << KW_WS << ". The line is ignored ..." << std::endl;
                    continue;
                }      
                
                words.erase(words.begin());
                int result = handle_property(raw_line, words, true);
                if (result == 0)
                {
                    all_lines += std::to_string(weight) + " " + std::get<0>(property.back()) + "\n";
                    value += weight * std::stold(std::get<1>(property.back()));
                    property.pop_back();
                }
            }
            
            all_lines = all_lines + KW_EWS;
            if (condition_stack == "")
            {
                std::ostringstream oss;
                oss << value;
                property.push_back(std::make_tuple(all_lines, oss.str()));
            }
            else
            {
                property.push_back(std::make_tuple(getCondStr(all_lines), getCondRslt(value)));
                condition_stack = "";
            }
        }
        else if (words[0] == KW_EWS)
        {
            std::cerr << std::endl << "[warning]: In line \"" << raw_line << "\":" << std::endl;
            std::cerr << "Missing a leading \"" << KW_WS << "\". The line is ignored ..." << std::endl;
            continue;
        }
        else
        {
            std::string op_name = words[0];
            int result = handle_property(raw_line, words, false);
            
            if (condition_stack != "" && keyword_condition.find(op_name) == keyword_condition.end())
            {
                if (result == 1)
                {
                    std::cerr << "The previous \"" << condition_stack << "\" is ignored ..." << std::endl;
                }
                else if (keyword_value.find(op_name) != keyword_value.end())
                {
                    std::string line = std::get<0>(property.back());
                    double value = std::stold(std::get<1>(property.back()));
                    property.pop_back();
                    property.push_back(std::make_tuple(getCondStr(line), getCondRslt(value)));
                }
                else
                {  
                    std::cerr << std::endl << "[warning]: In line \"" << raw_line << "\":" << std::endl;
                    std::cerr << op_name << " not supported for \'" << condition_stack << "\'. The previous \"" << condition_stack << "\" is ignored ..." << std::endl;
                }
                condition_stack = "";
            }
        }
    }
    if (condition_stack != "")
    {
        std::cerr << std::endl << "[warning]: A previous \"" << condition_stack << "\" not used." << std::endl;
    }
    
    // show results
    std::cout << std::endl;
    for (auto& item : property)
    {
        //std::cout << std::setw(30) << std::left << "\"" + std::get<0>(item) + "\":" << std::get<1>(item) << std::endl;
        std::cout << "\"" << std::get<0>(item) + "\":\n\t" << std::get<1>(item) << std::endl;
    }
}

/**Function*************************************************************

  Synopsis    [handle a line of a self-defined property.
               return 0 if succeed]

  Description []
               
  SideEffects []

  SeeAlso     []                  
                                 
***********************************************************************/
int Simulator::handle_property(std::string& raw_line, std::vector<std::string>& words, bool is_in_ws)
{
    std::string op_name = words[0];
    words.erase(words.begin());
    
    if (is_in_ws && keyword_value.find(op_name) == keyword_value.end())
    {
        std::cerr << std::endl << "[warning]: In line \"" << raw_line << "\":" << std::endl;
        std::cerr << op_name << " not supported inside \'" << KW_WS << "\'. The line is ignored ..." << std::endl;
        return 1;
    }
    
    if (op_name == KW_ASSIGN)                                        // ============================ KW_ASSIGN
    {
        if (words.size() < 2)
        {
            std::cerr << std::endl << "[warning]: In line \"" << raw_line << "\":" << std::endl;
            std::cerr << "Not enough arguments. The line is ignored ..." << std::endl;
            return 1;
        }
        
        std::string var_name = words[0];
        if (defined_var.find(var_name) != defined_var.end()) {
            std::cerr << std::endl << "[warning]: In line \"" << raw_line << "\":" << std::endl;
            std::cerr << var_name << " already defined. The line is ignored ..." << std::endl;
            return 1;
        }
        words.erase(words.begin());
        
        std::string op_name = words[0];
        if (keyword_truth.find(op_name) == keyword_truth.end()) {
            std::cerr << std::endl << "[warning]: In line \"" << raw_line << "\":" << std::endl;
            std::cerr << op_name << " not supported for \'" << KW_ASSIGN << "\'. The line is ignored ..." << std::endl;
            return 1;
        }
        words.erase(words.begin());
        
        DdNode* function = get_property(op_name, words);
        if (function != nullptr) 
            defined_var[var_name] = function;
        return 0;
    }
    else if (op_name == KW_DIST)                                     // ============================ KW_DIST
    {
        // parse
        std::unordered_map<int, int> qubits;
        for (std::string& word : words)
        {
            std::string tmp;
            std::stringstream word_ss(word);
            getline(word_ss, tmp, '[');
            getline(word_ss, tmp, ']');
            
            char *p;
            int index = strtol(tmp.c_str(), &p, 10);
            if (*p != 0 || tmp.length() == 0 || *p == tmp[0])
            {
                std::cerr << std::endl << "[warning]: In line \"" << raw_line << "\":" << std::endl;
                std::cerr << "Syntax \'" << word << "\' is not specifying a qubit. The line is ignored ..." << std::endl;
                return 1;
            }       
            if (qubits.find(index) != qubits.end())
            {
                std::cerr << std::endl << "[warning]: In line \"" << raw_line << "\":" << std::endl;
                std::cerr << "Qubit \'" << index << "\' is doubly defined. The line is ignored ..." << std::endl;
                return 1;
            }
            qubits[index] = qubits.size();
        }          
        if (qubits.size() == 0)
        {
            std::cerr << std::endl << "[warning]: In line \"" << raw_line << "\":" << std::endl;
            std::cerr << "No qubits are specified. The line is ignored ..." << std::endl;
            return 1;
        }
                
        // backup and reorder
        double oneroot2 = 1 / sqrt(2);
        double H_factor = pow(oneroot2, k%2);
        int nAnci_oneInt = ceil(log(r) / log(2)), nAnci_fourInt = ceil(log(w) / log(2)), nAnci = nAnci_oneInt + nAnci_fourInt, 
            nnAnci_fourInt = n + nAnci_fourInt, nVar = n + nAnci;
        
        int *permutation = new int[nVar];
        int indCount1 = 0;
        int indCount0 = qubits.size();  // number of measured qubits
        for (int i = 0; i < n; i++)
        {
            int index = Cudd_ReadInvPerm(manager, i);
            if (qubits.find(index) != qubits.end())
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
        for (int i = n; i < nVar; i++)
            permutation[i] = Cudd_ReadInvPerm(manager, i);
        int dum = Cudd_ShuffleHeap(manager, permutation);
        Node_Table.clear();    // need to re-calculate
        nodecount();
        
        // iterate
        int *assign = new int[qubits.size()];
        unsigned long long nEntries = pow(2, qubits.size());      // should timeout before overflow occurs in "nEntries"
        for (int i = 0; i < qubits.size(); ++i)                   // initialize assignment
            assign[i] = 0;
    
        std::map<std::string, long double> probs;
        for (unsigned long long i = 0; i < nEntries; i++)         // compute every entry
        {
            std::string real_assign(qubits.size(), '2');
            DdNode* root = bigBDD;
            Cudd_Ref(root);
            
            DdNode* prev = bigBDD;
            int prev_decision = 0;
            
            int current_position = Cudd_ReadPerm(manager, Cudd_NodeReadIndex(root));
            for (int i = 0; i < qubits.size(); ++i)
            {
                real_assign[ qubits[Cudd_ReadInvPerm(manager, i)] ] = assign[i] + '0';
                if ( current_position <= i )
                {
                    DdNode* child = Cudd_Child(manager, root, assign[i]);
                    Cudd_Ref(child);
                    Cudd_RecursiveDeref(manager, prev);
                    prev = root;
                    prev_decision = assign[i];
                    root = child;
                    current_position = Cudd_ReadPerm(manager, Cudd_NodeReadIndex(root));
                }
            }  
        
            int prev_position = Cudd_ReadPerm(manager, Cudd_NodeReadIndex(prev));
            long double value = measure_probability(prev, k/2, nVar, nAnci_fourInt, prev_decision) * H_factor * H_factor / pow(2, qubits.size() - 1 - prev_position);             
            probs[real_assign] = value;
             
            Cudd_RecursiveDeref(manager, root);
            Cudd_RecursiveDeref(manager, prev);
            full_adder_plus_1(n, assign);
        }
        delete[] assign;
        
        // store the value
        std::string result = "";
        for (auto item = probs.begin(); item != probs.end(); ++item)
        {
            std::ostringstream oss;
            oss << item->second;
            result += item->first + ": " + oss.str();
            if (item != std::prev(probs.end()))
                result += ", ";
        }
        std::string line = op_name;
        for (std::string& word : words)
            line = line + " " + word;
        property.push_back(std::make_tuple(line, result));
        return 0;
    }
    else if (op_name == KW_AMP)                                      // ============================ KW_AMP
    {
        // parse
        if (words.size() == 0 || words[0].length() != n)
        {
            std::cerr << std::endl << "[warning]: In line \"" << raw_line << "\":" << std::endl;
            std::cerr << "Illegal computational basis. The line is ignored ..." << std::endl;
            return 1;
        }
        
        int *assign = new int[n];
        for (int i = 0; i < n; i++)
        {
            assign[n - 1 - i] = words[0][i] - '0';
            if (words[0][i] != '0' && words[0][i] != '1')
            {
                std::cerr << std::endl << "[warning]: In line \"" << raw_line << "\":" << std::endl;
                std::cerr << "Illegal computational basis. The line is ignored ..." << std::endl;
                delete[] assign;
                return 1;
            }
        }
        
        // get value
        std::string amplitude_string = get_amplitude_string(assign);
        std::string line = op_name + " " + words[0];
        property.push_back(std::make_tuple(line, amplitude_string));
        return 0;
    }
    else if (op_name == KW_EXPT)                             // ============================ KW_EXPT
    {
        // parse  
        enum {EXPT, M0, M1};
        std::vector<char> basis;
        std::vector<int> type;
        for (int i = 0; i < words.back().length(); ++i)
        {
            if (words.back()[i] != 'i' && words.back()[i] != 'x' && words.back()[i] != 'y' && words.back()[i] != 'z')
            {
                std::cerr << std::endl << "[warning]: In line \"" << raw_line << "\":" << std::endl;
                std::cerr << "Unsupported Pauli basis \'" << words.back()[i] << "\'. This line is ignored..." << std::endl;
                return 1;
            }
            
            basis.push_back(words.back()[i]);
            if (words.back()[i] == 'i')
                type.push_back(EXPT);
            else if (i < words.back().length() - 1 && words.back()[i + 1] == '0')
            {
                type.push_back(M0);
                i++;
            }
            else if (i < words.back().length() - 1 && words.back()[i + 1] == '1')
            {
                type.push_back(M1);
                i++;
            }
            else
                type.push_back(EXPT);
            
        }
        
        std::vector<int> qubits;
        for (int i = 0; i < words.size() - 1; ++i)
        {
            std::string tmp;
            std::stringstream ss(words[i]);
            getline(ss, tmp, '[');
            getline(ss, tmp, ']');
            
            char *p;
            int index = strtol(tmp.c_str(), &p, 10);
            if (*p != 0 || tmp.length() == 0 || *p == tmp[0])
            {
                std::cerr << std::endl << "[warning]: In line \"" << raw_line << "\":" << std::endl;
                std::cerr << "Syntax \'" << tmp << "\' is not specifying a qubit. The line is ignored ..." << std::endl;
                return 1;
            }    
            
            qubits.push_back(index);
        }
        
        if (qubits.size() != basis.size())
        {
            std::cerr << std::endl << "[warning]: In line \"" << raw_line << "\":" << std::endl;
            std::cerr << "Number of qubits for calculating expection value does not match the length of Pauli string. This line is ignored..." << std::endl;
            return 1;
        }
        
        // backup
        int r_backup = r;
        int k_backup = k;
        DdNode **copy[w];
        for (int i = 0; i < w; i++)
        {
            copy[i] = new DdNode* [r];
            for (int j = 0; j < r; j++)
            {
                copy[i][j] = All_Bdd[i][j];
                Cudd_Ref(copy[i][j]);
            }
        }     
        DdNode *bigBDD_backup = bigBDD;
        Node_Table.clear();    // need to re-calculate
      
        // change basis
        for (int i = 0; i < qubits.size(); ++i)
        {        
            if (basis[i] == 'y') 
                Phase_shift_dagger(-2, qubits[i]);
            if (basis[i] == 'y' || basis[i] == 'x') 
                Hadamard(qubits[i]);
        } 
        
        // update bigBDD
        create_bigBDD();
        
        double oneroot2 = 1 / sqrt(2);
        double H_factor = pow(oneroot2, k%2);
        int nAnci_oneInt = ceil(log(r) / log(2)), nAnci_fourInt = ceil(log(w) / log(2)), nAnci = nAnci_oneInt + nAnci_fourInt, 
            nnAnci_fourInt = n + nAnci_fourInt, nVar = n + nAnci;
        
        int *permutation = new int[nVar];
        int indCount0 = 0;
        for (int i = 0; i < n; i++)
        {
            int index = Cudd_ReadInvPerm(manager, i);
            permutation[indCount0] = index;
            indCount0++;
        }
        for (int i = n; i < nVar; i++)
            permutation[i] = Cudd_ReadInvPerm(manager, i);
        int dum = Cudd_ShuffleHeap(manager, permutation);
        Node_Table.clear();    // need to re-calculate
        nodecount();
        
        // get value
        DdNode* expval_map = Cudd_ReadOne(manager);
        Cudd_Ref(expval_map);
        DdNode* collapse_map = Cudd_ReadOne(manager);
        Cudd_Ref(collapse_map);
        for (int i = 0; i < qubits.size(); ++i)
        {
            if (basis[i] == 'i')
                continue;
            
            if (type[i] == EXPT)
            {
                DdNode* tmp = Cudd_bddXnor(manager, expval_map, Cudd_Not(Cudd_bddIthVar(manager, qubits[i])));
                Cudd_Ref(tmp);
                Cudd_RecursiveDeref(manager, expval_map);
                expval_map = tmp;
            }
            else if (type[i] == M0)
            {
                DdNode* tmp = Cudd_bddAnd(manager, collapse_map, Cudd_Not(Cudd_bddIthVar(manager, qubits[i])));
                Cudd_Ref(tmp);
                Cudd_RecursiveDeref(manager, collapse_map);
                collapse_map = tmp;
            }
            else if (type[i] == M1)
            {
                DdNode* tmp = Cudd_bddAnd(manager, collapse_map, Cudd_bddIthVar(manager, qubits[i]));
                Cudd_Ref(tmp);
                Cudd_RecursiveDeref(manager, collapse_map);
                collapse_map = tmp;
            }
            else
                assert(false);
        }
        
        DdNode* expval_map_pos = Cudd_bddAnd(manager, expval_map, collapse_map);
        Cudd_Ref(expval_map_pos);
        DdNode* expval_map_neg = Cudd_bddAnd(manager, Cudd_Not(expval_map), collapse_map);
        Cudd_Ref(expval_map_neg);
        Cudd_RecursiveDeref(manager, expval_map);
        Cudd_RecursiveDeref(manager, collapse_map);
                
        double value = get_prob(expval_map_pos) - get_prob(expval_map_neg);
        Cudd_RecursiveDeref(manager, expval_map_pos);
        Cudd_RecursiveDeref(manager, expval_map_neg);
        
        std::string line = op_name;
        for (std::string& word : words)
            line = line + " " + word;
        std::ostringstream oss;
        oss << value;
        property.push_back(std::make_tuple(line, oss.str()));
        
        // restore
        for (int i = 0; i < w; i++)
        {
            for (int j = 0; j < r; j++)
                Cudd_RecursiveDeref(manager, All_Bdd[i][j]);
            delete[] All_Bdd[i];
        }
        r = r_backup;
        k = k_backup;
        for (int i = 0; i < w; i++)
        {
            All_Bdd[i] = new DdNode* [r];
            for (int j = 0; j < r; j++)
                All_Bdd[i][j] = copy[i][j];
            delete[] copy[i];
        }
        Cudd_RecursiveDeref(manager, bigBDD);
        bigBDD = bigBDD_backup;
        Node_Table.clear();    // need to re-calculate
        return 0;
    } 
    else if (keyword_truth.find(op_name) != keyword_truth.end())     // ============================ other supported
    {
        DdNode* function = get_property(op_name, words);
        double value = get_prob(function);
        Cudd_RecursiveDeref(manager, function);
        
        std::string line = op_name;
        for (std::string& word : words)
            line = line + " " + word;
        std::ostringstream oss;
        oss << value;
        property.push_back(std::make_tuple(line, oss.str()));
        return 0;
    }
    else if (keyword_condition.find(op_name) != keyword_condition.end())
    {
        if (condition_stack != "")
        {
            std::cerr << std::endl << "[warning]: In line \"" << raw_line << "\":" << std::endl;
            std::cerr << "Another \'" << condition_stack << "\' existing. The line is ignored ..." << std::endl;
            return 1;
        }
        
        for (int i = 0; i < 2; ++i)
        {
            char *p;       
            range[i] = strtof(words[i].c_str(), &p);
            if (*p != 0)
            {
                std::cerr << std::endl << "[warning]: In line \"" << raw_line << "\":" << std::endl;
                std::cerr << "Syntax \'" << words[i] << "\' is not supported for " << op_name << ". The line is ignored ..." << std::endl;
                return 1;
            }   
        }   
        
        condition_stack = op_name;
        return 0;
    }
    else                                                             // ============================ else
    {
        std::cerr << std::endl << "[warning]: In line \"" << raw_line << "\":" << std::endl;
        std::cerr << "Syntax \'" << op_name << "\' is not supported in this simulator. The line is ignored ..." << std::endl;
        return 1;
    }

}


/**Function*************************************************************

  Synopsis    [get the function of a self-defined property]

  Description []
               
  SideEffects [return nullptr if failed.
               returned node already referenced.]

  SeeAlso     []                  
                                 
***********************************************************************/
DdNode* Simulator::get_property(std::string& op_name, std::vector<std::string>& followings)
{
    if (op_name == KW_BOOL)                                // ============================ KW_BOOL
    {
        std::string function;
        for (std::string& word : followings) 
            function = function + word + " ";
        
        std::vector<std::string> func = boolean_parser(function);            
        return func2node(func);
    }
    if (op_name == KW_INTEQ || op_name == KW_INTNEQ || op_name == KW_INTGT || op_name == KW_INTLT ||
        op_name == KW_HWEQ  || op_name == KW_HWNEQ  || op_name == KW_HWGT  || op_name == KW_HWLT     )
    {
        DdNode* result;
        
        // special case
        int constant = stoi(followings.back());
        if (op_name == KW_INTEQ || op_name == KW_INTNEQ || op_name == KW_INTGT || op_name == KW_INTLT)
        {
            if (constant >= pow(2, followings.size() - 1))
            {
                result = Cudd_NotCond(Cudd_ReadOne(manager), (op_name == KW_INTEQ || op_name == KW_INTGT));
                Cudd_Ref(result);
                return result;
            }
        }
        if (op_name == KW_HWEQ  || op_name == KW_HWNEQ  || op_name == KW_HWGT  || op_name == KW_HWLT)
        {
            if (constant > followings.size() - 1)
            {
                result = Cudd_NotCond(Cudd_ReadOne(manager), (op_name == KW_HWEQ || op_name == KW_HWGT));
                Cudd_Ref(result);
                return result;
            }
        }
        if (constant < 0)
        {
            result = Cudd_NotCond(Cudd_ReadOne(manager), (op_name == KW_INTEQ || op_name == KW_INTLT || op_name == KW_HWEQ || op_name == KW_HWLT));
            Cudd_Ref(result);
            return result;
        }
        
        // parse
        std::vector<DdNode*> qubits;
        for (int i = 0; i < followings.size() - 1; ++i)
        {
            std::string tmp;
            std::stringstream ss(followings[i]);
            getline(ss, tmp, '[');
            getline(ss, tmp, ']');
            int ith_var = std::stoi(tmp);
            qubits.push_back(Cudd_bddIthVar(manager, ith_var));
            Cudd_Ref(qubits[i]);
        }
        
        // get node
        if (op_name == KW_INTEQ || op_name == KW_INTNEQ || op_name == KW_INTGT || op_name == KW_INTLT)
        {
            // constant nodes
            std::vector<DdNode*> constant_nodes(qubits.size(), nullptr);
            for (int i = qubits.size() - 1; i >= 0; --i)
            {
                constant_nodes[i] = Cudd_NotCond(Cudd_ReadOne(manager), (constant % 2 == 0));
                Cudd_Ref(constant_nodes[i]);
                constant = constant / 2;
            }
            
            // get result
            if (op_name == KW_INTEQ || op_name == KW_INTNEQ)
                result = Cudd_NotCond(node_equiv(qubits, constant_nodes), (op_name == KW_INTNEQ));
            else if (op_name == KW_INTGT)
                result = node_larger(qubits, constant_nodes);
            else // (op_name == KW_INTLT)
                result = node_larger(constant_nodes, qubits);
            
            // deRef
            for (int i = 0; i < qubits.size(); ++i)
            {
                Cudd_RecursiveDeref(manager, qubits[i]);
                Cudd_RecursiveDeref(manager, constant_nodes[i]);
            }
            return result;
        }
        
        else // (op_name == KW_HWEQ  || op_name == KW_HWNEQ  || op_name == KW_HWGT  || op_name == KW_HWLT)
        {
            // adder
            int n_adder = ceil(log(qubits.size() + 1) / log(2));
            std::vector<DdNode*> adder(n_adder, nullptr);
            for (int i = 0; i < n_adder; ++i)
            {
                adder[i] = Cudd_Not(Cudd_ReadOne(manager));
                Cudd_Ref(adder[i]);
            }
            for (int i = 0; i < qubits.size(); ++i)
            {
                DdNode* carry = qubits[i];
                Cudd_Ref(carry);
                for (int j = n_adder - 1; j >= 0; --j)
                {
                    DdNode* tmp1 = Cudd_bddAnd(manager, carry, adder[j]);
                    Cudd_Ref(tmp1);
                    DdNode* tmp2 = Cudd_bddXor(manager, carry, adder[j]);
                    Cudd_Ref(tmp2);
                    Cudd_RecursiveDeref(manager, carry);
                    Cudd_RecursiveDeref(manager, adder[j]);
                    carry = tmp1;
                    adder[j] = tmp2;
                }
                Cudd_RecursiveDeref(manager, carry);
            }
            
            // constant nodes
            std::vector<DdNode*> constant_nodes(n_adder, nullptr);
            for (int i = n_adder - 1; i >= 0; --i)
            {
                constant_nodes[i] = Cudd_NotCond(Cudd_ReadOne(manager), (constant % 2 == 0));
                Cudd_Ref(constant_nodes[i]);
                constant = constant / 2;
            }
        
            // get result
            if (op_name == KW_HWEQ || op_name == KW_HWNEQ)
                result = Cudd_NotCond(node_equiv(adder, constant_nodes), (op_name == KW_HWNEQ));
            else if (op_name == KW_HWGT)
                result = node_larger(adder, constant_nodes);
            else // (op_name == KW_HWLT)
                result = node_larger(constant_nodes, adder);
            
            // deRef
            for (int i = 0; i < qubits.size(); ++i)
                Cudd_RecursiveDeref(manager, qubits[i]);
            for (int i = 0; i < n_adder; ++i)
            {
                Cudd_RecursiveDeref(manager, adder[i]);
                Cudd_RecursiveDeref(manager, constant_nodes[i]);
            }
            return result;
        }
    }
    
    assert(false);
    return nullptr;
}
  
/**Function*************************************************************

  Synopsis    [get probability of a self-defined function]

  Description []
               
  SideEffects []

  SeeAlso     []

***********************************************************************/
double Simulator::get_prob(DdNode* function)
{
    DdNode *copy_bigBDD = bigBDD;
    Cudd_Ref(copy_bigBDD);
    
    DdNode* tmp = Cudd_bddAnd(manager, bigBDD, function);
    Cudd_Ref(tmp);
    Cudd_RecursiveDeref(manager, bigBDD);
    bigBDD = tmp;
    
    double oneroot2 = 1 / sqrt(2);
    double H_factor = pow(oneroot2, k%2);
    int nAnci_oneInt = ceil(log(r) / log(2)), nAnci_fourInt = ceil(log(w) / log(2)), nAnci = nAnci_oneInt + nAnci_fourInt, nnAnci_fourInt = n + nAnci_fourInt, nVar = n + nAnci;
    
    double value_0 = measure_probability(bigBDD, k/2, nVar, nAnci_fourInt, 0) * H_factor * H_factor;   
    double value_1 = measure_probability(bigBDD, k/2, nVar, nAnci_fourInt, 1) * H_factor * H_factor;
    
    Cudd_RecursiveDeref(manager, bigBDD);
    bigBDD = copy_bigBDD;
    return (value_0 + value_1);
}

/**Function*************************************************************

  Synopsis    [parse the result of property condition to string]

  Description []
               
  SideEffects []
 
  SeeAlso     []

***********************************************************************/

std::string Simulator::getCondStr(std::string& line)
{
    
    assert(condition_stack != "");
    if (condition_stack == KW_BTN)
    {
        std::ostringstream oss1, oss2;
        oss1 << range[0];
        oss2 << range[1];
        return oss1.str() + " <= " + line + " <= " + oss2.str();
    }
    if (condition_stack == KW_OOF)
    {
        std::ostringstream oss1, oss2;
        oss1 << range[0];
        oss2 << range[1];
        return "NOT " + oss1.str() + " <= " + line + " <= " + oss2.str();
    }
    if (condition_stack == KW_GEQ)
    {
        std::ostringstream oss1;
        oss1 << range[0];
        return line + " >= " + oss1.str();
    }
    if (condition_stack == KW_LEQ)
    {
        std::ostringstream oss1;
        oss1 << range[0];
        return line + " <= " + oss1.str();
    }
    assert(false);
    return "";
}

std::string Simulator::getCondRslt(double value)
{
    assert(condition_stack != "");
    if (condition_stack == KW_BTN)
        return (value >= range[0] && value <= range[1]) ? "true" : "false";
    if (condition_stack == KW_OOF)
        return !(value >= range[0] && value <= range[1]) ? "true" : "false";
    if (condition_stack == KW_GEQ)
        return (value >= range[0]) ? "true" : "false";
    if (condition_stack == KW_LEQ)
        return (value <= range[0]) ? "true" : "false";
    assert(false);
    return "";
}



/**Function*************************************************************

  Synopsis    [get probbility amplitude of a computational basis]

  Description []
               
  SideEffects []
 
  SeeAlso     []

***********************************************************************/
std::string Simulator::get_amplitude_string(int *assign)
{
    mpf_t re, im, sin_val, cos_val, tmp_float;
    mpf_init(re);
    mpf_init(im);
    mpf_init(sin_val);
    mpf_init(cos_val);
    mpf_init(tmp_float);    
    
    mpz_t int_value;
    mpz_init(int_value);    
    mpf_t int_value_as_float;
    mpf_init(int_value_as_float);
    
    mpf_t one_over_sqrt_2;        
    mpf_init(one_over_sqrt_2);
    mpf_sqrt_ui(one_over_sqrt_2, 2);
    mpf_div_ui(one_over_sqrt_2, one_over_sqrt_2, 2);
    // sqrt_val = 1/sqrt(2)
    
    mpz_t tmp_int;
    mpz_init(tmp_int);
    mpz_t two;
    mpz_init(two);
    mpz_set_str(two, "2", 10);
    mpf_t H_factor;
    mpf_init(H_factor);
    mpz_pow_ui(tmp_int, two, k / 2);
    mpf_set_z(H_factor, tmp_int);    
    if (k % 2 == 1)
        mpf_div(H_factor, H_factor, one_over_sqrt_2);
    
    DdNode *tmp;
    
    for (int j = 0; j < w; j++) // compute every complex value
    {
        std::string bitstring = "";
        
        mpz_init(int_value);
        for (int h = 0; h < r; h++) // compute every integer
        {
            tmp = Cudd_Eval(manager, All_Bdd[j][h], assign);
            Cudd_Ref(tmp);
            int oneEntry = !(Cudd_IsComplement(tmp));
            Cudd_RecursiveDeref(manager, tmp);
            if (oneEntry)
            {
                bitstring += "1";
            }
            else
            {
                bitstring += "0";
            }
            
        }
        std::reverse(bitstring.begin(), bitstring.end());
        bool isNeg = (bitstring[0]=='1');
        if (isNeg)
        {
            std::string inv_str = "1";
            for (int digitIndex = 0; digitIndex < r; digitIndex++)
            {
                inv_str += "0";
            }
            mpz_t raw_data;
            mpz_init_set_str(raw_data, bitstring.c_str(), 2);
            mpz_t inverter;
            mpz_init_set_str(inverter, inv_str.c_str(), 2);
            
            mpz_sub(int_value, raw_data, inverter);
        }
        else
        {
            mpz_set_str(int_value, bitstring.c_str(), 2);
        }
        
        switch (j)
        {
        case 3:
            mpf_set_d(cos_val, 1);
            mpf_set_d(sin_val, 0);
            break;
        case 2:
            mpf_set(cos_val, one_over_sqrt_2);
            mpf_set(sin_val, one_over_sqrt_2);
            break;
        case 1:
            mpf_set_d(cos_val, 0);
            mpf_set_d(sin_val, 1);
            break;
        case 0:
            mpf_neg(cos_val, one_over_sqrt_2); // -sqrt(2)/2
            mpf_set(sin_val, one_over_sqrt_2);
            break;
        default:
            std::cerr << "Warning: Unknown index (" << j << ")\n";
            break;
        }
                        
        mpf_set_z(int_value_as_float, int_value);
        
        mpf_mul(tmp_float, int_value_as_float, cos_val);
        mpf_add(re, re, tmp_float);
        mpf_mul(tmp_float, int_value_as_float, sin_val);
        mpf_add(im, im, tmp_float);
        // translate to re and im
    }
    
    mpf_div(tmp_float, re, H_factor);            
    long double final_re = mpf_get_d(tmp_float) * normalize_factor;
    mpf_div(tmp_float, im, H_factor);
    long double final_im = mpf_get_d(tmp_float) * normalize_factor;
    
    if ((final_re == 0)&&(final_im == 0))
        return "0";
    if (final_re == 0)
        return std::to_string(final_im) + "i";
    if (final_im == 0)
        return std::to_string(final_re);
    if (final_im < 0)
        return std::to_string(final_re) + std::to_string(final_im) + "i";
        
    return std::to_string(final_re) + "+" + std::to_string(final_im) + "i";
}

/**Function*************************************************************

  Synopsis    [get statevector string based on BDDs]

  Description []
               
  SideEffects []

  SeeAlso     []

***********************************************************************/
void Simulator::getStatevector()
{
    int *assign = new int[n];
    unsigned long long nEntries = pow(2, n);      // should timeout before overflow occurs in "nEntries"    
    for (int i = 0; i < n; i++)                   // initialize assignment
        assign[i] = 0;

    statevector = "[";
    for (unsigned long long i = 0; i < nEntries; i++) // compute every entry
    {
        long double final_re = 0;
        long double final_im = 0;
        bool isZero = 0;
        for (int j = 0; j < n; j++)
        {
            if (!measured_qubits_to_clbits[j].empty())
                if (assign[j] != stoi(measure_outcome.substr(n - 1 - j, 1)))
                {
                    isZero = 1;
                    break;
                }
        }
        
        std::string amplitude_string;
        if (isZero)
            amplitude_string = "0";
        else
            amplitude_string = get_amplitude_string(assign);
        
        statevector += amplitude_string;
        if (i != nEntries - 1)
            statevector = statevector + ", ";
        full_adder_plus_1(n, assign);
    }
    statevector += "]";

    delete[] assign;
}