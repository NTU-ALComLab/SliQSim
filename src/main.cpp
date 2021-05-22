#include <boost/program_options.hpp>
#include "Simulator.h"
#include "util_sim.h"

int main(int argc, char **argv)
{
    namespace po = boost::program_options;
    po::options_description description("Options");
    description.add_options()
    ("help", "produce help message")
    ("sim_qasm", po::value<std::string>()->implicit_value(""), "simulate qasm file string")
    ("seed", po::value<unsigned int>()->implicit_value(1), "seed for random number generator")
    ("print_info", "print simulation statistics such as runtime, memory, etc.")
    ("type", po::value<unsigned int>()->default_value(0), "the simulation type being executed.\n" 
                                                           "0: weak simulation (default option) where the sampled outcome(s) will be provided after the simulation. " 
                                                           "The number of outcomes being sampled can be set by argument \"shots\" (1 by default).\n"
                                                           "1: strong simulation where the resulting state vector will be shown after the simulation. "
                                                           "Note that this option should be avoided if the number of qubits is large since it will result in extremely long runtime.")
    ("shots", po::value<unsigned int>()->default_value(1), "the number of outcomes being sampled, " 
                                                          "this argument is only used when the " 
                                                          "simulation type is set to \"weak\"")
    ;
    
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, description), vm);
    po::notify(vm);

    if (vm.count("help") || argc == 1) 
    {
	    std::cout << description << std::endl;
	    return 1;
	}

    struct timeval t1, t2;
    double elapsedTime;

    int type = vm["type"].as<unsigned int>(), shots = vm["shots"].as<unsigned int>();
    if (type == 1)
    {
        shots = 1;
    }
    else
    {
        type = 0;
    }
    std::random_device rd;
    unsigned int seed;
    if (vm.count("seed")) 
		seed = vm["seed"].as<unsigned int>();
    else
        seed = rd();

    // start timer
    gettimeofday(&t1, NULL);

    assert(shots > 0);
    Simulator simulator(type, shots, seed);

    if (vm.count("sim_qasm"))
    {
        // read in file into a string
        std::stringstream strStream;
        if (vm["sim_qasm"].as<std::string>() == "")
        {
            strStream << std::cin.rdbuf();
        }
        else
        {
            std::ifstream inFile;
            inFile.open(vm["sim_qasm"].as<std::string>()); //open the input file
            strStream << inFile.rdbuf(); //read the file
        }
        
        std::string inFile_str = strStream.str(); //str holds the content of the file
        simulator.sim_qasm(inFile_str);
    }

    //end timer
    gettimeofday(&t2, NULL);
    elapsedTime = (t2.tv_sec - t1.tv_sec) * 1000.0;
    elapsedTime += (t2.tv_usec - t1.tv_usec) / 1000.0;

    double runtime = elapsedTime / 1000;
    size_t memPeak = getPeakRSS();
    if (vm.count("print_info"))
        simulator.print_info(runtime, memPeak);

    return 0;
}