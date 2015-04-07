#ifndef RANDOMBOOLEANNETWORK_HPP_
#define RANDOMBOOLEANNETWORK_HPP_

#include "ArrayDirectedGraph.hpp"
#include "Exception.hpp"
#include "FileFinder.hpp"
#include "OutputFileHandler.hpp"
#include <vector>
#include <map>
#include <bdd.h>
#include <iostream>
#include <fstream>

class RandomBooleanNetwork
{
private:

    /** The number of nodes in the network */
    unsigned mNodesNumber;

    /** The average of incoming edges for a single node */
    unsigned mAverageInputsPerNode; //K

    /** The graph of the boolean network */
    ArrayDirectedGraph* mpRbnGraph;

    /** The reverse transition function associated to the network (global) */
    bdd* mpReverseTransitionFunction;

    /**
     * Pointer to an array in which the i-th entry represent the function
     * associated with the i-th node.
     */
    bdd* mpNodeFunction;

    /**
     * Pointer to an array in which the i-th entry represent the function
     * associated with the i-th node at step n (= mpNodeFunction[i] ^ n).
     * Used for finding attractors.
     */
    bdd* mpNodeNthFunction;

    /**
     * Pointer to an array in which the i-th entry represent the
     * 2*i-th variable in the BDD. Remember that we have 2*mNodesNumber
     * variables. The evens represent the current state of the network
     * and the odds represent the next state.
     */
    bdd* mpVariables;

    /**
    * Pointer to an array in which the i-th entry represent the
    * 2*i-th+1 variable in the BDD. Remember that we have 2*mNodesNumber
    * variables. The evens represent the current state of the network
    * and the odds represent the next state.
    */
    bdd* mpNextVariables;

    /** Vector containing the attractors (as BDDs) */
    std::vector<bdd> mAttractors;

    /**
     * A vector of the same size of the mAttractors in which the
     * i-th entry represent the number of state of the i-th attractor
     */
    std::vector<unsigned> mAttractorLength;

    /**
     * Read a .cnet or .net file, build the graph and associate
     * a boolean function in every node.
     *
     * @param path the path of the input file
     */
    void createNetworkFromNetFile(const std::string path);

    /**
     * Read a .gml file and build the graph.
     *
     * @param path the path of the input file
     */
    void createGraphFromGmlFile(const std::string file_path);

    /**
     * Initialize buddy and setting some property.
     */
    void initBinaryDecisionDiagram() const;

    /**
     * Create a BDD boolean function and associate it to a given node
     *
     * @param node_id the ID of the node
     * @param canalyzing_function a bool indicating if the function must be
     * canalizing.
     */
    void createBooleanFunction(unsigned node_id, bool canalyzing_function);

    /**
     * Find the attractor id starting from a flipped state.
     *
     * @param flip_state bdd state
     * @param transition_function bdd transition function
     * @param set_variables bdd set of variables
     * @param replace_forward_assignment bddPair indicating the forward assignment
     */
    unsigned getStateAttractor(bdd flip_state, bdd transition_function, bdd set_variables, bddPair* replace_forward_assignment) const;

    /**
     * Given a set of states returning to they self, the method find the attractors
     * associated, and storage them in mAttractors.
     *
     * @param states_return_to_himself set of states returning to himself at step n
     * @param j step n
     * @param set_variables bdd set of variables
     * @param variables_pair bddPair indicating the assignment
     */
    void storageAttractors(bdd states_return_to_himself, unsigned j, bdd set_variables, bddPair* variables_pair);

    /**
     * This method compute delta_j+1, fixed delta_j. Transition function
     * T^j is memorized in mpNodeNthFunction[] and will be upgrated to T^j+1.
     *
     * @param variables_pair bddPair indicating the assignment
     *
     * @return bdd representing the states of delta_j+1.
     */
    bdd find_next_cycles(bddPair* variables_pair);

    /**
     * This method find the backward reachable states of a given bdd
     * 'states' (which can be a single or a set of states) using
     * the Reverse Transition Function.
     *
     * @param current_ring initial states
     * @param states set of states which is used to calculate F^-1(state)
     * iteratively
     * @param set_variables bdd set of variables
     * @param variables_pair bddPair indicating the assignment
     * @param steps_max number of steps max used
     */
    bdd find_backward_reachable_states(bdd current_ring, bdd states, bdd set_variables,
            bddPair* variables_pair, unsigned &steps_max);

    /**
     * Normalize a frequency matrix and obtain a probability distribution
     * in every row. This is the ATM as described here:
     * http://dx.plos.org/10.1371/journal.pone.0017703
     * http://dx.plos.org/10.1371/journal.pone.0017703.g001
     *
     * @return a stochastic matrix representing the ATM
     */
    std::vector<std::map<unsigned,double> > getStochasticMatrix(std::vector<std::map<unsigned,unsigned> >) const;

public:

    /**
     * Constructor 1: it generate the graph, and associate
     * a boolean function for every node of the network.
     *
     * @param nodes_number number of nodes in the network (N)
     * @param avarage_inputs_per_node average of incoming edges in the graph
     * @param scale_free topology parameter
     * @param probability_canalyzing_function probability to generate (random)
     * canalyzing functions for a node.
     */
    RandomBooleanNetwork(unsigned nodes_number, unsigned avarage_inputs_per_node,
            bool scale_free, double probability_canalyzing_function);

    /**
     * Constructor 2: it reads the graph from a file (.gml format), and
     * associate at random a boolean function for every node of the network.
     *
     * @param file_path path of the file.
     * @param probability_canalyzing_function probability to generate (random)
     * canalyzing functions for a node.
     */
    RandomBooleanNetwork(const std::string file_path, double probability_canalyzing_function);

    /**
     * Constructor 3: it reads the network (graph and functions) from a
     * file in .cnet or .net format.
     *
     * @param file_path path of the file.
     */
    RandomBooleanNetwork(const std::string file_path);

    /**
     * Destruptor: deallocate memory such as arrays and BDDs.
     */
    ~RandomBooleanNetwork();

    /**
     * Find attractors of the synchronous network using Zheng et al.
     * algorithm. http://dx.doi.org/10.1371/journal.pone.0060593
     *
     * It when it finish, mAttractors vector has a attractor for every
     * entry.
     */
    void findAttractors();

    /**
     * Method that induce noise in every bit of every state of the
     * attractors, and fill a frequency matrix in which rows and columns
     * are the attractors. Then, it normalises the frequency matrix.
     * http://dx.plos.org/10.1371/journal.pone.0017703
     *
     * @return a frequency matrix representing the ATN in frequency
     */
    std::vector<std::map<unsigned,double> > getAttractorMatrix()  const;

    /**
     * Save the network and if possible the attractors in a .net file.
     *
     * @param directory the name of the subfolder of testoutput.
     * @param filename the file name of the output file.
     */
    void printNetworkToNetFile(const std::string directory, const std::string filename) const;

    /**
     * Save the network and if possible the attractors in a .cnet file.
     *
     * @param directory the name of the subfolder of testoutput.
     * @param filename the name of the output file.
     */
    void printNetworkToCnetFile(const std::string directory, const std::string filename) const;

    /**
     * Save the network in a BoolNet (R package) file.
     *
     * @param directory the name of the subfolder of testoutput.
     * @param filename the name of the output file.
     */
    void printNetworkToBoolNetFile(const std::string directory, const std::string filename) const;

    /**
     * Save the network in a BooleanNet (py tool) file.
     *
     * @param directory the name of the subfolder of testoutput.
     * @param filename the file name of the output file.
     */
    void printNetworkToBooleanNetFile(const std::string directory, const std::string filename) const;

    /**
     * Save the graph (functions are lost) in a .gml file.
     *
     * @param directory the name of the subfolder of testoutput.
     * @param filename the name of the output file.
     */
    void printGraphToGmlFile(const std::string directory, const std::string filename) const;

    /**
     * Save the graph (functions are lost) in a .sif file.
     *
     * @param directory the name of the subfolder of testoutput.
     * @param filename the file name of the output file.
     */
    void printGraphToSifFile(const std::string directory, const std::string filename) const;

    /**
     * Save the graph (functions are lost) in a .dot file.
     *
     * @param directory the name of the subfolder of testoutput.
     * @param filename the file name of the output file.
     */
    void printGraphToDotFile(const std::string directory, const std::string filename) const;

    /**
     * print the graph in the console.
     */
    void printNetwork() const;

    /**
     * getter of mNodesNumber.
     *
     * @return mNodesNumber
     */
    unsigned getNodesNumber() const;

    /**
     * getter of mAverageInputsPerNode.
     *
     * @return mAverageInputsPerNode
     */
    unsigned getAvarageInputsPerNode() const;

    /**
     * @return the number of attractors founded
     */
    unsigned getAttractorsNumber() const;

    /**
     * getter of mAttractorLength.
     *
     * @return mAttractorLength.
     */
    std::vector<unsigned> getAttractorLength() const;
};

#endif /* RANDOMBOOLEANNETWORK_HPP_ */
