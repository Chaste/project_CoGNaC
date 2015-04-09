#ifndef THRESHOLDERGODICSETDIFFERENTIATIONTREE_HPP_
#define THRESHOLDERGODICSETDIFFERENTIATIONTREE_HPP_

#include <vector>
#include <map>
#include <set>
#include "DifferentiationTree.hpp"
#include "RandomBooleanNetwork.hpp"
#include "Exception.hpp"
#include "LinearSystem.hpp"
#include "PetscTools.hpp"
#include "PetscVecTools.hpp"
#include "PetscMatTools.hpp"


class ThresholdErgodicSetDifferentiationTree
{
private:
    /** Stochastic matrix representing the ATN */
    std::vector<std::map<unsigned,double> > mStochasticMatrix;

    /**
     * A vector of the same size of the mAttractors in which the
     * i-th entry represent the number of state of the i-th attractor
     */
    std::vector<unsigned> mAttractorLength;

    /** A pointer to a Boolean Network */
    RandomBooleanNetwork* mpBooleanNetwork;

    /** A pointer to a DifferentiationTree */
    //DifferentiationTree* mpDifferentiationTree;

    /**
     * Method to compute a differentiation tree starting from a
     * stochastic matrix and the attractor lengths.
     */
    DifferentiationTree* computeDifferentiationTree() const;

    /**
     * Return all the different values >= 0 conteined in the stochastic
     * matrix.
     *
     * @return a set of double, all the possible threshold candidates.
     */
    std::set<double> getThresholdValues() const;

    /**
     * Normalize a pruned matrix in an ergodic matrix
     *
     * @param matrix pruned matrix
     */
    void normalizePrunedMatrix(std::vector<std::map<unsigned,double> > &matrix) const;

    /**
     * Get a pruned matrix from the stochastic matrix
     *
     * @param threshold value
     *
     * @return pruned matrix
     */
    std::vector<std::map<unsigned,double> > getPrunedMatrix(double threshold) const;

    /**
     * Get the strongly connected components from a matrix (see as
     * Adjacency table of a graph where exist an edge if the value
     * matrix[i][j] is > 0), using the Tarjan's algorithm.
     *
     * @param matrix
     *
     * @return set of sets ssc of the graph.
     */
    std::set<std::set<unsigned> > getStronglyConnectedComponents(std::vector<std::map<unsigned,double> > graph) const;

    /**
     * Get the strongly connected components from a matrix (see as
     * Adjacency table of a graph where exist an edge if the value
     * matrix[i][j] is > 0), using the Tarjan's algorithm.
     *
     * @param matrix
     *
     * @return set of sets ssc of the graph.
     */
    std::set<std::set<unsigned> > getTerminalStronglyConnectedComponents(std::vector<std::map<unsigned,double> > graph) const;

    /**
     * Helper of the getStronglyConnectedComponents() method. Find the
     * components using depth first search.
     *
     * @param u
     * @param graph
     * @param lowlink
     * @param used
     * @param stack
     * @param time
     * @param components
     */
    void findComponentDepthFirstSearch(unsigned u, std::vector<std::map<unsigned,double> > graph, unsigned* lowlink,
            bool* used, std::vector<unsigned> &stack, unsigned &time, std::set<std::set<unsigned> > &components) const;

    /**
     * Find the stationary distribution of the TES using the
     * exact method.
     *
     * @param tes_map_matrix
     * @param component
     *
     * @return a stochastic vector containing the stationary distribution
     */
    std::vector<double> findStationaryDistribution(std::vector<std::map<unsigned, double> > tes_map_matrix, std::set<unsigned> component) const;

    /**
     * The length of the cell cycle of the TES is obtained by weighted sum
     * http://dx.plos.org/10.1371/journal.pone.0097272
     * It will be normalized in the DifferentiationTreeBasedCellCycleModel
     * on the based of a user defined time scale link.
     *
     * @param stationary_distribution stationary probability vector
     * @param component
     *
     * @return length of the cell cycle model
     */
    double getCellCycleLength(std::vector<double> stationary_distribution, std::set<unsigned> component) const;

    /**
     * This method assigns probabilities and a cell cycle length
     * to each node of the tree.
     *
     * @param a differentiation tree.
     */
    void assignProbabilitiesAndCellCycleLengths(DifferentiationTree* differentiation_tree) const;

public:

    /**
     * Constructor 1: compute the differentiation tree from a stochastic
     * matrix and a vector of attractors lengths.
     *
     * @param m_stochastic_matrix
     * @param m_attractor_length
     */
    ThresholdErgodicSetDifferentiationTree(std::vector<std::map<unsigned,double> > m_stochastic_matrix,
            std::vector<unsigned> m_attractor_length);

    /**
     * Constructor 2: create a RandomBooleanNetwork object and compute
     * the differentiation tree from the boolean network. This constructor
     * call the constructor 1 of the RandomBooleanNetwork class.
     *
     * @param nodes_number number of nodes in the network (N)
     * @param avarage_inputs_per_node average of incoming edges in the graph
     * @param scale_free topology parameter
     * @param probability_canalyzing_function probability to generate (random)
     * canalyzing functions for a node.
     */
    ThresholdErgodicSetDifferentiationTree(unsigned nodes_number, unsigned avarage_input_number_per_node,
            bool scale_free, double probability_canalyzing_function);

    /**
     * Constructor 3: create a RandomBooleanNetwork object and compute
     * the differentiation tree from the boolean network. This constructor
     * call the constructor 2 of the RandomBooleanNetwork class.
     *
     * @param file_path path of the file.
     * @param probability_canalyzing_function probability to generate (random)
     * canalyzing functions for a node.
     */
    ThresholdErgodicSetDifferentiationTree(const std::string file_path, double probability_canalyzing_function);

    /**
     * Constructor 4: create a RandomBooleanNetwork object and compute
     * the differentiation tree from the boolean network. This constructor
     * call the constructor 3 of the RandomBooleanNetwork class.
     *
     * @param file_path path of the file.
     *
     */
    ThresholdErgodicSetDifferentiationTree(const std::string file_path);

    /**
     * Distruptor: delete the RandomBooleanNetwork object and
     * the DifferentiationTree object if the are not NULL.
     */
    ~ThresholdErgodicSetDifferentiationTree();

    /**
     * getter of mpRandomBooleanNetwork.
     *
     * @return mpRandomBooleanNetwork.
     */
    const RandomBooleanNetwork* getBooleanNetwork() const;

    /**
     * getter of mAttractorLength.
     *
     * @return mAttractorLength.
     */
    std::vector<unsigned> getAttractorLength() const;

    /**
     * getter of mStochasticMatrix.
     *
     * @return mStochasticMatrix.
     */
    std::vector<std::map<unsigned,double> > getStochasticMatrix() const;

    /**
     * Return a pointer to a DifferentiationTree. This tree is
     * created using the stochastic matrix and the attractor length
     * as threshold-dependent ATN.
     * GestoDifferent plugin for Cytoscape do the same.
     * doi: 10.1093/bioinformatics/bts726
     *
     * In addition, we invoke colorise() method on the tree.
     *
     * We don't impose const to this object because users
     * can set their own time duration, or colorise the
     * tree.
     *
     * User should delete the object explicitly after use.
     *
     * @return a pointer to a new DifferentiationTree,.
     */
    DifferentiationTree* getDifferentiationTree() const;

    /**
	 * Test if the network is an ergodic set (or composed only
	 * by ergodic sets).
	 *
	 * @param matrix
	 * @param components: a set containing the terminal components.
     *
     * @return bool true the matrix is a TES
     */
    bool isThresholdErgodicSet(std::vector<std::map<unsigned,double> > matrix, const std::set<std::set<unsigned> > components) const;

    /**
	 * Save the stochastic matrix and the lengths of the attractors in a .dat file.
	 *
	 * @param directory the name of the subfolder of testoutput.
	 * @param filename the file name of the output file.
     */
    void printStochasticMatrixAndAttractorLengthsToDatFile(std::string directory, std::string filename) const;

};

#endif /* THRESHOLDERGODICSETDIFFERENTIATIONTREE_HPP_ */
