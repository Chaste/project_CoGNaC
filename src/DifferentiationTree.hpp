#ifndef DIFFERENTIATIONTREE_HPP_
#define DIFFERENTIATIONTREE_HPP_

#include "DifferentiationTreeNode.hpp"
#include "Exception.hpp"
#include "OutputFileHandler.hpp"
#include <list>
#include <map>
#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>

class DifferentiationTree
{
private:
	/** A Boolean value indicating if the tree has a fake root */
	bool mFakeRoot;

    /** Vertices of the tree */
    std::vector<DifferentiationTreeNode*> mVertices;

    /** Set of Leaves*/
    std::set<unsigned> mLeaves;

    /**
     * A vector of set, in which the i-th entry represent the set
     * of vertices at level i-th.
     */
    std::vector<std::set<unsigned> > mLevelNodes;

    /** A vector in which for every node it is given the level */
    std::vector<unsigned> mNodesLevel;

    /**
     * Recursive method for check the isomorphism between trees.
     * Taken from "The design and analysis of computer algorithms",
     * Aho, Hopcroft, Ullman, pp. 84 - 85 - 86.
     *
     * Levels starting zero (root). Note that this order is inverse
     * respect the Aho, Hopcroft, Ullman algorithm.
     *
     * @param step: a decreasing unsigned integer. Algorithm stops
     *                  if (step == 0) or some of the recursive call
     *                  return false. It can be viewed as the current
     *                  level of the trees, starting from the lower level
     *                  (leaves) to the zero level (root);
     * @param tree1: a pointer to tree1 (const);
     * @param tree2: a pointer to tree2 (const);
     * @param levelNodes1: a reference to the vector levelNodes of the
     *                  tree1. We take it as reference for efficiency (with const keyword).
     *                  In the i-th level it give a set of node ids
     *                  belong to the i-th level (root belong to level zero, not one)
     * @param levelNodes2: a reference to the vector levelNodes of the
     *                  tree2.
     * @param sub_tree_type1: data structure representing sub tree
     * @param sub_tree_type2: data structure representing sub tree
     * @param tuples2: vector of tuples
     * @param tuples2: vector of tuples
     */
    bool checkIsomorphism(
            unsigned step,
            const DifferentiationTree* tree1,
            const DifferentiationTree* tree2,
            const std::vector<std::set<unsigned> >& levelNodes1,
            const std::vector<std::set<unsigned> >& levelNodes2,
            std::map<unsigned, unsigned>& sub_tree_type1,
            std::map<unsigned, unsigned>& sub_tree_type2,
            std::vector<std::list<unsigned> >& tuples1,
            std::vector<std::list<unsigned> >& tuples2
            ) const;
public:
    /**
     * Constructor 1: empty
     */
    DifferentiationTree();
    /**
     * Constructor 2:
     *
     * @param cell cycle length: the duration (in hours) of the root's cell cycle .
     */
    DifferentiationTree(double cell_cycle_length);

    /**
     * Constructor 3:
     *
     * @param component: a set containing the terminal strongly connected components
     * of the set.
     */
    DifferentiationTree(std::set<unsigned> component);

    /**
     * Constructor 4:
     *
     * @param cell_cycle_length: the duration (in hours) of the root's cell cycle.
     * @param component_states: a set containing the scc of the TES.
     * @param stationary_distribution: a probability vector where the i-th entry represent the i-th scc.
     */
    DifferentiationTree(
            double cell_cycle_length,
            std::set<unsigned> component_states,
            std::vector<double> stationary_distribution,
	        std::vector<double> differentiation_probability
            );
    /**
     * Destruptor: deallocate all the memory deleting every element
     * in the vector (mVertices).
     */

    ~DifferentiationTree();

    /**
     * Recursive method used for find in the tree the node which have
     * a the component states set which contains every element
     * of child_component_states set, and it is the smaller
     * set with this property in the tree.
     *
     * @param starting_node (starting from the root)
     * @param child_component_states
     */
    unsigned searchParentNode(std::set<unsigned> child_component_states,
            unsigned starting_node = 0) const;

    /**
     * If the tree has a fake root, it returns a vector of differentiation trees
     * where the roots of these trees are its children. These trees have not a fake
     * root, so can be used during simulations.
     *
     * @return a vector of pointers to DifferentiationTree objects which
     * should be deleted by the user at the end.
     */
    std::vector<DifferentiationTree*> SplitTreeFromFakeRoot() const;

    /**
     * This method normalize the CellCycleLength attribute of every
     * node in the tree based on a user-defined hours mean.
     * It is the second and last step after
     * ThresholdEgodicSetDifferentiationTree::getCellCycleLength.
     *
     * http://dx.plos.org/10.1371/journal.pone.0097272
     *
     * @param hours_mean
     */
    void normaliseLength(double hours_mean);

    /**
     * This method change the mColour attribute of every node in the tree.
     * It is used for Paraview visualization, colouring with different colours
     * stem cells (root), transit cells and fully differentiated cells
     * (leaves). Colour values are in the range [0, 1].
     */
    void colorise();

    /**
     * Set colour to a node in the tree
     *
     * @param node_id
     * @param colour
     */
    void setColour(unsigned node_id, double colour);

    /**
     * Add a child to a given parent id
     */
    void addChild(unsigned parent, DifferentiationTreeNode* node);

    /**
     * Add a new empty child to a given parent id.
     *
     * @param parent
     */
    void addNewChild(unsigned parent);

    /**
     * Add a new child setting its cell cycle length
     * to a given parent id.
     *
     * @param parent
     * @param cell_cycle_length
     */
    void addNewChild(unsigned parent, double cell_cycle_length);

    /**
     * Add a new child setting its component
     * to a given parent id.
     *
     * @param parent
     * @param component
     */
    void addNewChild(unsigned parent, std::set<unsigned> component);

    /**
     * Add a new child setting its cell cycle length,
     * component states and probability distribution
     * to a given parent id.
     *
     * @param parent
     * @param cell_cycle_length
     * @param component_states
     * @param stationary_distribution
     * @param differentiation_probability
     */
    void addNewChild(
            unsigned parent,
            double cell_cycle_length,
            std::set<unsigned> component_states,
            std::vector<double> stationary_distribution,
            std::vector<double> differentiation_probability
            );

    /**
     * Add node at the tree. The correct parent is found checking the
     * components states of the node.
     *
     * @param node
     */
    void addNode(DifferentiationTreeNode* node);

    /**
     * Print the tree in the console output.
     */
    void printDifferentiationTree() const;

    /**
     * Save the differentiation tree in a .gml file.
     *
     * @param directory the name of the subfolder of testoutput.
     * @param filename the file name of the output file.
     */
    void printDifferentiationTreeToGmlFile(std::string directory, std::string filename) const;

    /**
     * Check if this tree and the user_defined_tree are isomorphic.
     *
     * @param user_defined_tree
     *
     * @return true if the trees are isomorphic
     */
    bool topologyTreeCompare(const DifferentiationTree* user_defined_tree) const;

    /**
     * @return the root of the tree
     */
    const DifferentiationTreeNode* getRoot() const;

    /**
     * Get the node_id-th node of the tree.
     *
     * @param node_id
     *
     * @return the node_id-th node if it exists.
     */
    const DifferentiationTreeNode* getNode(unsigned node_id) const;

    /**
     * @return the number of vertices of the tree.
     */
    unsigned size () const;

    /**
     * @return true the tree has a fake root.
     */
    bool hasFakeRoot() const;

    /**
     * Getter of mLeaves
     *
     * @return mLevaves
     */
    std::set<unsigned> getLeaves () const;

    /**
     * Getter of mLevelNodes
     *
     * @return mLevelNodes
     */
    std::vector<std::set<unsigned> > getLevelNodes () const;

    /**
     * Getter of mNodesLevel
     *
     * @return mNodesLevel
     */
    std::vector<unsigned> getNodesLevel() const;

    /**
     * Setter of mFakeRoot
     *
     * @param fake_root
     */
    void setFakeRoot(bool fake_root);

    /**
     * Set a threshold value to a node
     *
     * @param node_id
     * @param threshold
     */
    void setThreshold(unsigned node_id, double threshold);

    /**
     * Set a stationary distribution and cell cycle length to a node
     *
     * @param node_id
     * @param stationary_distribution
     * @param cell_cycle_length
     */
    void setStationaryDistributionAndCellCycleLength(
    		unsigned node_id,
			std::vector<double> stationary_distribution,
			double cell_cycle_length);
    /**
	 * Set a differentiation probability distribution to a node
	 *
	 * @param differentiation_probability
	 */
    void setDifferentiationProbability(
        		unsigned node_id,
    			std::vector<double> differentiation_probability);


};

#endif /* DIFFERENTIATIONTREE_HPP_ */
