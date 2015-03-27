#ifndef DIFFERENTIATIONTREENODE_HPP_
#define DIFFERENTIATIONTREENODE_HPP_

#include <set>
#include <vector>
#include <iostream>
#include "Exception.hpp"

class DifferentiationTreeNode
{
private:
    /** Colour used for paraview simulations. It is a value in [0,1]*/
    double mColour;

    /** Cell cycle length in hours.*/
    double mCellCycleLength;

    /** Id of the parent */
    unsigned mParent;

    /** Components of the node (TES) */
    std::set<unsigned> mComponentStates;

    /** Children of the nodes */
    std::vector<unsigned> mChildren;

    /** Probability distribution. It must be the same size of mChildren*/
    std::vector<double> mStationaryDistribution;

public:
    /**
     * Constructor 1: empty node. Used only for topology tree compare.
     */
    DifferentiationTreeNode();

    /**
     * Constructor 2: only cell_cycle_length as parameter. Used for
     * used-defined tree for DifferentiationTreeBasedCellCycleModel
     * simulations.
     *
     * @param cell_cycle_length
     */
    DifferentiationTreeNode(double cell_cycle_length);

    /**
     * Constructor 3. Used by ThresholdErgodicSetDifferentiationTree
     *
     * @param cell_cycle_length
     * @param component_states
     * @param stationary_distribution
     */
    DifferentiationTreeNode(double cell_cycle_length, std::set<unsigned> component_states,
            std::vector<double> stationary_distribution);

    /** Distructor */
    ~DifferentiationTreeNode();

    /**
     * Add a child to this node.
     *
     * @param child
     */
    void addChild(unsigned child);

    /**
     * Remove the last child of this node.
     */
    void removeLastChild();

    /**
     * Check if this node is a leaf.
     *
     * @return true if it has not children
     */
    bool isLeaf() const;

    /**
     * @return cardinality of mComponentStates
     */
    unsigned getNumberOfComponentStates() const;

    /**
     * @return number of children
     */
    unsigned getNumberOfChildren() const;

    /**
     * Getter of mColour
     *
     * @return mColour
     */
    double getColour() const;

    /**
     * Getter of mCellCycleLength
     *
     * @return mCellCycleLength
     */
    double getCellCycleLength() const;

    /**
     * Getter of mParent
     *
     * @return mParent
     */
    unsigned getParent() const;

    /**
     * Getter of mComponentStates
     *
     * @return mComponentStates
     */
    std::set<unsigned> getComponentStates() const;

    /**
     * Getter of mChildren
     *
     * @return mChildren
     */
    std::vector<unsigned> getChildren() const;

    /**
     * Getter of mStationaryDistribution
     *
     * @return mStationaryDistribution
     */
    std::vector<double> getStationaryDistribution() const;

    /**
     * Setter of mStationaryDistribution
     *
     * @param distribution
     */
    void setStationaryDistribution(std::vector<double> distribution);

    /**
     * Setter of mColour
     *
     * @param colour
     */
    void setColour(double colour);

    /**
     * Setter of mCellCycleLength
     *
     * @param mCellCycleLength
     */
    void setCellCycleLength(double length);

    /**
     * Setter of mParent
     *
     * @param parent
     */
    void setParent(unsigned parent);
};

#endif /* DIFFERENTIATIONTREENODE_HPP_ */
