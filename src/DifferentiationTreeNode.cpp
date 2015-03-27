#include "DifferentiationTreeNode.hpp"
#include <iostream>
#include <cassert>

DifferentiationTreeNode::DifferentiationTreeNode() :
mColour (0.0),
mCellCycleLength(0.0),
mParent(0)
{
}

DifferentiationTreeNode::DifferentiationTreeNode(double cell_cycle_length) :
mColour (0.0),
mCellCycleLength(cell_cycle_length),
mParent(0)
{
    if (mCellCycleLength < 0.0) EXCEPTION("Cell cycle length must be non-negative.");
}

DifferentiationTreeNode::DifferentiationTreeNode(double cell_cycle_length,
            std::set<unsigned> component_states, std::vector<double> stationary_distribution) :
mColour (0.0),
mCellCycleLength(cell_cycle_length),
mParent(0),
mComponentStates(component_states),
mStationaryDistribution(stationary_distribution)
{
    if (mCellCycleLength < 0.0) EXCEPTION("Cell cycle length must be non-negative.");
    if (mStationaryDistribution.size() != mComponentStates.size())
        EXCEPTION("Number of components must be equal to the size of the vector of probabilities.");
}

DifferentiationTreeNode::~DifferentiationTreeNode()
{
}

void DifferentiationTreeNode::addChild(unsigned child)
{
    mChildren.push_back(child);
}

void DifferentiationTreeNode::setParent(unsigned parent)
{
    mParent = parent;
}

void DifferentiationTreeNode::removeLastChild()
{
    if (mChildren.size() > 0)
    {
        mChildren.pop_back();
    }
}

bool DifferentiationTreeNode::isLeaf() const
{
    return (mChildren.size() == 0);
}
unsigned DifferentiationTreeNode::getNumberOfComponentStates() const
{
    return mComponentStates.size();
}

unsigned DifferentiationTreeNode::getNumberOfChildren() const
{
    return mChildren.size();
}

double DifferentiationTreeNode::getColour() const
{
    return mColour;
}

double DifferentiationTreeNode::getCellCycleLength() const
{
    return mCellCycleLength;
}
unsigned DifferentiationTreeNode::getParent() const
{
    return mParent;
}
std::set<unsigned> DifferentiationTreeNode::getComponentStates() const
{
    return mComponentStates;
}
std::vector<unsigned> DifferentiationTreeNode::getChildren() const
{
    return mChildren;
}

std::vector<double> DifferentiationTreeNode::getStationaryDistribution() const
{
    return mStationaryDistribution;
}

void DifferentiationTreeNode::setStationaryDistribution(std::vector<double> distribution)
{
    mStationaryDistribution = distribution;
}

void DifferentiationTreeNode::setColour(double colour)
{
    mColour = colour;
}

void DifferentiationTreeNode::setCellCycleLength(double length)
{
    mCellCycleLength = length;
}
