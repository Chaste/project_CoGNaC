#ifndef TESTDIFFERENTIATIONTREE_HPP_
#define TESTDIFFERENTIATIONTREE_HPP_

#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "SmartPointers.hpp"

#include "DifferentiationTree.hpp"
#include "DifferentiationTreeNode.hpp"
#include <set>
#include <vector>

class TestDifferentiationTree : public CxxTest::TestSuite
{
public:
    void testDifferentiationTree()
    {
        std::set<unsigned> set0;
        set0.insert(0);
        set0.insert(1);
        std::set<unsigned> set1;
        set1.insert(0);
        std::set<unsigned> set2;
        set2.insert(1);
        std::vector<double> dist0;
        dist0.push_back(.5);
        dist0.push_back(.5);
        std::vector<double> dist1;
        dist1.push_back(1);

        DifferentiationTree* tree1 = new DifferentiationTree(8, set0, dist0);
        TS_ASSERT(tree1->getLeaves().size() == 1);

        unsigned parent = tree1->searchParentNode(set1);
        tree1->addNewChild(parent, 4, set1, dist1);
        TS_ASSERT(tree1->getLeaves().size() == 1);

        parent = tree1->searchParentNode(set2);
        tree1->addNewChild(parent, 3, set2, dist1);
        TS_ASSERT(tree1->getLeaves().size() == 2);

        std::vector<unsigned> nodes_level = tree1->getNodesLevel();
        for (unsigned i=0;i<tree1->size(); i++)
        {
            if (i==0)
            {
                TS_ASSERT(nodes_level.at(i)  == 0);
            }
            else
            {
                TS_ASSERT(nodes_level.at(i)  == 1);
            }
        }
        std::vector<std::set<unsigned> > level_nodes = tree1->getLevelNodes();
        for (unsigned i=0; i<level_nodes.size(); i++)
        {
            if (i==0)
            {
                TS_ASSERT(level_nodes.at(i).size()  == 1);
            }
            else
            {
                TS_ASSERT(level_nodes.at(i).size()  == 2);
            }
        }
        const DifferentiationTreeNode* root = tree1->getRoot();
        TS_ASSERT(root->getNumberOfChildren() == 2);
        tree1->normaliseLength(1.5);

        TS_ASSERT_DELTA(tree1->getNode(0)->getCellCycleLength(), 2.4, 1e-4);
        TS_ASSERT_DELTA(tree1->getNode(1)->getCellCycleLength(), 1.2, 1e-4);
        TS_ASSERT_DELTA(tree1->getNode(2)->getCellCycleLength(), 0.9, 1e-4);

        tree1->colorise();
        TS_ASSERT_DELTA(tree1->getNode(0)->getColour(), 0.0, 1e-4);
        TS_ASSERT_DELTA(tree1->getNode(1)->getColour(), 1.0, 1e-4);
        TS_ASSERT_DELTA(tree1->getNode(2)->getColour(), 0.9, 1e-4);

        delete tree1;
    }
    void testComparingTrees1()
    {
        DifferentiationTree* tree1 = new DifferentiationTree();
        tree1->addNewChild(0);
        tree1->addNewChild(0);
        tree1->addNewChild(2);
        tree1->addNewChild(2);
        tree1->addNewChild(3);
        tree1->addNewChild(3);
        tree1->addNewChild(4);

        DifferentiationTree* tree2 = new DifferentiationTree();
        tree2->addNewChild(0);
        tree2->addNewChild(0);
        tree2->addNewChild(2);
        tree2->addNewChild(2);
        tree2->addNewChild(2);
        tree2->addNewChild(3);
        tree2->addNewChild(4);

        TS_ASSERT_EQUALS(tree1->topologyTreeCompare(tree2),tree2->topologyTreeCompare(tree1));
        TS_ASSERT_EQUALS(tree1->topologyTreeCompare(tree2),false);
        delete tree1;
        delete tree2;

    }

    void testComparingTrees2()
    {
        DifferentiationTree* tree1 = new DifferentiationTree();
        tree1->addNewChild(0);
        tree1->addNewChild(0);
        tree1->addNewChild(1);
        tree1->addNewChild(1);
        tree1->addNewChild(1);
        tree1->addNewChild(2);
        tree1->addNewChild(2);
        tree1->addNewChild(2);
        tree1->addNewChild(6);
        tree1->addNewChild(6);
        tree1->addNewChild(8);
        tree1->addNewChild(8);

        DifferentiationTree* tree2 = new DifferentiationTree();
        tree2->addNewChild(0);
        tree2->addNewChild(0);
        tree2->addNewChild(1);
        tree2->addNewChild(1);
        tree2->addNewChild(1);
        tree2->addNewChild(3);
        tree2->addNewChild(3);
        tree2->addNewChild(4);
        tree2->addNewChild(4);
        tree2->addNewChild(2);
        tree2->addNewChild(2);
        tree2->addNewChild(2);

        TS_ASSERT_EQUALS(tree1->topologyTreeCompare(tree2),tree2->topologyTreeCompare(tree1));
        TS_ASSERT_EQUALS(tree1->topologyTreeCompare(tree2),true);
        delete tree1;
        delete tree2;
    }

    void testComparingTrees3()
    {
        DifferentiationTree* tree1 = new DifferentiationTree();
        tree1->addNewChild(0);
        tree1->addNewChild(0);
        tree1->addNewChild(1);
        tree1->addNewChild(1);
        tree1->addNewChild(1);
        tree1->addNewChild(2);
        tree1->addNewChild(2);
        tree1->addNewChild(2);
        tree1->addNewChild(6);
        tree1->addNewChild(6);
        tree1->addNewChild(8);
        tree1->addNewChild(8);
        tree1->addNewChild(2);

        DifferentiationTree* tree2 = new DifferentiationTree();
        tree2->addNewChild(0);
        tree2->addNewChild(0);
        tree2->addNewChild(1);
        tree2->addNewChild(1);
        tree2->addNewChild(1);
        tree2->addNewChild(3);
        tree2->addNewChild(3);
        tree2->addNewChild(4);
        tree2->addNewChild(4);
        tree2->addNewChild(2);
        tree2->addNewChild(2);
        tree2->addNewChild(2);
        tree2->addNewChild(0);
        TS_ASSERT_EQUALS(tree1->topologyTreeCompare(tree2),tree2->topologyTreeCompare(tree1));
        TS_ASSERT_EQUALS(tree1->topologyTreeCompare(tree2),false);
        delete tree1;
        delete tree2;
    }


};
#endif /* TESTDIFFERENTIATIONTREE_HPP_ */
