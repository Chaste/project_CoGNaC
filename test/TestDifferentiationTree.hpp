#ifndef TESTDIFFERENTIATIONTREE_HPP_
#define TESTDIFFERENTIATIONTREE_HPP_

/*
 * = Testing the class {{{DifferentiationTree}}} =
 *
 * == Introduction ==
 *
 * EMPTYLINE
 *
 * This class is used to test that the classes
 * {{{DifferentiationTree}}} and {{{DifferentiationTreeNode}}}
 * are implemented correctly.
 *
 * == Including header files ==
 *
 * EMPTYLINE
 *
 * We begin by including the necessary header files.
*/

#include <cxxtest/TestSuite.h>
#include <set>
#include <vector>

#include "DifferentiationTree.hpp"
#include "DifferentiationTreeNode.hpp"

class TestDifferentiationTree : public CxxTest::TestSuite
{
public:
	/*
	 * == Testing simple methods and functions. ==
	 *
	 * EMPTYLINE
	 *
	 */
    void testDifferentiationTree()
    {
    	/* First of all we instantiate three sets and two probability
    	 * distributions, which we use for the tree declaration later. */

        std::set<unsigned> set0;
        set0.insert(0);
        set0.insert(1);
        std::set<unsigned> set1;
        set1.insert(0);
        std::set<unsigned> set2;
        set2.insert(1);
        std::vector<double> prob_distribution1;
        prob_distribution1.push_back(.5);
        prob_distribution1.push_back(.5);
        std::vector<double> prob_distribution2;
        prob_distribution2.push_back(1);

        /* We instantiate a tree, whose root has associated the ergodic set {0,1}: */
        DifferentiationTree* tree1 = new DifferentiationTree(8, set0, prob_distribution1, prob_distribution1);
        TS_ASSERT(tree1->getLeaves().size() == 1);

        /* We add two nodes in the tree, with ergodic set {0} and {1} respectively. Then we check that
         * each node has been added correctly, verifying the topological properties of the tree. */
        unsigned parent = tree1->searchParentNode(set1);
        tree1->addNewChild(parent, 4, set1, prob_distribution2, prob_distribution2);
        TS_ASSERT(tree1->getLeaves().size() == 1);

        parent = tree1->searchParentNode(set2);
        tree1->addNewChild(parent, 3, set2, prob_distribution2, prob_distribution2);
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

        /* We also test the {{{normaliseLength}}} and {{{colorise}}} method. */
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
	/*
	 * == Comparing trees. ==
	 *
	 * EMPTYLINE
	 *
	 */
    void testComparingTrees()
    {
    	/* We instantiate six trees and test if the method {{{topologyTreeCompare}}}
    	 * runs correctly. Remember that the parameter in function addNewChild is
    	 * always the parent id in this case. */
        DifferentiationTree* tree1 = new DifferentiationTree();
        tree1->addNewChild(0);//1
        tree1->addNewChild(0);//2
        tree1->addNewChild(2);//3
        tree1->addNewChild(2);//4
        tree1->addNewChild(3);//5
        tree1->addNewChild(3);//6
        tree1->addNewChild(4);//7

        DifferentiationTree* tree2 = new DifferentiationTree();
        tree2->addNewChild(0);//1
        tree2->addNewChild(0);//2
        tree2->addNewChild(2);//3
        tree2->addNewChild(2);//4
        tree2->addNewChild(2);//5
        tree2->addNewChild(3);//6
        tree2->addNewChild(4);//7

        TS_ASSERT_EQUALS(tree1->topologyTreeCompare(tree2),tree2->topologyTreeCompare(tree1));
        TS_ASSERT_EQUALS(tree1->topologyTreeCompare(tree2),false);
        delete tree1;
        delete tree2;

        DifferentiationTree* tree3 = new DifferentiationTree();
        tree3->addNewChild(0);//1
        tree3->addNewChild(0);//2
        tree3->addNewChild(1);//3
        tree3->addNewChild(1);//4
        tree3->addNewChild(1);//5
        tree3->addNewChild(2);//6
        tree3->addNewChild(2);//7
        tree3->addNewChild(2);//8
        tree3->addNewChild(6);//9
        tree3->addNewChild(6);//10
        tree3->addNewChild(8);//11
        tree3->addNewChild(8);//12

        DifferentiationTree* tree4 = new DifferentiationTree();
        tree4->addNewChild(0);//1
        tree4->addNewChild(0);//2
        tree4->addNewChild(1);//3
        tree4->addNewChild(1);//4
        tree4->addNewChild(1);//5
        tree4->addNewChild(3);//6
        tree4->addNewChild(3);//7
        tree4->addNewChild(4);//8
        tree4->addNewChild(4);//9
        tree4->addNewChild(2);//10
        tree4->addNewChild(2);//11
        tree4->addNewChild(2);//12

        TS_ASSERT_EQUALS(tree3->topologyTreeCompare(tree4),tree4->topologyTreeCompare(tree3));
        TS_ASSERT_EQUALS(tree3->topologyTreeCompare(tree4),true);
        delete tree3;
        delete tree4;

        DifferentiationTree* tree5 = new DifferentiationTree();
        tree5->addNewChild(0);//1
        tree5->addNewChild(0);//2
        tree5->addNewChild(1);//3
        tree5->addNewChild(1);//4
        tree5->addNewChild(1);//5
        tree5->addNewChild(2);//6
        tree5->addNewChild(2);//7
        tree5->addNewChild(2);//8
        tree5->addNewChild(6);//9
        tree5->addNewChild(6);//10
        tree5->addNewChild(8);//11
        tree5->addNewChild(8);//12
        tree5->addNewChild(2);//13

        DifferentiationTree* tree6 = new DifferentiationTree();
        tree6->addNewChild(0);//1
        tree6->addNewChild(0);//2
        tree6->addNewChild(1);//3
        tree6->addNewChild(1);//4
        tree6->addNewChild(1);//5
        tree6->addNewChild(3);//6
        tree6->addNewChild(3);//7
        tree6->addNewChild(4);//8
        tree6->addNewChild(4);//9
        tree6->addNewChild(2);//10
        tree6->addNewChild(2);//11
        tree6->addNewChild(2);//12
        tree6->addNewChild(0);//13

        TS_ASSERT_EQUALS(tree5->topologyTreeCompare(tree6),tree6->topologyTreeCompare(tree5));
        TS_ASSERT_EQUALS(tree5->topologyTreeCompare(tree6),false);
        delete tree5;
        delete tree6;
    }
};
#endif /* TESTDIFFERENTIATIONTREE_HPP_ */
