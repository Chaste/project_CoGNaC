#ifndef TESTGRAPHNODE_HPP_
#define TESTGRAPHNODE_HPP_

/*
 * = Testing the class {{{GraphNode}}} =
 *
 * == Introduction ==
 *
 * EMPTYLINE
 *
 * This class is used to test that the class {{{GraphNode}}}
 * is implemented correctly.
 *
 * == Including header files ==
 * We begin by including the necessary header files.
*/

#include <cxxtest/TestSuite.h>

#include "GraphNode.hpp"

class TestGraphNode  : public CxxTest::TestSuite
{
public:
	/*
	 * == Testing methods and functions. ==
	 *
	 * EMPTYLINE
	 */
	void testGraphNode()
    {
		/* We start instantiating four nodes. */
        GraphNode node0(0);
        GraphNode node1(1);
        GraphNode node2(2);
        GraphNode node3(3);

        /* We add some input edge to each node, obtaining a graph. This graph
         * is not managed directly using {{{ArrayDirectedGraph}}} class, but
         * but used in order to test data saved on each node. */
        node0.addInputEdgeById(node1.getId());
        node0.addInputEdgeById(node3.getId());
        node1.addInputEdgeById(node1.getId());
        node2.addInputEdgeById(node3.getId());

        /* We test if edges exist. */
        TS_ASSERT(node0.existsEdgeFrom(node1.getId()));
        TS_ASSERT(!node3.existsEdgeFrom(node1.getId()));
        TS_ASSERT(node3.getIncomingVerticesSize() == 0);

        /* We define an array of node ids and associate
         * only five of these nodes as incoming arcs of the node 3.
         */
        unsigned array[]= {89,2,3,4,5,6,7,8,9};
        unsigned* p_a = array;
        node3.setIncomingVerticesId(p_a, 5);
        TS_ASSERT(node3.getIncomingVerticesSize() == 5);

        /* Test if {{{removeInputEdgeById}}} works correctly: */
        node3.removeInputEdgeById(node2.getId());
        TS_ASSERT(node3.getIncomingVerticesSize() == 4);

        node3.removeInputEdgeById(100);
        TS_ASSERT(node3.getIncomingVerticesSize() == 4);

        /* Also we test if the method {{{setId}}} works as expected. */
        node3.setId(99);
        TS_ASSERT(node3.getId() == 99);

        /* Finally, we test the method {{{getIncomingVerticesId}}}
         * by evaluating it on the node 1.
         */
        std::vector<unsigned> inc_vertices = node1.getIncomingVerticesId();
        TS_ASSERT(inc_vertices.size() > 0);
        TS_ASSERT(inc_vertices.at(0) == 1);
    }
};


#endif /* TESTGRAPHNODE_HPP_ */
