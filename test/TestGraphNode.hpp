#ifndef TESTGRAPHNODE_HPP_
#define TESTGRAPHNODE_HPP_

#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "SmartPointers.hpp"

#include "GraphNode.hpp"


class TestGraphNode  : public CxxTest::TestSuite
{
public:
	void testGraphNode()
    {
        GraphNode node0(0);
        GraphNode node1(1);
        GraphNode node2(2);
        GraphNode node3(3);

        node0.addInputEdgeById(node1.getId());
        node0.addInputEdgeById(node3.getId());
        node1.addInputEdgeById(node1.getId());
        node2.addInputEdgeById(node3.getId());

        TS_ASSERT(node0.existsEdgeFrom(node1.getId()));
        TS_ASSERT(!node3.existsEdgeFrom(node1.getId()));
        TS_ASSERT(node3.getIncomingVerticesSize() == 0);

        unsigned array[]= {89,2,3,4,5,6,7,8,9};
        unsigned* p_a = array;
        node3.setIncomingVerticesId(p_a, 5);
        TS_ASSERT(node3.getIncomingVerticesSize() == 5);

        node3.removeInputEdgeById(node2.getId());
        TS_ASSERT(node3.getIncomingVerticesSize() == 4);

        node3.removeInputEdgeById(100);
        TS_ASSERT(node3.getIncomingVerticesSize() == 4);
        node3.setId(99);
        TS_ASSERT(node3.getId() == 99);

        std::vector<unsigned> inc_vertices = node1.getIncomingVerticesId();
        TS_ASSERT(inc_vertices.size() > 0);
        TS_ASSERT(inc_vertices.at(0) == 1);
    }
};


#endif /* TESTGRAPHNODE_HPP_ */
