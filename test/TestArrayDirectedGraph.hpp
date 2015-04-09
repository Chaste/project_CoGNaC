#ifndef TESTARRAYDIRECTEDGRAPH_HPP_
#define TESTARRAYDIRECTEDGRAPH_HPP_

/*
 * = Testing the class {{{ArrayDirectedGraph}}} =
 * This class is used to test that the class {{{ArrayDirectedGraph}}} is implemented correctly.
 *
 * === Including header files ===
 * We begin by including the necessary header files.
*/

#include <cxxtest/TestSuite.h>

#include "ArrayDirectedGraph.hpp"
#include "GraphNode.hpp"


class TestArrayDirectedGraph  : public CxxTest::TestSuite
{
public:
	/*
	 * == Testing methods and functions. ==
	 */
    void testArrayDirectedGraph()
    {
    	/* We construct and initialise a simple graph with 4 nodes.: */
        ArrayDirectedGraph graph(4);
        TS_ASSERT(graph.getSize() == 4);

        /* We insert five edges, testing that the addition of each edge is done properly (each addition returns a {{{true}}} value): */
        TS_ASSERT(graph.addEdgeById(1,0));
        TS_ASSERT(graph.addEdgeById(0,0));
        TS_ASSERT(graph.addEdgeById(0,3));
        TS_ASSERT(graph.addEdgeById(3,3));
        TS_ASSERT(graph.addEdgeById(2,3));

        /* Also we test that if we try to insert an already added arc, it returns false: */
        TS_ASSERT(!graph.addEdgeById(1,0));

        /* We take a node and check if the incoming edges are the same we have added before. */
        unsigned input_nodes[] = {1,0};
        std::vector<unsigned> inc_vertices = graph.getIncomingVerticesById(0);
        TS_ASSERT_EQUALS(graph.getIncomingVerticesNumberById(0), 2u);
        for (unsigned i=0;i<2;i++)
        {
            TS_ASSERT_EQUALS(inc_vertices[i],input_nodes[i]);
        }

        /* We also add two edges to a node using {{setIncomingVerticesById}}}.
         * and test that the addition is done correctly. */
        unsigned* inputs = new unsigned[2];
        inputs[0] = 1;
        inputs[1] = 2;

        graph.setIncomingVerticesById(inputs, 2, 3);
        input_nodes[0] = 1;
        input_nodes[1] = 2;
        inc_vertices = graph.getIncomingVerticesById(3);
        TS_ASSERT_EQUALS(graph.getIncomingVerticesNumberById(3), 2u);
        for (unsigned i=0;i<2;i++)
        {
            TS_ASSERT_EQUALS(inc_vertices[i],input_nodes[i]);
        }
        delete[] inputs;

        /* We also test the function {{{getInputEdges}}} and {{{getOutputEdges}}}: */
        inputs = graph.getInputEdges();
        TS_ASSERT(inputs[0] == 2);
        delete[] inputs;

        inputs = graph.getOutputEdges();
        TS_ASSERT(inputs[2] == 1);
        delete[] inputs;

        /* We sort the graph in decreasing order respect number of
         * input nodes, starting from node zero. Then, we test the
         * number of incoming arcs of each node. */
        graph.sortGraph();
        TS_ASSERT(graph.getIncomingVerticesNumberById(0) == 2);
        TS_ASSERT(graph.getIncomingVerticesNumberById(1) == 2);
        TS_ASSERT(graph.getIncomingVerticesNumberById(2) == 0);
        TS_ASSERT(graph.getIncomingVerticesNumberById(3) == 0);

        /* Finally, we generate a random graph and we sort it.
         * Then, we test the sort property.
         */
        ArrayDirectedGraph randomGraph(6);
        randomGraph.randomNetworkGenerator(2);
        randomGraph.sortGraph();

        for (unsigned i=0; i<randomGraph.getSize() - 1; i++)
        {
            TS_ASSERT(randomGraph.getIncomingVerticesNumberById(i) >=
                    randomGraph.getIncomingVerticesNumberById(i+1));
        }

    }
};

#endif /* TESTARRAYDIRECTEDGRAPH_HPP_ */
