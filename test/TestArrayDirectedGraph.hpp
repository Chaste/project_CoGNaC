#ifndef TESTARRAYDIRECTEDGRAPH_HPP_
#define TESTARRAYDIRECTEDGRAPH_HPP_

#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "SmartPointers.hpp"

#include "ArrayDirectedGraph.hpp"
#include "GraphNode.hpp"


class TestArrayDirectedGraph  : public CxxTest::TestSuite
{
public:
    void testArrayDirectedGraph()
    {
        ArrayDirectedGraph graph(4);
        TS_ASSERT(graph.getSize() == 4);

        TS_ASSERT(graph.addEdgeById(1,0));
        TS_ASSERT(graph.addEdgeById(0,0));
        TS_ASSERT(graph.addEdgeById(0,3));
        TS_ASSERT(graph.addEdgeById(3,3));
        TS_ASSERT(graph.addEdgeById(2,3));

        TS_ASSERT(!graph.addEdgeById(1,0));

        unsigned input_nodes[] = {1,0};
        std::vector<unsigned> inc_vertices = graph.getIncomingVerticesById(0);
        TS_ASSERT_EQUALS(graph.getIncomingVerticesNumberById(0), 2u);
        for (unsigned i=0;i<2;i++)
        {
            TS_ASSERT_EQUALS(inc_vertices[i],input_nodes[i]);
        }

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

        inputs = graph.getInputEdges();
        TS_ASSERT(inputs[0] == 2);
        delete[] inputs;

        inputs = graph.getOutputEdges();
        TS_ASSERT(inputs[2] == 1);
        delete[] inputs;


        graph.sortGraph();
        TS_ASSERT(graph.getIncomingVerticesNumberById(0) == 2);
        TS_ASSERT(graph.getIncomingVerticesNumberById(1) == 2);
        TS_ASSERT(graph.getIncomingVerticesNumberById(2) == 0);
        TS_ASSERT(graph.getIncomingVerticesNumberById(3) == 0);

        ArrayDirectedGraph randomGraph(6);
        randomGraph.randomNetworkGenerator(2);
        randomGraph.sortGraph();

        for (unsigned i=0; i<5; i++)
        {
            TS_ASSERT(randomGraph.getIncomingVerticesNumberById(i) >=
                    randomGraph.getIncomingVerticesNumberById(i+1));
        }

    }
};

#endif /* TESTARRAYDIRECTEDGRAPH_HPP_ */
