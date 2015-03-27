#ifndef ARRAYDIRECTEDGRAPH_HPP_
#define ARRAYDIRECTEDGRAPH_HPP_

#include "GraphNode.hpp"
#include "Exception.hpp"
#include "OutputFileHandler.hpp"
#include <vector>
#include <iostream>
#include <fstream>


class ArrayDirectedGraph
{
private:
    /** Array of GraphNode. It contains every node in the Graph */
    GraphNode** mpVertices;

    /** Array that is a counter of input edges for every node */
    unsigned* mpInputEdges;

    /** Array that is a counter of output edges for every node */
    unsigned* mpOutputEdges;

    /** Number of vertices of the graph */
    unsigned mSize;

    /** Number of Edges of the graph */
    unsigned mNumberOfEdges;
public:

    /**
     * Costructor: create the data structures
     *
     * @param size number of vertices of the graph
     */
    ArrayDirectedGraph(unsigned size);

    /**
     * Destruptor: deallocate nodes and data structures
     */
    ~ArrayDirectedGraph();

    /**
     * Completely random graph generator.
     *
     * @params average_inputs_number_per_node
     */
    void randomNetworkGenerator(unsigned average_inputs_number_per_node);

    /**
     * Albert Barabasi random graph generator.
     *
     * @params average_inputs_number_per_node
     */
    void albertBarabasiGenerator(unsigned average_inputs_number_per_node);

    /**
     * Erdos Renyi random graph generator.
     *
     * @params average_inputs_number_per_node
     */
    void erdosRenyiGenerator(unsigned average_inputs_number_per_node);

    /**
     * Add an edge between two nodes.
     *
     * @param input id
     * @param output id
     *
     * @return true if the edge doesn't exist yet
     */
    bool addEdgeById(unsigned input, unsigned output);

    /**
     * Method for sort the graph decreasing by the number
     * of input vertices. Node zero will have the maximum
     * of input vertices in the graph, and so on.
     */
    void sortGraph();

    /**
     * Get the incoming vertices of a given node as a vector of ids.
     *
     * @param node id of the node
     *
     * @return a vector of incoming vertices ids
     */
    std::vector<unsigned> getIncomingVerticesById(unsigned node) const;

    /**
     * Get the incoming vertices of a given node as a array of ids.
     * User must delete the object explicitly after use.
     *
     * @param node id of the node
     *
     * @return an array of incoming vertices ids
     */
    unsigned* getIncomingVerticesArrayIdById(unsigned node) const;

    /**
     * Get the number of incoming vertices of a given node.
     *
     * @param node id of the node
     *
     * @return the number of incoming vertices
     */
    unsigned getIncomingVerticesNumberById(unsigned node) const;

    /**
     * Set incoming vertices of a given node, using an id array.
     *
     * @param p_inputs_array array of ids
     * @param vertices_number number of vertices in the array
     * @param node id of the node.
     */
    void setIncomingVerticesById(const unsigned* p_inputs_array,const unsigned vertices_number,const unsigned node);

    /**
     * Save graph in a .gml file.
     *
     * @param directory the name of the subfolder of testoutput.
     * @param filename the file name of the output file.
     */
    void printGraphToGmlFile(const std::string directory, const std::string filename) const;

    /**
     * Save graph in a .sif file.
     *
     * @param directory the name of the subfolder of testoutput.
     * @param filename the file name of the output file.
     */
    void printGraphToSifFile(const std::string directory, const std::string filename) const;

    /**
     * Save graph in a .dot file.
     *
     * @param directory the name of the subfolder of testoutput.
     * @param filename the file name of the output file.
     */
    void printGraphToDotFile(const std::string directory, const std::string filename) const;

    /**
     * Print graph in the console output.
     */
    void printGraph() const;

    /**
     * Getter of mSize
     *
     * @return mSize
     */
    unsigned getSize() const;

    /**
     * Getter of mInputEdges. Returns a copy of mInputEdges.
     * User must delete the object explicitly after use.
     *
     * @return mInputEdges
     */
    unsigned* getInputEdges() const;

    /**
     * Getter of mOutputEdges. Returns a copy of mInputEdges.
     * User must delete the object explicitly after use.
     *
     * @return mOutputEdges
     */
    unsigned* getOutputEdges() const;
};

#endif /* ARRAYDIRECTEDGRAPH_HPP_ */
