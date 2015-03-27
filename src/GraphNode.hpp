#ifndef GRAPHNODE_HPP_
#define GRAPHNODE_HPP_

#include "Exception.hpp"
#include <vector>
#include <algorithm>
#include <set>

class GraphNode
{
private:
    /** Id of the node */
    unsigned mId;

    /** Incoming vertices of the node */
    std::vector<unsigned> mIncomingVertices;

public:
    /**
     * Constructor: create a GraphNode object
     */
    GraphNode(unsigned node_id);

    /** Destruptor */
    ~GraphNode();

    /**
     * Add an input to this node.
     *
     * @param input_node_id
     *
     * @return true if the edge doesn't exist yet.
     */
    bool addInputEdgeById(unsigned input_node_id);

    /**
     * Remove an input to this node.
     *
     * @param input_node_id
     *
     * @return true if the edge already exists.
     */
    bool removeInputEdgeById(unsigned input_node_id);

    /**
     * Check if exists an input to this node from another node.
     *
     * @param input_node_id
     *
     * @return true if the edge exists.
     */
    bool existsEdgeFrom(unsigned input_node_id) const;

    /**
     * Print node information in the console output.
     */
    void printNode() const;

    /**
     * Getter of mId
     *
     * @return mId
     */
    unsigned getId() const;

    /**
     * @return the number of incoming vertices.
     */
    unsigned getIncomingVerticesSize() const;

    /**
     * Getter of mIncomingVertices
     *
     * @return mIncomingVertices
     */
    std::vector<unsigned> getIncomingVerticesId() const;

    /**
     * Getter of mIncomingVertices. Returns a copy of mIncomingVertices as
     * array.
     * User must delete the object explicitly after use.
     *
     * @return incoming vertices as array
     */
    unsigned* getInputVerticesId() const;

    /**
     * Setter of mId
     *
     * @param id
     */
    void setId(unsigned id);

    /**
     * Setter the incoming vertices from an array.
     *
     * @param p_input_vertices
     * @param vertices_number
     */
    void setIncomingVerticesId(const unsigned* p_input_vertices, const unsigned vertices_number);
};

#endif /* GRAPHNODE_HPP_ */
