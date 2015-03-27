#include "GraphNode.hpp"
#include <boost/shared_ptr.hpp>

GraphNode::GraphNode(unsigned node_id) :
mId (node_id)
{
}

GraphNode::~GraphNode()
{
}

bool GraphNode::addInputEdgeById(unsigned input_id){
    if (std::find(mIncomingVertices.begin(), mIncomingVertices.end(), input_id)
            == mIncomingVertices.end())
    {
        mIncomingVertices.push_back(input_id);
        return true;
    }
    return false;
}

bool GraphNode::removeInputEdgeById(unsigned input_node_id){
    std::vector<unsigned>::iterator position =
                std::find(mIncomingVertices.begin(), mIncomingVertices.end(), input_node_id);
    if (position != mIncomingVertices.end())
    {
        mIncomingVertices.erase(position);
        return true;
    }
    return false;
}

bool GraphNode::existsEdgeFrom(unsigned input_node_id) const{
    return std::find (
            mIncomingVertices.begin(), mIncomingVertices.end(),
            input_node_id ) != mIncomingVertices.end();
}

void GraphNode::printNode() const
{
    std::cout << "ID: " << mId << "\n";
    for (unsigned i=0; i<mIncomingVertices.size(); i++)
    {
        std::cout << "--. input node: " << mIncomingVertices.at(i) << "\n";
    }
}

unsigned GraphNode::getId() const{
    return mId;
}

unsigned GraphNode::getIncomingVerticesSize() const{
    return mIncomingVertices.size();
}

std::vector<unsigned> GraphNode::getIncomingVerticesId() const{
    return mIncomingVertices;
}

unsigned* GraphNode::getInputVerticesId() const{
    if (mIncomingVertices.size() > 0){
        unsigned* p_input_vertices = new unsigned[mIncomingVertices.size()];
        std::copy(mIncomingVertices.begin(), mIncomingVertices.end(), p_input_vertices);
        return p_input_vertices;
    }
    return NULL;
}

void GraphNode::setId(unsigned id)
{
    mId = id;
}

void GraphNode::setIncomingVerticesId(
        const unsigned* p_input_vertices,
        const unsigned vertices_number
        )
{
    if (!p_input_vertices)
        EXCEPTION("Error, pointer is NULL.");    mIncomingVertices.clear();
    //we use set for avoid duplicate incoming edges.
    std::set<unsigned> set_vertices (
            p_input_vertices,
            p_input_vertices+vertices_number);
    std::set<unsigned>::iterator set_it;
    for (set_it = set_vertices.begin(); set_it != set_vertices.end(); ++set_it)
    {
        mIncomingVertices.push_back(*set_it);
    }

}
