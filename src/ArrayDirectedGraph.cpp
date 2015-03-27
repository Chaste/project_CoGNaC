#include "ArrayDirectedGraph.hpp"
#include "RandomNumberGenerator.hpp"
#include <cassert>
#include <math.h>
#include <iostream>

ArrayDirectedGraph::ArrayDirectedGraph(unsigned size) :
mSize (size)
{
   if (mSize == 0) EXCEPTION("Size of the graph must be > 0.");
   mNumberOfEdges = 0;
   mpVertices = new GraphNode*[mSize];
   mpInputEdges = new unsigned[mSize];
   mpOutputEdges = new unsigned[mSize];
   for (unsigned index=0; index<mSize; index++){
       mpVertices[index] = new GraphNode(index);
       mpInputEdges[index] = 0;
       mpOutputEdges[index] = 0;
   }
}

ArrayDirectedGraph::~ArrayDirectedGraph()
{
    for (unsigned i=0; i<mSize; i++)
    {
        delete mpVertices[i];
    }
    delete[] mpVertices;
    delete[] mpInputEdges;
    delete[] mpOutputEdges;
}

void ArrayDirectedGraph::randomNetworkGenerator(unsigned average_inputs_number_per_node){
    for (unsigned node=0; node < mSize; node ++)
    {
        unsigned edge_counter = 0;
        while (edge_counter < average_inputs_number_per_node)
        {
            bool added = false;
            while (!added)
            {
                unsigned target_node = RandomNumberGenerator::Instance()->randMod(mSize);
                if ((unsigned) target_node != node) //avoiding self edge
                {
                    if (RandomNumberGenerator::Instance()->ranf() > 0.5)
                    {
                        added = addEdgeById(node, target_node);
                    } else {
                        added = addEdgeById(target_node, node);
                    }
                }
            }
            edge_counter++;
        }
    }
}

void ArrayDirectedGraph::albertBarabasiGenerator(unsigned average_inputs_number_per_node)
{
    unsigned initial_nodes = average_inputs_number_per_node + 1;
    for(unsigned i = 0; i < initial_nodes - 1; i++)
    {
        for(unsigned j = i + 1; j < initial_nodes; j++)
        {
            addEdgeById(i,j);
            addEdgeById(j,i);
        }
    }
    for(unsigned i = initial_nodes; i < mSize; i++)
    {
        for(unsigned j = 0; j < average_inputs_number_per_node; j++)
        {
            unsigned node_0 = 0;
            unsigned node_1 = 0;
            do
            {
                unsigned preferential_attachment = 0;
                double random_number = RandomNumberGenerator::Instance()->ranf();
                double degree = (mpInputEdges[preferential_attachment] + mpOutputEdges[preferential_attachment])
                        * 1.0 / mNumberOfEdges;
                while(random_number > degree){
                    preferential_attachment++;
                    degree += (mpInputEdges[preferential_attachment] + mpOutputEdges[preferential_attachment])
                                * 1.0 / mNumberOfEdges;
                }
                if(RandomNumberGenerator::Instance()->ranf() > 0.5)
                {
                    node_0 = i;
                    node_1 = preferential_attachment;
                } else
                {
                    node_0 = preferential_attachment;
                    node_1 = i;
                }
            } while(!addEdgeById(node_0, node_1));
        }
    }
}

void ArrayDirectedGraph::erdosRenyiGenerator(unsigned average_inputs_number_per_node)
{
    /*
     * It can generate self edges.
     */
    for(unsigned i = 0; i < average_inputs_number_per_node * mSize; i++)
    {
        while (!addEdgeById(RandomNumberGenerator::Instance()->randMod(mSize),RandomNumberGenerator::Instance()->randMod(mSize)));
    }
}

bool ArrayDirectedGraph::addEdgeById(unsigned input_id, unsigned output_id){
    if (! (input_id < mSize && output_id < mSize)) EXCEPTION("Source ID or target ID out of range.");
    if (mpVertices[output_id]->addInputEdgeById(input_id)){
        mpInputEdges[output_id]++;
        mpOutputEdges[input_id]++;
        mNumberOfEdges ++;
        return true;
    }
    return false;
}

void ArrayDirectedGraph::sortGraph(){
    /* sort the graph in decreasing order respect number of
     * input nodes, starting from node zero.
     * We use bubble sort, because size of the graph is
     * usually not too big.
     */
    unsigned i=0;
    for (i=0; i<mSize - 1; i++)
    {
        for (unsigned j=i+1; j<mSize; j++)
        {
            if(mpInputEdges[i] < mpInputEdges[j])
            {
                GraphNode* temp_node = mpVertices[i];
                mpVertices[i] = mpVertices[j];
                mpVertices[j] = temp_node;


                unsigned temp_num = mpInputEdges[i];
                mpInputEdges[i] = mpInputEdges[j];
                mpInputEdges[j] = temp_num;

                temp_num = mpOutputEdges[i];
                mpOutputEdges[i] = mpOutputEdges[j];
                mpOutputEdges[j] = temp_num;
            }
        }
        mpVertices[i]->setId(i);
    }
    mpVertices[i]->setId(i);
}

void ArrayDirectedGraph::printGraphToGmlFile(std::string directory, std::string filename) const
{
    if (directory.empty())
        EXCEPTION("Directory name not valid.");

    if (filename.size() < 5 || filename.compare(filename.size()-4,4,".gml") != 0)
        EXCEPTION("File path not valid. It must terminate with '.gml' extension.");

    OutputFileHandler handler(directory,false);
   	out_stream p_file = handler.OpenOutputFile(filename);

	*p_file << "graph [\n";
	*p_file << " comment \"Graph generated by Chaste\"\n";
	*p_file << " directed 1\n";
	unsigned square_position = 0;
	float x = 0.0, y = 0.0, step = 80.0, square_root = 0.0;
	square_root = sqrt((float) mSize);
	if (ceilf(square_root) == square_root) square_position = (unsigned) square_root;
	else square_position = (unsigned) square_root + 1;
	/* for each node */
	for (unsigned i=0; i < mSize; i++)
	{
		if (i > 0 && i % square_position == 0)
		{
			y += step;
			x = 0.0;
		}
		*p_file << " node [\n";
		*p_file << "  id " << i+1 << "\n";
		*p_file << "  label \"" << i+1 << "\"\n";
		*p_file << "  graphics [\n   x " << x << "\n   y " << y << "\n";
		*p_file << "   w 40.0\n   h 40.0\n   type \"ellipse\"\n  ]\n";
		*p_file << "  LabelGraphics [\n   text \"" << i+1 << "\"\n  ]\n";
		*p_file << " ]\n";
		x += step;
	}
	/* for each edge (we take the inputs of every edge) */
	for (unsigned i=0; i < mSize; i++)
	{
		if (mpVertices[i]->getIncomingVerticesSize() > 0)
		{
			std::vector<unsigned> input_vertices = mpVertices[i]->getIncomingVerticesId();
			for (unsigned j=0; j<input_vertices.size(); j++)
			{
				*p_file << " edge [\n";
				*p_file << "  source " << input_vertices[j] + 1 << "\n";
				*p_file << "  target " << i+1 << "\n";
				*p_file << "  label \"" << input_vertices[j] + 1 << " -> "<< i+1 << "\"\n";
				*p_file << "  graphics [\n   arrow \"last\"\n  ]\n";
				*p_file << " ]\n";
			}
		}

	}
	*p_file << "]";
	p_file->close();
}

void ArrayDirectedGraph::printGraphToSifFile(std::string directory, std::string filename) const
{
    if (directory.empty())
        EXCEPTION("Directory name not valid.");

    if (filename.size() < 5 || filename.compare(filename.size()-4,4,".sif") != 0)
        EXCEPTION("File path not valid. It must terminate with '.sif' extension.");

    OutputFileHandler handler(directory,false);
   	out_stream p_file = handler.OpenOutputFile(filename);

	for (unsigned i=0; i < mSize; i++)
	{
		if (mpVertices[i]->getIncomingVerticesSize() > 0)
		{
			std::vector<unsigned> input_vertices = mpVertices[i]->getIncomingVerticesId();
			for (unsigned j=0; j<input_vertices.size(); j++)
			{
			   *p_file << "node" << input_vertices[j] + 1 << "\tDirectedEdge\tnode" << i+1 << "\n";
			}
		}

	}
	p_file->close();
}

void ArrayDirectedGraph::printGraphToDotFile(std::string directory, std::string filename) const
{
    if (directory.empty())
        EXCEPTION("Directory name not valid.");

    if (filename.size() < 5 || filename.compare(filename.size()-4,4,".dot") != 0)
        EXCEPTION("File path not valid. It must terminate with '.dot' extension.");

    OutputFileHandler handler(directory,false);
   	out_stream p_file = handler.OpenOutputFile(filename);

	*p_file << "digraph boolean_network_graph {\n";
	for (unsigned i=0; i < mSize; i++)
	{
		if (mpVertices[i]->getIncomingVerticesSize() > 0)
		{
			std::vector<unsigned> input_vertices = mpVertices[i]->getIncomingVerticesId();
			for (unsigned j=0; j<input_vertices.size(); j++)
			{
				*p_file << " " << input_vertices[j] + 1 << " -> " << i+1 << ";\n";
			}
		}
	}
	*p_file << "}";
	p_file->close();
}

void ArrayDirectedGraph::printGraph() const
{
    /* console output */
    for(unsigned node = 0; node < mSize; node++)
    {
        std::cout << "Node ID: " << mpVertices[node]->getId() << "\n";
        if (mpVertices[node]->getIncomingVerticesSize() > 0)
        {
            std::vector<unsigned> input_vertices = mpVertices[node]->getIncomingVerticesId();
            for (unsigned input_node=0; input_node<input_vertices.size(); input_node++)
            {
                std::cout << "---> Input node: " << input_vertices[input_node] << "\n";
            }
        }
    }
}

unsigned ArrayDirectedGraph::getSize() const
{
    return mSize;
}

unsigned* ArrayDirectedGraph::getInputEdges() const
{
    /* Return a copy of mpInputEdges */
    unsigned* input_edges = new unsigned[mSize];
    std::copy(mpInputEdges, mpInputEdges+mSize, input_edges);

    return input_edges;
}

unsigned* ArrayDirectedGraph::getOutputEdges() const
{
    /* Return a copy of mpOutputEdges */
    unsigned* output_edges = new unsigned[mSize];
    std::copy(mpOutputEdges, mpOutputEdges+mSize, output_edges);

    return output_edges;
}

std::vector<unsigned> ArrayDirectedGraph::getIncomingVerticesById(unsigned node) const{
    if (node >= mSize)
        EXCEPTION("Error, node id must be < size of the graph.");
    return mpVertices[node]->getIncomingVerticesId();
}

unsigned* ArrayDirectedGraph::getIncomingVerticesArrayIdById(unsigned node) const{
    if (node >= mSize)
        EXCEPTION("Error, node id must be < size of the graph.");
    return mpVertices[node]->getInputVerticesId();
}

unsigned ArrayDirectedGraph::getIncomingVerticesNumberById(unsigned node) const{
    if (node >= mSize)
        EXCEPTION("Error, node id must be < size of the graph.");
    return mpVertices[node]->getIncomingVerticesSize();
}

void ArrayDirectedGraph::setIncomingVerticesById(
        const unsigned* p_inputs_array,
        const unsigned vertices_number,
        const unsigned node
        )
{
    if (node >= mSize)
        EXCEPTION("Error, node id must be < size of the graph.");
    if (!p_inputs_array)
        EXCEPTION("Error, pointer is NULL.");
    mpVertices[node]->setIncomingVerticesId(p_inputs_array, vertices_number);
}
