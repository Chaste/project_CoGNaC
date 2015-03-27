#include "RandomBooleanNetwork.hpp"
#include "RandomNumberGenerator.hpp"
#include <math.h>
#include <cassert>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

/**
 * Default bddallsathandler defined by buddy for allsat()
 */
void allsatHandlerPrint(char *varset,int size);

/**
 * bddallsathandler defined for print states (as string)
 * in the stringBuffer global variable for allsat(). Used
 * for printing network in .net or .cnet files.
 */
void allsatHandlerPrintToStringBuffer(char *varset,int size);

/** String buffer used for print functions */
std::vector<std::string> stringBuffer;

RandomBooleanNetwork::RandomBooleanNetwork(unsigned nodes_number, unsigned avarage_inputs_per_node,
        bool scale_free, double probability_canalyzing_function) :
mNodesNumber(nodes_number),
mAverageInputsPerNode(avarage_inputs_per_node)
{
    if (mNodesNumber == 0 || mAverageInputsPerNode == 0 || mAverageInputsPerNode >= mNodesNumber)
        EXCEPTION("Constructor parameters not valid.");
    if (probability_canalyzing_function < 0.0 || probability_canalyzing_function > 1.0)
        EXCEPTION("Range of a probability must be between 0 and 1.");

    mpRbnGraph = new ArrayDirectedGraph(mNodesNumber);
    if (scale_free)
    {
       mpRbnGraph->albertBarabasiGenerator(mAverageInputsPerNode);
    } else {
       //mpRbnGraph->randomNetworkGenerator(mAverageInputsPerNode);
       mpRbnGraph->erdosRenyiGenerator(mAverageInputsPerNode);
    }
    //mpRbnGraph->sortGraph();
    initBinaryDecisionDiagram();

    mpReverseTransitionFunction = new bdd();
    *mpReverseTransitionFunction = bddtrue;
    mpNodeFunction = new bdd[mNodesNumber];
    mpNodeNthFunction = new bdd[mNodesNumber];
    mpVariables = new bdd[mNodesNumber];
    mpNextVariables = new bdd[mNodesNumber];
    for (unsigned i=0; i<mNodesNumber; i++)
    {
        mpVariables[i] = bdd_ithvar(i*2);
        mpNextVariables[i] = bdd_ithvar(i*2 + 1);
    }
    for (unsigned node_id=0; node_id < mNodesNumber;node_id++){
        createBooleanFunction(node_id, RandomNumberGenerator::Instance()->ranf() <= probability_canalyzing_function);
    }
}

RandomBooleanNetwork::RandomBooleanNetwork(const std::string file_path,
        double probability_canalyzing_function)
{
    if (probability_canalyzing_function < 0.0 || probability_canalyzing_function > 1.0)
        EXCEPTION("Range of a probability must be between 0 and 1.");
    if (file_path.size() > 4)
    {
        try
        {
            if (file_path.compare(file_path.size()-4,4,".gml") == 0)
            {
                createGraphFromGmlFile(file_path);
            } else EXCEPTION("File format is not correct.");

            initBinaryDecisionDiagram();
            mpReverseTransitionFunction = new bdd();
            *mpReverseTransitionFunction = bddtrue;
            mpNodeFunction = new bdd[mNodesNumber];
            mpNodeNthFunction = new bdd[mNodesNumber];
            mpVariables = new bdd[mNodesNumber];
            mpNextVariables = new bdd[mNodesNumber];

            for (unsigned i=0; i<mNodesNumber; i++)
            {
                mpVariables[i] = bdd_ithvar(i*2);
                mpNextVariables[i] = bdd_ithvar(i*2 + 1);
            }

            for (unsigned node_id=0; node_id < mNodesNumber;node_id++){
                createBooleanFunction(node_id, RandomNumberGenerator::Instance()->ranf() <= probability_canalyzing_function);
            }

        } catch (Exception& e)
        {
            EXCEPTION("Error reading the file.");
        }
    } else EXCEPTION("Error in the file path.");
}

RandomBooleanNetwork::RandomBooleanNetwork(const std::string file_path)
{
    if (file_path.size() > 4)
    {
        try {
            if (file_path.compare(file_path.size()-4,4,".net") == 0
                    ||  file_path.compare(file_path.size()-4,4,"cnet") == 0)
            {
                createNetworkFromNetFile(file_path);
            }
            else EXCEPTION("File format is not correct.");
        } catch (Exception& e)
        {
            EXCEPTION("Error reading the file.");
        }
    } else EXCEPTION("Error in the file path.");
}

RandomBooleanNetwork::~RandomBooleanNetwork()
{
    delete mpReverseTransitionFunction;
    delete[] mpVariables;
    delete[] mpNextVariables;
    delete[] mpNodeFunction;
    delete[] mpNodeNthFunction;
    delete mpRbnGraph;

    mAttractors.clear();
    mpReverseTransitionFunction = NULL;
    mpVariables = NULL;
    mpNextVariables = NULL;
    mpNodeFunction = NULL;
    mpNodeNthFunction = NULL;
    mpRbnGraph = NULL;
    //mAttractors = NULL;

    //bdd_done();
    RandomNumberGenerator::Instance()->Destroy();

}

void RandomBooleanNetwork::createNetworkFromNetFile(const std::string file_path)
{
    std::string line;
    std::ifstream input_file (file_path.c_str());
    unsigned vertices_number = 0;
    std::vector<std::string> strs;
    bool vertices_line_flag = false;

    try
    {
        if (!input_file.is_open()) EXCEPTION("Not able to open the file.");
        do
        {
            if (!std::getline(input_file,line))
                EXCEPTION("Error vertices line not found.");
            boost::algorithm::trim(line);
            if (line.empty()) continue;
            boost::algorithm::to_lower(line);
            boost::split(strs, line, boost::is_any_of("\t "), boost::algorithm::token_compress_on);
            if (strs.size() != 0)
            {
                if (strs.at(0).at(0) == '#') continue;
                else if (strs.size() != 2) EXCEPTION("Error vertices line not found.");
                else
                {
                    if (strs.at(0).compare(".v") == 0)
                    {
                        vertices_number = boost::lexical_cast<unsigned>(strs.at(1));
                        if (vertices_number > 0) vertices_line_flag = true;
                        else EXCEPTION("Number of vertices not valid.");
                    }

                }
            }
        } while (!vertices_line_flag);
    } catch (const Exception&)
    {
        input_file.close();
        throw;
    }

    /* Alloc memory */
    mNodesNumber = vertices_number;
    mpRbnGraph = new ArrayDirectedGraph(mNodesNumber);
    initBinaryDecisionDiagram();
    mpReverseTransitionFunction = new bdd();
    *mpReverseTransitionFunction = bddtrue;

    mpNodeFunction = new bdd[mNodesNumber];
    mpNodeNthFunction = new bdd[mNodesNumber];
    mpVariables = new bdd[mNodesNumber];
    mpNextVariables = new bdd[mNodesNumber];

    for (unsigned i=0; i<mNodesNumber; i++)
    {
        mpVariables[i] = bdd_ithvar(i*2);
        mpNextVariables[i] = bdd_ithvar(i*2 + 1);
    }
    /* End Alloc */

    try
    {
        unsigned counter_node_id = 0, id = 0;
        while (counter_node_id < mNodesNumber)
        {
            bool node_line_flag = false;
            do
            {
                if (!std::getline(input_file,line))
                    EXCEPTION("The input file is not structured in the correct way.");
                boost::algorithm::trim(line);
                if (line.empty()) continue;
                boost::algorithm::to_lower(line);
                boost::split(strs, line, boost::is_any_of("\t "), boost::algorithm::token_compress_on);
                if (strs.size() != 0)
                {
                    if (strs.at(0).at(0) == '#') continue;
                    else if (strs.size() < 3) EXCEPTION("The input file is not structured in the correct way.");
                    else
                    {
                        if (strs.at(0).compare(".n") == 0)
                        {
                            id = boost::lexical_cast<unsigned>(strs.at(1));
                            id--;
                            if (id == counter_node_id) node_line_flag = true;
                            else EXCEPTION("Number of vertices not valid.");
                        }
                    }
                }
            } while (!node_line_flag);
            counter_node_id++;
            unsigned number_of_incoming_vertices = boost::lexical_cast<unsigned>(strs.at(2));
            if (number_of_incoming_vertices + 3 != strs.size())
                EXCEPTION("Error reading node " + boost::lexical_cast<std::string>(counter_node_id) + " in node declaration.");
            bdd function = bddfalse;
            bdd reverse_function = bddfalse;

            if (number_of_incoming_vertices == 0)
            {
                do
                {
                    if (!std::getline(input_file,line))
                    {
                        if (counter_node_id == mNodesNumber) node_line_flag = false;
                        else EXCEPTION("Error reading node " + boost::lexical_cast<std::string>(counter_node_id) + ". File ended.");

                    }
                    boost::algorithm::trim(line);
                    if (line.empty())
                    {
                        node_line_flag = false;
                        continue;
                    }
                    boost::algorithm::to_lower(line);
                    boost::split(strs, line, boost::is_any_of("\t "), boost::algorithm::token_compress_on);
                    if (strs.size() != 0)
                    {
                        if (strs.at(0).at(0) == '#') continue;
                        else if (strs.at(0).compare(".n") == 0) EXCEPTION ("Error reading node " + boost::lexical_cast<std::string>(counter_node_id) +
                                ". Found a node line before a empty line.");
                        else
                        {
                            if (strs.at(0).compare("0") == 0) continue;
                            else if (strs.at(0).compare("1") == 0)
                            {
                                function |= bddtrue;
                                reverse_function |= bddtrue;
                            } else EXCEPTION ("Error reading node " + boost::lexical_cast<std::string>(counter_node_id) +
                                    ". Unknown character found.");
                        }
                    } else node_line_flag = false;
                } while (node_line_flag);
            }
            else
            {
                unsigned j = 0;
                std::vector<unsigned> input_ids;
                while (j < number_of_incoming_vertices)
                {
                    unsigned id_input = boost::lexical_cast<unsigned>(strs.at(j+3));
                    id_input--;
                    if (id_input > mNodesNumber)
                        EXCEPTION ("Error reading node " + boost::lexical_cast<std::string>(counter_node_id) +
                                    ". in node declaration. Input node id value is greater or" +
                                    " equal to the number of vertices.");
                    if (mpRbnGraph->addEdgeById(id_input,id))
                        input_ids.push_back(id_input);
                    else EXCEPTION ("Error reading node " + boost::lexical_cast<std::string>(counter_node_id) +
                            ". in node declaration. Input node id value is not correct.");
                    j++;
                }
                bool eof = false;
                bool empty_line = false;
                do
                {
                    if (std::getline(input_file,line))
                    {
                        boost::algorithm::trim(line);
                        if (line.empty())
                        {
                            empty_line = true;
                            continue;
                        }
                        boost::algorithm::to_lower(line);
                        boost::split(strs, line, boost::is_any_of("\t "), boost::algorithm::token_compress_on);
                        if (strs.size() != 0)
                        {
                            if (strs.at(0).at(0) == '#') continue;
                            else if (strs.size() != 2 || strs.at(0).size() != number_of_incoming_vertices)
                                EXCEPTION ("Error reading node " + boost::lexical_cast<std::string>(counter_node_id) + " in the function definition.");
                            else
                            {
                                //we don't care about false values, because the function in an OR of TRUE inputs.
                                if (strs.at(1).compare("0") == 0) continue;
                                else if (strs.at(1).compare("1") == 0)
                                {
                                    j = 0;
                                    bdd inner_function = bddtrue;
                                    bdd reverse_inner_function = bddtrue;
                                    while (j < number_of_incoming_vertices)
                                    {
                                        char input_value = strs.at(0).at(j);
                                        switch (input_value)
                                        {
                                            case '0':
                                                //inner_function = inner_function AND (NOT x_j)
                                                inner_function &= !mpVariables[input_ids.at(j)];
                                                reverse_inner_function &= !mpNextVariables[input_ids.at(j)];
                                                break;
                                            case '1':
                                                //inner_function = inner_function AND (x_j)
                                                inner_function &= mpVariables[input_ids.at(j)];
                                                reverse_inner_function &= mpNextVariables[input_ids.at(j)];
                                                break;
                                            case '-': break; //x_j is free for this input.
                                            default:
                                                EXCEPTION ("Error reading node " + boost::lexical_cast<std::string>(counter_node_id) + " in the function definition:"
                                                        + " unknown symbol.");
                                        }
                                        j++;
                                    }
                                    function |= inner_function;
                                    reverse_function |= reverse_inner_function;
                                } else EXCEPTION ("Error reading node " + boost::lexical_cast<std::string>(counter_node_id) + " in the function definition.");
                            }
                        }
                    } else eof = true;
                } while (!empty_line && !eof);
                if (eof && counter_node_id < mNodesNumber) EXCEPTION ("Error reading the file. Some vertices is missing.");
            }
            *mpReverseTransitionFunction &= bdd_apply(mpVariables[id], reverse_function, bddop_biimp);
            mpNodeFunction[id] = function;
            mpNodeNthFunction[id] = function;
        }
    } catch(const Exception& e)
    {
        delete mpReverseTransitionFunction;
        delete[] mpVariables;
        delete[] mpNextVariables;
        delete[] mpNodeFunction;
        delete[] mpNodeNthFunction;
        delete mpRbnGraph;

        mAttractors.clear();
        mpReverseTransitionFunction = NULL;
        mpVariables = NULL;
        mpNextVariables = NULL;
        mpNodeFunction = NULL;
        mpNodeNthFunction = NULL;
        mpRbnGraph = NULL;

        bdd_done();
        input_file.close();
        throw;
    }
}

void RandomBooleanNetwork::createGraphFromGmlFile(const std::string file_path)
{
    std::string line;
    std::ifstream input_file (file_path.c_str());
    bool graph_found = false;
    bool node_found = false;
    bool edge_found = false;
    bool source_found = false;
    bool stop = false;
    unsigned vertices_number = 0;
    std::vector<std::string> strs;
    std::map<std::string, unsigned> node_map;
    try
    {
        if (!input_file.is_open()) EXCEPTION("Not able to open the file."); //EXCEPTION
        do
        {
            /*DIRECTED 1 */
            if (!std::getline(input_file,line))
                EXCEPTION("The input file is not structured in the correct way.");
            boost::algorithm::trim(line);
            boost::algorithm::to_lower(line);
            boost::split(strs, line, boost::is_any_of("\t "), boost::algorithm::token_compress_on);
            if (strs.size() != 0)
            {
                if (strs.size() > 2) EXCEPTION("The input file is not structured in the correct way.");
                else
                {
                    if (strs.at(0).compare("graph") == 0)
                    {
                        graph_found = true;
                    }
                    else EXCEPTION("The input file is not well-structured.");
                }
            }
        }
        while (!graph_found);

        do
        {
            if (!std::getline(input_file,line))
            {
                stop = true;
                continue;
            }
            boost::algorithm::trim(line);
            boost::algorithm::to_lower(line);
            boost::split(strs, line, boost::is_any_of("\t "), boost::algorithm::token_compress_on);
            if (strs.size() != 0)
            {
                if (strs.size() < 3)
                {
                    if (strs.at(0).compare("node") == 0)
                    {
                        node_found = true;
                    }
                    else if (strs.at(0).compare("edge") == 0)
                    {
                        edge_found = true;
                    }
                    else if (node_found && strs.at(0).compare("id")==0)
                    {
                        if (node_map.find(strs.at(1)) == node_map.end())
                        {
                            node_map.insert(std::pair<std::string, unsigned>(strs.at(1), vertices_number));
                            vertices_number++;
                            node_found = false;
                        }
                        else EXCEPTION("The input file is not well-structured.");
                    }
                }
            }
        } while (!edge_found && !stop);
    } catch (Exception& e)
    {
        input_file.close();
        throw;
        //TERMINATE(e.GetMessage());
    }
    if (vertices_number > 0)
    {
        mNodesNumber = vertices_number;
        mpRbnGraph = new ArrayDirectedGraph(mNodesNumber);
        if (!stop)
        {
            try
            {
                unsigned source;
                unsigned target;
                do
                {
                    boost::algorithm::trim(line);
                    boost::algorithm::to_lower(line);
                    boost::split(strs, line, boost::is_any_of("\t "), boost::algorithm::token_compress_on);
                    if (strs.size() != 0)
                    {
                        if (strs.size() < 3)
                        {
                            if (strs.at(0).compare("edge") == 0)
                            {
                                edge_found = true;
                            }
                            else if (edge_found && strs.size() == 2 && strs.at(0).compare("source") == 0)
                            {
                                source_found = true;
                                std::map<std::string, unsigned>::iterator map_it;
                                map_it = node_map.find(strs.at(1));
                                if (map_it != node_map.end())
                                {
                                    source = map_it->second;
                                }
                                else EXCEPTION("The input file is not well-structured.");
                            }
                            else if (source_found && strs.size() == 2 && strs.at(0).compare("target")==0)
                            {
                                std::map<std::string, unsigned>::iterator map_it;
                                map_it = node_map.find(strs.at(1));
                                if (map_it != node_map.end())
                                {
                                    target = map_it->second;
                                    mpRbnGraph->addEdgeById(source, target);
                                    edge_found = false;
                                    source_found = false;
                                }
                                else EXCEPTION("The input file is not well-structured.");
                            }
                        }
                    }
                } while (std::getline(input_file,line));
                if (edge_found || source_found) EXCEPTION("The input file is not well-structured.");
            } catch (Exception& e)
            {
                delete mpRbnGraph;
                input_file.close();
                throw;
                //TERMINATE(e.GetMessage());
            }
        }
    } else EXCEPTION("The input file is not well-structured.");
}

void RandomBooleanNetwork::initBinaryDecisionDiagram() const
{
    //Graph may be ordered by number of input vertices calling graph.sortGraph().
    /*
    unsigned max_input_number = mpRbnGraph->getIncomingVerticesNumberById(0);
    unsigned cache_size = 0;

    if (max_input_number < 50) {
        if (mNodesNumber <= 5)
        {
            cache_size = 1000;
        } else if (mNodesNumber <= 10)
        {
            cache_size = 2000;
        } else if (mNodesNumber <= 15)
        {
            cache_size = 6000;
        } else if (mNodesNumber <= 20)
        {
            cache_size = 20000;
        } else if (mNodesNumber <= 25)
        {
            cache_size = 40000;
        } else
        {
            cache_size = 200000;
        }
        if (mNodesNumber <= 30)
        {
            cache_size = 200000;
        } else if (mNodesNumber <= 35)
        {
            cache_size = 1000000;
        } else
        {
            cache_size = 2000000;
        }
    }
    bdd_init(cache_size*10, cache_size);
    */
    bdd_setvarnum(mNodesNumber * 2);
    bdd_setcacheratio(64);
    //bdd_setmaxincrease(node_number);
    bdd_varblockall();
}

void RandomBooleanNetwork::createBooleanFunction(unsigned node_id, bool canalyzing_function)
{
    assert(node_id < mNodesNumber);
    bdd function = bddfalse;
    bdd reverse_function = bddfalse;
    unsigned inputs_number = mpRbnGraph->getIncomingVerticesNumberById(node_id);

    if (inputs_number < 2)
    {
        if (inputs_number == 1)
        {
            std::vector<unsigned> input_vertices = mpRbnGraph->getIncomingVerticesById(node_id);
            bool is_activator = RandomNumberGenerator::Instance()->ranf() > 0.5;
            if (is_activator)
            {
                function |= mpVariables[input_vertices[0]];
                reverse_function |= mpNextVariables[input_vertices[0]];
            }
            else
            {
                function |= !mpVariables[input_vertices[0]];
                reverse_function |= !mpNextVariables[input_vertices[0]];
            }
        }
        else
        {
            mpRbnGraph->addEdgeById(node_id, node_id);
            function |= mpVariables[node_id];
            reverse_function |= mpNextVariables[node_id];
        }
    }
    else
    {
        unsigned* input_vertices = mpRbnGraph->getIncomingVerticesArrayIdById(node_id);
        for (unsigned end=inputs_number-1; end>0; end--)
        {
          // Pick a random integer from {0,..,end}
          unsigned k = RandomNumberGenerator::Instance()->randMod(end+1);
          unsigned temp = input_vertices[end];
          input_vertices[end] = input_vertices[k];
          input_vertices[k] = temp;
        }
        mpRbnGraph->setIncomingVerticesById(input_vertices, inputs_number, node_id);
        if (canalyzing_function)
        {
            //Nested canalyzing function (NCF)
            for (unsigned i=0; i<inputs_number; i++)
            {
                bool canalyzed_output = RandomNumberGenerator::Instance()->ranf() > 0.5;
                if (canalyzed_output) //Output == 1
                {
                    bool canalyzed_input = RandomNumberGenerator::Instance()->ranf() > 0.5;
                    if (canalyzed_input) //Input == 1
                    {
                        function |= mpVariables[i];
                        reverse_function |= mpNextVariables[i];
                    }
                    else //Input == 0
                    {
                        function |= !mpVariables[i];
                        reverse_function |= !mpNextVariables[i];
                    }
                }
            }
        }
        else
        {
            unsigned inputs_count = 0;
            bool is_activator_dominated = RandomNumberGenerator::Instance()->ranf() > 0.5;
            bool has_activator = false;
            bool has_inhibitor = false;
            bdd activator_functions = bddfalse;
            bdd reverse_activator_functions = bddfalse;
            bdd inhibitor_functions = bddfalse;
            bdd reverse_inhibitor_functions = bddfalse;
            while (inputs_count < inputs_number)
            {
                bdd inner_function = bddtrue;
                bdd reverse_inner_function = bddtrue;
                bool activator_function = RandomNumberGenerator::Instance()->ranf() > 0.5;
                unsigned number_of_function_inputs = RandomNumberGenerator::Instance()->randMod(inputs_number - inputs_count) + 1;
                unsigned number_of_activators = RandomNumberGenerator::Instance()->randMod(number_of_function_inputs + 1);
                for (unsigned i=0; i<number_of_activators; i++)
                {
                    inner_function &= mpVariables[input_vertices[i + inputs_count]];
                    reverse_inner_function &= mpNextVariables[input_vertices[i + inputs_count]];
                }
                for (unsigned i=0; i<(number_of_function_inputs - number_of_activators); i++)
                {
                    inner_function &= !mpVariables[input_vertices[i + inputs_count + number_of_activators]];
                    reverse_inner_function &= !mpNextVariables[input_vertices[i + inputs_count + number_of_activators]];
                }
                if (activator_function)
                {
                    has_activator = true;
                    activator_functions |= inner_function;
                    reverse_activator_functions |= reverse_inner_function;
                }
                else
                {
                    has_inhibitor = true;
                    inhibitor_functions |= inner_function;
                    reverse_inhibitor_functions |= reverse_inner_function;
                }
                inputs_count += number_of_function_inputs;
            }
            if (has_inhibitor && has_activator) {
                if (is_activator_dominated)
                {
                    function |= activator_functions | (!inhibitor_functions);
                    reverse_function |= reverse_activator_functions | (!reverse_inhibitor_functions);
                }
                else
                {
                    function |= activator_functions & (!inhibitor_functions);
                    reverse_function |= reverse_activator_functions & (!reverse_inhibitor_functions);
                }
            }
            else
            {
                if (has_activator) {
                    function |= activator_functions;
                    reverse_function |= reverse_activator_functions;
                }
                if (has_inhibitor) {
                    function |= !inhibitor_functions;
                    reverse_function |= !reverse_inhibitor_functions;
                }
            }
        }
        delete[] input_vertices;
        input_vertices = NULL;
    }
    *mpReverseTransitionFunction &= bdd_apply(mpVariables[node_id], reverse_function, bddop_biimp);
    mpNodeFunction[node_id] = function;
    mpNodeNthFunction[node_id] = function;
    bdd_reorder(2);

}

bdd RandomBooleanNetwork::find_next_cycles(bddPair* variables_pair){
    assert(mNodesNumber > 0);
    //This method compute delta_j+1, fixed delta_j. Transition function
    //T^j is memorized in mpNodeNthFunction[] and will be upgrated to T^j+1.
    for(unsigned i = 0; i < mNodesNumber; i++){
        mpNodeNthFunction[i] = bdd_veccompose(mpNodeNthFunction[i], variables_pair);
    }

    bdd new_states = bddtrue;
    for (unsigned i = 0; i < mNodesNumber; i++) {
        new_states &= bdd_apply(mpVariables[i], mpNodeNthFunction[i], bddop_biimp);
    }
    return new_states;

}

bdd RandomBooleanNetwork::find_backward_reachable_states(bdd current_ring, bdd states,
        bdd set_variables, bddPair* variables_pair, unsigned &steps_max)
{
    //This method will find the backward reachable states using the Reverse Transition Function
    //if a given bdd 'states' which can be a single or a set of states.
    assert(mpReverseTransitionFunction);

    unsigned length = 0;
    bool ring_found = false;
    do
    {
        states = bdd_appex(states, *mpReverseTransitionFunction, bddop_and, set_variables);
        states = bdd_replace(states, variables_pair);
        if (current_ring == (current_ring | states)) {
            ring_found = true;
        }
        else
        {
            length++;
            current_ring |= states;
        }
    } while (states != bddfalse && !ring_found && current_ring != bddtrue);
    if (length > steps_max)
    {
        steps_max = length;
    }
    return current_ring;
}

void RandomBooleanNetwork::findAttractors()
{
    assert(mNodesNumber > 0);

    if (mAttractors.empty()) {
        mAttractorLength.clear();
        unsigned j = 1;
        unsigned steps_max = 1; //such as forall x, F^step_max (x) = y where y is a state of an attractor.
        bddPair* replace_backward_assignment;
        bddPair* replace_forward_assignment;
        replace_backward_assignment = bdd_newpair();
        replace_forward_assignment = bdd_newpair();

        bdd explored_states = bddfalse; //totalR
        bdd states_return_to_themself = bddtrue; //I

        int* p_variables_id = new int[mNodesNumber];
        int* p_next_variables_id = new int[mNodesNumber];

        int index = 0;
        for(unsigned i = 0; i < mNodesNumber; i++){
            states_return_to_themself &= bdd_apply(mpVariables[i], mpNodeNthFunction[i], bddop_biimp);
            p_variables_id[i] = index*2;
            p_next_variables_id[i] = index*2 +1;
            index++;
        }
        bdd_setpairs(replace_backward_assignment, p_next_variables_id, p_variables_id, mNodesNumber);
        bdd_setbddpairs(replace_forward_assignment, p_variables_id, mpNodeNthFunction, mNodesNumber);
        bdd set_variables = bdd_makeset(p_variables_id, mNodesNumber);

        delete[] p_variables_id;
        delete[] p_next_variables_id;
        p_variables_id = NULL;
        p_next_variables_id = NULL;

        while(explored_states != bddtrue)
        {
            if(j > 1)
            {
                states_return_to_themself = find_next_cycles(replace_forward_assignment);
                if (states_return_to_themself == bddtrue)
                {
                    storageAttractors((!explored_states), j, set_variables, replace_backward_assignment);
                    break;
                }

                states_return_to_themself = states_return_to_themself - explored_states;
                if(states_return_to_themself == bddfalse){
                    j++;
                    continue;
                }
            }
            j++;
            bdd current_ring = bddfalse;
            storageAttractors(states_return_to_themself, j-1, set_variables, replace_backward_assignment);
            current_ring |= states_return_to_themself;
            explored_states |= find_backward_reachable_states(current_ring,
                            states_return_to_themself,set_variables, replace_backward_assignment, steps_max);
        }
        if (steps_max > j-1)
        {
            steps_max -= j-1;
            do {
                for(unsigned i = 0; i < mNodesNumber; i++){
                    mpNodeNthFunction[i] = bdd_veccompose(mpNodeNthFunction[i], replace_forward_assignment);
                }
                steps_max --;
            } while (steps_max > 0);
        }
    } // ELSE ALREADY FOUNDED!
    /*
    for (unsigned i=0; i<mAttractors.size(); i++)
    {
        std::cout<< i << ". ";
        bdd_allsat(mAttractors.at(i), allsatHandlerPrint);
    }
    */
}

void RandomBooleanNetwork::storageAttractors(bdd states_return_to_themself, unsigned j,
        bdd set_variables, bddPair* variables_pair)
{
    while (states_return_to_themself != bddfalse)
    {
        bdd attractor_states = bdd_satoneset(states_return_to_themself, set_variables, bddtrue);
        bdd backward_states = attractor_states;
        do
        {
            states_return_to_themself -= backward_states;
            backward_states = bdd_appex(backward_states, *mpReverseTransitionFunction, bddop_and, set_variables);
            backward_states = bdd_replace(backward_states, variables_pair);
            backward_states &= states_return_to_themself;
            if (backward_states == bddfalse)
            {
                break;
            }
            attractor_states |= backward_states;
        }
        while (states_return_to_themself != (states_return_to_themself - backward_states));
        mAttractors.push_back(attractor_states);
        mAttractorLength.push_back(j);
    }
    //std::cout<<"Attractors of size " << j << "\n";
}

std::vector<std::map<unsigned,double> > RandomBooleanNetwork::getAttractorMatrix() const
{
    assert(!mAttractors.empty());

    bdd transition_function = bddtrue;
    bddPair* replace_forward_assignment = bdd_newpair();

    std::vector<std::map<unsigned,unsigned> > frequency_attractor_matrix;

    int* p_variables_id = new int[mNodesNumber];
    int* p_next_variables_id = new int[mNodesNumber];

    for (unsigned i=0; i < mNodesNumber; i++)
    {
        transition_function &= bdd_apply(mpNextVariables[i], mpNodeNthFunction[i], bddop_biimp);
        p_variables_id[i] = (int)i*2;
        p_next_variables_id[i] = (int) i*2 +1;
    }

    bdd set_variables = bdd_makeset(p_variables_id, mNodesNumber);
    bdd_setpairs(replace_forward_assignment, p_next_variables_id, p_variables_id, mNodesNumber);
    delete[] p_variables_id;
    delete[] p_next_variables_id;

    for (unsigned i=0; i < mAttractors.size(); i++)
    {
        std::map<unsigned,unsigned> map;
        frequency_attractor_matrix.push_back(map);
        std::map<unsigned,unsigned>::iterator iterator;
        bdd current_attractor = mAttractors.at(i);
        do
        {
            bdd state = bdd_satoneset(current_attractor, set_variables, bddtrue);
            current_attractor -= state;
            for (unsigned j=0; j<mNodesNumber; j++)
            {
                bdd flip_state = bdd_compose(state, bdd_nithvar(j*2), j*2);
                unsigned position = getStateAttractor(flip_state, transition_function, set_variables, replace_forward_assignment);

                iterator = frequency_attractor_matrix.at(i).find(position);
                if (iterator != frequency_attractor_matrix.at(i).end())
                {
                    iterator->second++;
                }
                else
                {
                    frequency_attractor_matrix.at(i).insert(std::pair<unsigned,unsigned>(position,1));
                }
            }
        } while (current_attractor != bddfalse);
    }

    std::vector<std::map<unsigned,double> > stochastic_matrix = getStochasticMatrix(frequency_attractor_matrix);

    return stochastic_matrix;
}

unsigned RandomBooleanNetwork::getStateAttractor(
            bdd flip_state, bdd transition_function,
            bdd set_variables, bddPair* replace_forward_assignment) const
{
    assert(!mAttractors.empty());

    flip_state = bdd_appex(flip_state, transition_function, bddop_and, set_variables);
    flip_state = bdd_replace(flip_state, replace_forward_assignment);
    bool found = false;
    unsigned index = 0;
    do
    {
        if ((flip_state & mAttractors.at(index)) != bddfalse)
        {
            found = true;
        }
        index++;
    } while (!found && index < mAttractors.size());
    return (index - 1);
}

std::vector<std::map<unsigned,double> > RandomBooleanNetwork::getStochasticMatrix(
        std::vector<std::map<unsigned,unsigned> > frequency_matrix) const
{
    assert(!mAttractors.empty());

    std::vector<std::map<unsigned,double> >  stochastic_matrix;
    for (unsigned row=0; row< mAttractors.size(); row++)
    {
        std::map<unsigned,double> map;
        stochastic_matrix.push_back(map);
        unsigned row_sum = 0;
        std::map<unsigned,unsigned>::iterator iterator;
        for (iterator=frequency_matrix.at(row).begin(); iterator!=frequency_matrix.at(row).end(); ++iterator)
        {
            row_sum += iterator->second;
        }
        for (iterator=frequency_matrix.at(row).begin(); iterator!=frequency_matrix.at(row).end(); ++iterator)
        {
            stochastic_matrix.at(row).insert(std::pair<unsigned,double>(iterator->first,
                    ((double) iterator->second / (double) row_sum)));
        }
    }
    return stochastic_matrix;
}

void RandomBooleanNetwork::printNetworkToNetFile(const std::string directory,
		const std::string filename) const
{
    assert(mpRbnGraph);

    if (directory.empty())
        EXCEPTION("Directory name not valid.");

    if (filename.size() < 5 || filename.compare(filename.size()-4,4,".net") != 0)
        EXCEPTION("File path not valid. It must terminate with '.net' extension.");

    OutputFileHandler handler(directory,false);
	out_stream p_file = handler.OpenOutputFile(filename);

	*p_file << "#This file has been generated by Chaste Random Boolean Network tool\n";
	*p_file << ".v " << mNodesNumber << "\n\n";

	if (mAttractors.size()>0)
	{
		*p_file << "# As a result of simulation, we get the following "
				<< mAttractors.size() << " attractors:\n";
		for (unsigned i=0; i<mAttractors.size(); i++)
		{
			stringBuffer.clear();
			bdd_allsat(mAttractors.at(i), allsatHandlerPrintToStringBuffer);
			*p_file << "#\n";
			for (unsigned j=0; j<stringBuffer.size(); j++)
			{
				*p_file << "# " << stringBuffer.at(j) << "\n";
			}
			*p_file << "# Attractor " << i + 1 << " is of length " << mAttractorLength.at(i) << "\n";
		}
		*p_file << "\n";
	}


	for (unsigned i=0; i < mNodesNumber; i++)
	{
		unsigned inputs_number = mpRbnGraph->getIncomingVerticesNumberById(i);
		*p_file << ".n " << i+1 << " " << inputs_number;
		std::vector<unsigned> incoming_nodes = mpRbnGraph->getIncomingVerticesById(i);
		for (unsigned j=0; j < inputs_number; j++)
		{
			*p_file << " " << incoming_nodes[j] + 1;
		}
		*p_file << "\n";

		if (inputs_number == 0)
		{
			if (mpNodeFunction[i] == bddtrue) *p_file << "1\n";
			else *p_file << "0\n";
		}
		else
		{
			stringBuffer.clear();
			bdd_allsat(mpNodeFunction[i], allsatHandlerPrintToStringBuffer);
			for (unsigned j=0; j < stringBuffer.size(); j++)
			{
				std::string inputs = stringBuffer.at(j);
				for (unsigned k=0; k<inputs_number; k++)
				{
					*p_file << inputs.at(incoming_nodes[k]);
				}
				*p_file << " 1\n";
			}
		}
		*p_file << std::endl;
	}
	p_file->close();

}

void RandomBooleanNetwork::printNetworkToCnetFile(const std::string directory,
		const std::string filename) const
{
    assert(mpRbnGraph);

    if (directory.empty())
        EXCEPTION("Directory name not valid.");

    if (filename.size() < 6 || filename.compare(filename.size()-5,5,".cnet") != 0)
        EXCEPTION("File path not valid. It must terminate with '.cnet' extension.");

    OutputFileHandler handler(directory,false);
	out_stream p_file = handler.OpenOutputFile(filename);

	*p_file << "#This file has been generated by Chaste Random Boolean Network tool\n";
	*p_file << ".v " << mNodesNumber << "\n\n";

	if (!mAttractors.empty())
	{
		*p_file << "# As a result of simulation, we get the following "
				<< mAttractors.size() << " attractors:\n";
		for (unsigned i=0; i<mAttractors.size(); i++)
		{
			stringBuffer.clear();
			bdd_allsat(mAttractors.at(i), allsatHandlerPrintToStringBuffer);
			*p_file << "#\n";
			for (unsigned j=0; j<stringBuffer.size(); j++)
			{
				*p_file << "# " << stringBuffer.at(j) << "\n";
			}
			*p_file << "# Attractor " << i + 1 << " is of length " << mAttractorLength.at(i) << "\n";
		}
		*p_file << "\n";
	}


	for (unsigned i=0; i < mNodesNumber; i++)
	{
		unsigned inputs_number = mpRbnGraph->getIncomingVerticesNumberById(i);
		*p_file << ".n " << i+1 << " " << inputs_number;
		std::vector<unsigned> incoming_nodes = mpRbnGraph->getIncomingVerticesById(i);
		for (unsigned j=0; j < inputs_number; j++)
		{
			*p_file << " " << incoming_nodes[j] + 1;
		}
		*p_file << "\n";

		if (inputs_number == 0)
		{
			if (mpNodeFunction[i] == bddtrue) *p_file << "1\n";
			else *p_file << "0\n";
		}
		else
		{
			stringBuffer.clear();
			bdd_allsat(mpNodeFunction[i], allsatHandlerPrintToStringBuffer);
			for (unsigned j=0; j < stringBuffer.size(); j++)
			{
				std::string inputs = stringBuffer.at(j);
				for (unsigned k=0; k<inputs_number; k++)
				{
					*p_file << inputs.at(incoming_nodes[k]);
				}
				*p_file << " 1\n";
			}
			stringBuffer.clear();
			bdd_allsat(bddtrue - mpNodeFunction[i], allsatHandlerPrintToStringBuffer);
			for (unsigned j=0; j < stringBuffer.size(); j++)
			{
				std::string inputs = stringBuffer.at(j);
				for (unsigned k=0; k<inputs_number; k++)
				{
					*p_file << inputs.at(incoming_nodes[k]);
				}
				*p_file << " 0\n";
			}
		}
		*p_file << std::endl;
	}
	p_file->close();
}

void RandomBooleanNetwork::printNetworkToBoolNetFile(const std::string directory,
		const std::string filename) const
{
    assert(mpRbnGraph);

    if (directory.empty())
        EXCEPTION("Directory name not valid.");

    if (filename.size() < 5 || filename.compare(filename.size()-4,4,".txt") != 0)
        EXCEPTION("File path not valid. It must terminate with '.txt' extension.");

    OutputFileHandler handler(directory,false);
	out_stream p_file = handler.OpenOutputFile(filename);

	*p_file << "targets, factors\n";
	for (unsigned i=0; i < mNodesNumber; i++)
	{
		unsigned inputs_number = mpRbnGraph->getIncomingVerticesNumberById(i);
		*p_file << "Gene" << i+1 << ", ";
		std::vector<unsigned> incoming_nodes = mpRbnGraph->getIncomingVerticesById(i);

		if (mpNodeFunction[i] == bddtrue) *p_file << "! Gene" << i+1 << " | Gene" << i+1;
		else if (mpNodeFunction[i] == bddfalse) *p_file << "! Gene" << i+1 << " & Gene" << i+1;
		else
		{
			stringBuffer.clear();
			bdd_allsat(mpNodeFunction[i], allsatHandlerPrintToStringBuffer);
			for (unsigned j=0; j < stringBuffer.size(); j++)
			{
				if (stringBuffer.size()!= 1)
				{
					if (j==0) *p_file << "(";
					else *p_file << " | (";
				}
				std::string inputs = stringBuffer.at(j);
				bool inserted = false;
				for (unsigned k=0; k<inputs_number; k++)
				{
					if (inserted)
					{
						if (inputs.at(incoming_nodes[k]) == '1')
							*p_file << " & Gene"<< incoming_nodes[k] + 1;
						else if (inputs.at(incoming_nodes[k]) == '0')
							*p_file << " & ! Gene"<< incoming_nodes[k] + 1;
					}
					else
					{
						if (inputs.at(incoming_nodes[k]) == '1')
						{
							*p_file << "Gene" << incoming_nodes[k] + 1;
							inserted = true;
						}
						else if (inputs.at(incoming_nodes[k]) == '0')
						{
							*p_file << "! Gene"<< incoming_nodes[k] + 1;
							inserted = true;
						}
					}

				}
				if (stringBuffer.size()!= 1) *p_file << ") ";
			}
		}
		*p_file << std::endl;
	}
	*p_file << std::endl;
	p_file->close();
}

void RandomBooleanNetwork::printNetworkToBooleanNetFile(const std::string directory,
		const std::string filename) const
{
    assert(mpRbnGraph);

    if (directory.empty())
        EXCEPTION("Directory name not valid.");

    if (filename.size() < 5 || filename.compare(filename.size()-4,4,".txt") != 0)
        EXCEPTION("File path not valid. It must terminate with '.txt' extension.");

    OutputFileHandler handler(directory,false);
	out_stream p_file = handler.OpenOutputFile(filename);

	*p_file << "#This file has been generated by Chaste Random Boolean Network tool\n\n";
	for (unsigned i=0; i < mNodesNumber; i++)
	{
		unsigned inputs_number = mpRbnGraph->getIncomingVerticesNumberById(i);
		*p_file << "1: Gene" << i+1 << "* = ";
		std::vector<unsigned> incoming_nodes = mpRbnGraph->getIncomingVerticesById(i);

		if (mpNodeFunction[i] == bddtrue) *p_file << "not Gene" << i+1 << " or Gene" << i+1 <<"\n";
		else if (mpNodeFunction[i] == bddfalse) *p_file << "not Gene" << i+1 << " and Gene" << i+1 <<"\n";
		else
		{
			stringBuffer.clear();
			bdd_allsat(mpNodeFunction[i], allsatHandlerPrintToStringBuffer);
			for (unsigned j=0; j < stringBuffer.size(); j++)
			{
				if (stringBuffer.size()!= 1)
				{
					if (j==0) *p_file << "(";
					else *p_file << " or (";
				}
				std::string inputs = stringBuffer.at(j);
				bool inserted = false;
				for (unsigned k=0; k<inputs_number; k++)
				{
					if (inserted)
					{
						if (inputs.at(incoming_nodes[k]) == '1')
							*p_file << " and Gene"<< incoming_nodes[k] + 1;
						else if (inputs.at(incoming_nodes[k]) == '0')
							*p_file << " and not Gene"<< incoming_nodes[k] + 1;
					}
					else
					{
						if (inputs.at(incoming_nodes[k]) == '1')
						{
							*p_file << "Gene" << incoming_nodes[k] + 1;
							inserted = true;
						}
						else if (inputs.at(incoming_nodes[k]) == '0')
						{
							*p_file << "not Gene"<< incoming_nodes[k] + 1;
							inserted = true;
						}
					}

				}
				if (stringBuffer.size()!= 1) *p_file << ") ";
			}
		}
	*p_file << std::endl;
	}
	*p_file << std::endl;
	p_file->close();
}

void RandomBooleanNetwork::printGraphToGmlFile(const std::string directory,
		const std::string filename) const
{
    assert(mpRbnGraph);

    mpRbnGraph->printGraphToGmlFile(directory, filename);
}

void RandomBooleanNetwork::printGraphToSifFile(const std::string directory,
		const std::string filename) const
{
    assert(mpRbnGraph);

    mpRbnGraph->printGraphToSifFile(directory, filename);
}

void RandomBooleanNetwork::printGraphToDotFile(const std::string directory,
		const std::string filename) const
{
    assert(mpRbnGraph);

    mpRbnGraph->printGraphToDotFile(directory, filename);
}

void RandomBooleanNetwork::printNetwork() const
{
    assert(mpRbnGraph);

    mpRbnGraph->printGraph();
    std::cout << "--------- \n";
}

unsigned RandomBooleanNetwork::getNodesNumber() const{
    return mNodesNumber;
}

unsigned RandomBooleanNetwork::getAvarageInputsPerNode() const{
    return mAverageInputsPerNode;
}

unsigned RandomBooleanNetwork::getAttractorsNumber() const{
    return mAttractors.size();
}

std::vector<unsigned> RandomBooleanNetwork::getAttractorLength() const
{
    return mAttractorLength;
}

void allsatHandlerPrint(char *varset,int size)
{
  for (int v=0; v<size; v=v+2)
  {
     std::cout << (varset[v] < 0 ? '-' : (char)('0' + varset[v]));
  }
  std::cout << std::endl;
}

void allsatHandlerPrintToStringBuffer(char *varset,int size)
{
    std::string input_values;

    for (int v=0; v<size; v=v+2)
    {
        if (varset[v] < 0) input_values += '-';
        else input_values += (char)('0' + varset[v]);
    }
    stringBuffer.push_back(input_values);
}
