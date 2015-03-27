#include "ThresholdErgodicSetDifferentiationTree.hpp"
#include <limits>

ThresholdErgodicSetDifferentiationTree::ThresholdErgodicSetDifferentiationTree(
        std::vector<std::map<unsigned,double> > stochastic_matrix,
        std::vector<unsigned> attractor_length) :
mStochasticMatrix(stochastic_matrix),
mAttractorLength(attractor_length),
mpBooleanNetwork(NULL)
{
    if (stochastic_matrix.size() != mAttractorLength.size())
        EXCEPTION("Rows (and columns) in stochastic matrix must be equal to the number of attractors");
}

ThresholdErgodicSetDifferentiationTree::ThresholdErgodicSetDifferentiationTree(unsigned nodes_number, unsigned avarage_input_number_per_node,
            bool scale_free, double probability_canalyzing_function)
{
    mpBooleanNetwork = new RandomBooleanNetwork(nodes_number,avarage_input_number_per_node,
            scale_free,probability_canalyzing_function);
    mpBooleanNetwork->findAttractors();
    mStochasticMatrix = mpBooleanNetwork->getAttractorMatrix();
    mAttractorLength = mpBooleanNetwork->getAttractorLength();
}

ThresholdErgodicSetDifferentiationTree::ThresholdErgodicSetDifferentiationTree(std::string file_path, double probability_canalyzing_function)
{
    mpBooleanNetwork = new RandomBooleanNetwork(file_path, probability_canalyzing_function);
    mpBooleanNetwork->findAttractors();
    mStochasticMatrix = mpBooleanNetwork->getAttractorMatrix();
    mAttractorLength = mpBooleanNetwork->getAttractorLength();
}

ThresholdErgodicSetDifferentiationTree::ThresholdErgodicSetDifferentiationTree(std::string file_path)
{
    mpBooleanNetwork = new RandomBooleanNetwork(file_path);
    mpBooleanNetwork->findAttractors();
    mStochasticMatrix = mpBooleanNetwork->getAttractorMatrix();
    mAttractorLength = mpBooleanNetwork->getAttractorLength();
}


ThresholdErgodicSetDifferentiationTree::~ThresholdErgodicSetDifferentiationTree()
{
    if (mpBooleanNetwork) delete mpBooleanNetwork;
}

DifferentiationTree* ThresholdErgodicSetDifferentiationTree::computeDifferentiationTree() const
{
    DifferentiationTree* differentiation_tree;
    unsigned step = 0;
    std::set<double> thresholds = getThresholdValues();
    unsigned number_of_previous_components = 0;
    //unsigned number_of_tes = 0;
    double threshold = 0.0;
    std::vector<std::map<unsigned,double> > current_stochastic_matrix(mStochasticMatrix);
    bool thresholds_ended = false;
    do
    {
        std::set<std::set<unsigned> > strongly_connected_components = getStronglyConnectedComponents(current_stochastic_matrix);
        unsigned* component = new unsigned[mStochasticMatrix.size()];
        std::set<std::set<unsigned> >::iterator it_components;
        std::set<unsigned>::iterator it_nodes;
        unsigned counter_set = 0;
        for (it_components = strongly_connected_components.begin(); it_components!=strongly_connected_components.end(); ++it_components)
        {
            for (it_nodes = it_components->begin(); it_nodes!=it_components->end(); ++it_nodes)
            {
                component[*it_nodes] = counter_set;
            }
            counter_set++;
        }
        bool isTes = isThresholdErgodicSet(current_stochastic_matrix, component);
        if (step==0)
        {
            if (strongly_connected_components.size() == 1)
            {
                std::set<unsigned> component_states(*strongly_connected_components.begin());
                std::vector<double> stationary_distribution = findStationaryDistribution(current_stochastic_matrix, component_states);
                double cell_cycle_length = this->getCellCycleLength(stationary_distribution, component_states);
                differentiation_tree = new DifferentiationTree(cell_cycle_length, component_states, stationary_distribution);
                number_of_previous_components = 1;
            }
            else
            {
                //Usually it is an error, but we create a "false root".
                if (isTes)
                {
                    std::set<unsigned> component_states;
                    for (unsigned i=0; i<mStochasticMatrix.size(); i++)
                    {
                        component_states.insert(i);
                    }
                    std::vector<double> stationary_distribution = findStationaryDistribution(current_stochastic_matrix, component_states);
                    double cell_cycle_length = this->getCellCycleLength(stationary_distribution, component_states);
                    differentiation_tree = new DifferentiationTree(cell_cycle_length, component_states, stationary_distribution);
                    for (it_components = strongly_connected_components.begin(); it_components!=strongly_connected_components.end(); ++it_components)
                    {
                        std::set<unsigned> component_states(*it_components);
                        unsigned parent = differentiation_tree->searchParentNode(component_states);
                        if (parent < differentiation_tree->size())
                        {
                            //std::set<unsigned> component_states(*strongly_connected_components.begin());
                            std::vector<double> stationary_distribution = findStationaryDistribution(current_stochastic_matrix, component_states);
                            double cell_cycle_length = this->getCellCycleLength(stationary_distribution, component_states);
                            differentiation_tree->addNewChild(parent, cell_cycle_length,
                                    component_states, stationary_distribution);
                        }
                    }
                    number_of_previous_components = strongly_connected_components.size();
                }
            }
        }
        else
        {
            if (isTes && number_of_previous_components != strongly_connected_components.size())
            {
                if (number_of_previous_components == 0)
                {
                    std::set<unsigned> component_states;
                    for (unsigned i=0; i<mStochasticMatrix.size(); i++)
                    {
                        component_states.insert(i);
                    }
                    std::vector<double> stationary_distribution = findStationaryDistribution(current_stochastic_matrix, component_states);
                    double cell_cycle_length = this->getCellCycleLength(stationary_distribution, component_states);
                    differentiation_tree = new DifferentiationTree(cell_cycle_length, component_states, stationary_distribution);
                }
                else
                {
                    for (it_components = strongly_connected_components.begin(); it_components!=strongly_connected_components.end(); ++it_components)
                    {
                        std::set<unsigned> component_states(*it_components);
                        unsigned parent = differentiation_tree->searchParentNode(component_states);
                        if (parent < differentiation_tree->size())
                        {
                            //std::set<unsigned> component_states(*strongly_connected_components.begin());
                            std::vector<double> stationary_distribution = findStationaryDistribution(current_stochastic_matrix, component_states);
                            double cell_cycle_length = this->getCellCycleLength(stationary_distribution, component_states);
                            differentiation_tree->addNewChild(parent, cell_cycle_length, component_states, stationary_distribution);
                        }
                    }
                }
                number_of_previous_components = strongly_connected_components.size();
            }
        }
        delete[] component;
        if (thresholds.size()>0)
        {
            threshold = *thresholds.begin();
            thresholds.erase(thresholds.begin());
            current_stochastic_matrix = getPrunedMatrix(threshold);
        }
        else
        {
            thresholds_ended = true;
        }
        step++;
    } while(!thresholds_ended && number_of_previous_components!=mStochasticMatrix.size());
    return differentiation_tree;
}

std::set<double> ThresholdErgodicSetDifferentiationTree::getThresholdValues() const
{
    std::set<double> thresholds;
    for (unsigned row=0; row< mStochasticMatrix.size(); row++)
    {
        std::map<unsigned,double>::const_iterator iterator;
        for (iterator=mStochasticMatrix.at(row).begin(); iterator!=mStochasticMatrix.at(row).end(); ++iterator)
        {
            thresholds.insert(iterator->second);
        }
    }
    return thresholds;
}

std::vector<std::map<unsigned,double> > ThresholdErgodicSetDifferentiationTree::getPrunedMatrix(
        double threshold) const
{
    std::vector<std::map<unsigned,double> > pruned_matrix;
    for (unsigned row=0; row< mStochasticMatrix.size(); row++)
    {
        std::map<unsigned,double> map;
        pruned_matrix.push_back(map);
        std::map<unsigned,double>::const_iterator iterator;
        for (iterator=mStochasticMatrix.at(row).begin(); iterator!=mStochasticMatrix.at(row).end(); ++iterator)
        {
            if (iterator->second > threshold)
            {
                pruned_matrix.at(row).insert(std::pair<unsigned,double>(iterator->first,iterator->second));
            }
        }
    }
    normalizePrunedMatrix(pruned_matrix);
    return pruned_matrix;
}

void ThresholdErgodicSetDifferentiationTree::normalizePrunedMatrix(
        std::vector<std::map<unsigned,double> > &matrix) const
{
    for (unsigned row=0; row< matrix.size(); row++)
    {
        double row_sum = 0;
        std::map<unsigned,double>::iterator iterator;
        for (iterator=matrix.at(row).begin(); iterator!=matrix.at(row).end(); ++iterator)
        {
            row_sum += iterator->second;
        }
        if (row_sum == 0)
        {
            matrix.at(row)[row] = 1.0; //using [] notation we don't care about previous existing elements.
        } else
        {
            for (iterator=matrix.at(row).begin(); iterator!=matrix.at(row).end(); ++iterator)
            {
                matrix.at(row).at(iterator->first) = iterator->second / row_sum;
            }
        }
    }
}

bool ThresholdErgodicSetDifferentiationTree::isThresholdErgodicSet(
        std::vector<std::map<unsigned,double> > matrix,
        const unsigned* component
        ) const
{
    bool components_connected = false;
    for(unsigned row = 0; row < matrix.size() && !components_connected; row++)
    {
        std::map<unsigned,double>::iterator iterator;
        for (iterator=matrix.at(row).begin(); iterator!=matrix.at(row).end() && !components_connected; ++iterator)
        {
            if (component[row] != component[iterator->first])
            {
                components_connected = true;
            }
        }
    }
    return !components_connected;
}
std::set<std::set<unsigned> > ThresholdErgodicSetDifferentiationTree::getStronglyConnectedComponents(
        std::vector<std::map<unsigned,double> > graph) const
{
    unsigned number_of_nodes = graph.size();
    unsigned time=0;
    unsigned* lowlink = new unsigned[number_of_nodes];
    bool* used = new bool[number_of_nodes];
    for (unsigned i=0; i<number_of_nodes; i++)
    {
        used[i] = false;
        lowlink[i] = 0;
    }
    std::vector<unsigned> stack;

    std::set<std::set<unsigned> > components;
    for (unsigned u=0; u<number_of_nodes; u++)
    {
        if (!used[u])
        {
            findComponentDepthFirstSearch(u, graph, lowlink,
                                    used, stack, time, components);
        }
    }
    delete[] used;
    delete[] lowlink;

    return components;
}

void ThresholdErgodicSetDifferentiationTree::findComponentDepthFirstSearch(
        unsigned u, std::vector<std::map<unsigned,double> > graph, unsigned* pLowlink,
        bool* pUsed, std::vector<unsigned> &rStack, unsigned &rTime, std::set<std::set<unsigned> > &rComponents
        ) const
{
    pLowlink[u] = rTime++;
    pUsed[u] = true;
    rStack.push_back(u);
    bool is_component_root = true;

    std::map<unsigned, double>::iterator iterator;
    for (iterator = graph.at(u).begin(); iterator != graph.at(u).end(); ++iterator)
    {
        unsigned v = iterator->first;
        if (!pUsed[v])
        {
            findComponentDepthFirstSearch(v, graph, pLowlink,
                        pUsed, rStack, rTime, rComponents);
        }
        if (pLowlink[u] > pLowlink[v])
        {
            pLowlink[u] = pLowlink[v];
            is_component_root = false;
        }
    }
    if (is_component_root)
    {
        std::set<unsigned> component;
        unsigned k = 0;
        do {
            k = rStack.back();
            rStack.pop_back();
            component.insert(k);
            pLowlink[k] = std::numeric_limits<unsigned int>::max();
        } while (k != u);
        rComponents.insert(component);
    }
}

std::vector<double> ThresholdErgodicSetDifferentiationTree::findStationaryDistribution(
        std::vector<std::map<unsigned, double> > tes_map_matrix,
        std::set<unsigned> component) const
{
    //Mat ergodic_matrix;
    unsigned size = component.size();
    if (size == 1)
    {
        std::vector<double> stationary_distribution(1,1.0);
        return stationary_distribution;
    }
    LinearSystem ls(size, size);
    ls.SetMatrixIsConstant(true);

    std::set<unsigned>::iterator set_iterator;
    unsigned row = 0;
    for (set_iterator = component.begin(); set_iterator!=component.end(); ++set_iterator)
    {
        if (row != size - 1)
        {
            ls.SetRhsVectorElement(row,0.0);
        }
        else
        {
            ls.SetRhsVectorElement(row,1.0);
        }
        unsigned state_value = *set_iterator;
        std::map<unsigned, double> map_row = tes_map_matrix.at(state_value);
        std::map<unsigned, double>::iterator map_iterator = map_row.begin();
        for(unsigned col = 0; col<size; col++)
        {
            if (col == size-1)
            {
                ls.SetMatrixElement(col, row, 1);
            }
            else
            {
                if (map_iterator != map_row.end() && map_iterator->first == col)
                {
                    if (col == row)
                    {
                        ls.SetMatrixElement(col, row, map_iterator->second - 1.0);
                    }
                    else
                    {
                        ls.SetMatrixElement(col, row, map_iterator->second);
                    }
                    map_iterator++;
                }
                else
                {
                    ls.SetMatrixElement(col, row, 0.0);
                }
            }
        }
        row++;
    }
    ls.AssembleFinalLinearSystem();
    /* Visualization */
    //ls.DisplayMatrix();
    //ls.DisplayRhs();
    /* End visualization */
    Vec solution_vector;
    solution_vector = ls.Solve();
    double* p_solution_elements_array;
    VecGetArray(solution_vector, &p_solution_elements_array);
    std::vector<double> stationary_distribution(p_solution_elements_array, p_solution_elements_array + size);
    VecRestoreArray(solution_vector, &p_solution_elements_array);
    PetscTools::Destroy(solution_vector);
    //delete p_solution_elements_array;
    return stationary_distribution;
}

double ThresholdErgodicSetDifferentiationTree::getCellCycleLength(
        std::vector<double> stationary_distribution, std::set<unsigned> component) const
{
    double length = 0;
    std::set<unsigned>::iterator set_iterator;
    unsigned i=0;
    for (set_iterator = component.begin(); set_iterator != component.end();++set_iterator)
    {
        assert(*set_iterator < mAttractorLength.size());
        length += ((double) mAttractorLength.at(*set_iterator)) * stationary_distribution.at(i);
        i++;
    }
    return length;
}

const RandomBooleanNetwork* ThresholdErgodicSetDifferentiationTree::getBooleanNetwork() const
{
    return mpBooleanNetwork;
}

std::vector<unsigned> ThresholdErgodicSetDifferentiationTree::getAttractorLength() const
{
    return mAttractorLength;
}
std::vector<std::map<unsigned,double> > ThresholdErgodicSetDifferentiationTree::getStochasticMatrix() const
{
    return mStochasticMatrix;
}
DifferentiationTree* ThresholdErgodicSetDifferentiationTree::getDifferentiationTree() const
{
    DifferentiationTree* differentiation_tree = computeDifferentiationTree();
    differentiation_tree->colorise();

    return differentiation_tree;
}
