#include "ThresholdErgodicSetDifferentiationTree.hpp"
#include <limits>
#include <iostream>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

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
	if (file_path.size() > 4)
	{
		try {
			if (file_path.compare(file_path.size()-4,4,".net") == 0
					||  file_path.compare(file_path.size()-4,4,"cnet") == 0)
			{
				mpBooleanNetwork = new RandomBooleanNetwork(file_path);
				mpBooleanNetwork->findAttractors();
				mStochasticMatrix = mpBooleanNetwork->getAttractorMatrix();
				mAttractorLength = mpBooleanNetwork->getAttractorLength();
			}
			else if (file_path.compare(file_path.size()-4,4,".dat") == 0)
			{
				std::vector<std::map<unsigned,double> > matrix;
				std::vector<unsigned> attractor_length;

			    std::vector<std::string> strs;
				std::string line;
				std::ifstream input_file (file_path.c_str());
				try
				{
					if (!input_file.is_open()) EXCEPTION("Not able to open the file.");
					unsigned matrix_length = 0;
					if (!std::getline(input_file,line)) throw;
					boost::algorithm::trim(line);
					boost::algorithm::to_lower(line);
					boost::split(strs, line, boost::is_any_of(","), boost::algorithm::token_compress_on);
					matrix_length = strs.size();
					if (strs.size() != 0)
					{
						double value = 0.0;
						unsigned i=0;
						do {
							if (i!=0)
							{
								if (!std::getline(input_file,line)) throw;
								boost::algorithm::trim(line);
								boost::algorithm::to_lower(line);
								boost::split(strs, line, boost::is_any_of(","), boost::algorithm::token_compress_on);
								if (strs.size() != matrix_length) throw;
							}
							unsigned j=0;
							std::map<unsigned,double> map;
							matrix.push_back(map);
							do {
								value = boost::lexical_cast<double>(strs.at(j));
								if (value != 0.0)
								{
									matrix.at(i).insert(std::pair<unsigned,double>(j,value));
								}
								j++;
							} while (j < matrix_length);
							i++;
						} while (i < matrix_length);
					} else throw;
					if (!std::getline(input_file,line)) throw;
					boost::algorithm::trim(line);
					if (!line.empty()) throw;
					if (!std::getline(input_file,line)) throw;
					boost::algorithm::trim(line);
					boost::algorithm::to_lower(line);
					boost::split(strs, line, boost::is_any_of(","), boost::algorithm::token_compress_on);
					if (strs.size() != matrix_length) throw;
					unsigned j=0;
					do {
						attractor_length.push_back(boost::lexical_cast<unsigned>(strs.at(j)));
						j++;
					} while (j < matrix_length);
					mStochasticMatrix = matrix;
					mAttractorLength = attractor_length;
					mpBooleanNetwork = NULL;
				} catch (const Exception&)
				{
					input_file.close();
					throw;
				}
		} else EXCEPTION("File format is not correct.");
		} catch (Exception& e)
		{
			EXCEPTION("Error reading the file.");
		}
	} else EXCEPTION("Error in the file path.");
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
    double threshold = 0.0;
    std::vector<std::map<unsigned,double> > current_stochastic_matrix(mStochasticMatrix);
    bool thresholds_ended = false;
    do
    {

        std::set<std::set<unsigned> > terminal_sccs = getTerminalStronglyConnectedComponents(current_stochastic_matrix);
        std::set<std::set<unsigned> >::iterator it_components;

        bool isTes = isThresholdErgodicSet(current_stochastic_matrix, terminal_sccs);
        if (step==0)
        {
			std::set<unsigned> component_states;
			for (unsigned i=0; i<mStochasticMatrix.size(); i++)
			{
				component_states.insert(i);
			}
			differentiation_tree = new DifferentiationTree(component_states);
			differentiation_tree->setThreshold(0,0.0);
			if (!isTes) differentiation_tree->setFakeRoot(true);
			number_of_previous_components = 1;
			if (isTes && terminal_sccs.size() > 1)
			{
				for (it_components = terminal_sccs.begin(); it_components!=terminal_sccs.end(); ++it_components)
				{
					std::set<unsigned> component_states(*it_components);
					unsigned parent = differentiation_tree->searchParentNode(component_states);
					if (parent < differentiation_tree->size() && differentiation_tree->getNode(parent)->getComponentStates()!=component_states)
					{
						differentiation_tree->addNewChild(parent, component_states);
						differentiation_tree->setThreshold(differentiation_tree->size()-1,threshold);
					}
				}
				number_of_previous_components = terminal_sccs.size();
			}
        }
        else
		{
			if (isTes && number_of_previous_components != terminal_sccs.size())
			{
				for (it_components = terminal_sccs.begin(); it_components!=terminal_sccs.end(); ++it_components)
				{
					std::set<unsigned> component_states(*it_components);
					unsigned parent = differentiation_tree->searchParentNode(component_states);
					if (parent < differentiation_tree->size() && differentiation_tree->getNode(parent)->getComponentStates()!=component_states)
					{
						differentiation_tree->addNewChild(parent, component_states);
						differentiation_tree->setThreshold(differentiation_tree->size()-1,threshold);
					}
				}
				number_of_previous_components = terminal_sccs.size();
			}
		}
		if (thresholds.size()>0)
		{
			threshold = *thresholds.begin();
			if (threshold != 1.0)
			{
				thresholds.erase(thresholds.begin());
				current_stochastic_matrix = getPrunedMatrix(threshold);
				bool rows_equals_to_zero = false;
				for (unsigned row=0; row< current_stochastic_matrix.size() && !rows_equals_to_zero; row++)
				{
					double row_sum = 0;
					std::map<unsigned,double>::iterator iterator;
					rows_equals_to_zero = current_stochastic_matrix.at(row).empty();
					for (iterator=current_stochastic_matrix.at(row).begin(); iterator!=current_stochastic_matrix.at(row).end() && !rows_equals_to_zero; ++iterator)
					{
						row_sum += iterator->second;
					}
					if (row_sum < 0.99)
					{
						rows_equals_to_zero = true;
					}
				}
				thresholds_ended = rows_equals_to_zero;
			} else	thresholds_ended = true;
		}
		else
		{
			thresholds_ended = true;
		}
		step++;
	} while(!thresholds_ended);
    assignProbabilitiesAndCellCycleLengths(differentiation_tree);
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
        if (row_sum != 0)
        {
        	for (iterator=matrix.at(row).begin(); iterator!=matrix.at(row).end(); ++iterator)
            {
                matrix.at(row).at(iterator->first) = iterator->second / row_sum;
            }
        }
    }
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

std::set<std::set<unsigned> > ThresholdErgodicSetDifferentiationTree::getTerminalStronglyConnectedComponents(
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

    std::set<std::set<unsigned> >::iterator it_components;
    std::set<unsigned>::iterator it_nodes;
    for (it_components = components.begin(); it_components != components.end(); ++it_components)
    {
    	bool is_a_terminal_scc = true;

    	for (it_nodes = it_components->begin(); it_nodes != it_components->end() && is_a_terminal_scc; ++it_nodes)
		{
    		std::map<unsigned,double> map = graph[*it_nodes];
    		std::map<unsigned,double>::iterator it_map;

    		for(it_map = map.begin();it_map!=map.end() && is_a_terminal_scc; ++it_map)
    		{
    			//each reachable node from a member of a terminal scc must be a member of the terminal scc
    			if (it_components->find(it_map->first) == it_components->end())
    			{
    				is_a_terminal_scc = false;
    			}
    		}
		}
    	if (!is_a_terminal_scc)
    	{
    		components.erase(*it_components);
    	}
    }

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

bool ThresholdErgodicSetDifferentiationTree::isThresholdErgodicSet(
        std::vector<std::map<unsigned,double> > matrix,
        const std::set<std::set<unsigned> > components
        ) const
{
    unsigned attractors_in_terminal_components = 0;
    std::set<std::set<unsigned> >::const_iterator comp_iterator;
    for(comp_iterator = components.begin();
    		comp_iterator != components.end();
    		++comp_iterator)
    {
    	attractors_in_terminal_components += comp_iterator->size();
    }
    return attractors_in_terminal_components==matrix.size();
}

void ThresholdErgodicSetDifferentiationTree::assignProbabilitiesAndCellCycleLengths(
		DifferentiationTree* differentiation_tree) const
{
	for (unsigned i=0; i<differentiation_tree->size(); i++)
	{
		if (i==0 && differentiation_tree->hasFakeRoot()) continue;
		double threshold = differentiation_tree->getNode(i)->getThreshold();
	    std::vector<std::map<unsigned,double> > current_stochastic_matrix(mStochasticMatrix);
	    if (threshold < 1.0)
		{
			current_stochastic_matrix = getPrunedMatrix(threshold);
			std::set<unsigned> component_states = differentiation_tree->getNode(i)->getComponentStates();
			std::vector<double> stationary_distribution = findStationaryDistribution(current_stochastic_matrix, component_states);
			double cell_cycle_length = this->getCellCycleLength(stationary_distribution, component_states);
			std::vector<double> stoc_differentiation_vector;
			std::vector<unsigned> children = differentiation_tree->getNode(i)->getChildren();
			std::set<unsigned>::iterator child_iterator;
			for (unsigned child_count = 0; child_count<children.size();++child_count)
			{
				double probability = 0.0;

				std::set<unsigned> component_states_child = differentiation_tree->getNode(children.at(child_count))->getComponentStates();
				std::set<unsigned>::iterator set_iterator;
				for (set_iterator = component_states_child.begin();
						set_iterator!=component_states_child.end();
						set_iterator++)
				{
					probability += stationary_distribution.at(*set_iterator);
				}
				stoc_differentiation_vector.push_back(probability);
			}
			differentiation_tree->setStationaryDistributionAndCellCycleLength(
					i,
					stationary_distribution,
					cell_cycle_length);
			differentiation_tree->setDifferentiationProbability(
					i, stoc_differentiation_vector);
		}
	    else
	    {
	    	std::set<unsigned> component_states = differentiation_tree->getNode(i)->getComponentStates();
	    	std::vector<double> stationary_distribution = findStationaryDistribution(current_stochastic_matrix, component_states);
	    	double cell_cycle_length = this->getCellCycleLength(stationary_distribution, component_states);
			differentiation_tree->setStationaryDistributionAndCellCycleLength(
					i,
					stationary_distribution,
					cell_cycle_length);
			std::vector<double> stoc_differentiation_vector;
			/** This is a fake probability. */
			//stoc_differentiation_vector.push_back(1.0);
			//differentiation_tree->setDifferentiationProbability(i, stoc_differentiation_vector);
	    }
	}
}

void ThresholdErgodicSetDifferentiationTree::printStochasticMatrixAndAttractorLengthsToDatFile(
		std::string directory, std::string filename) const
{
	if (directory.empty())
		EXCEPTION("Directory name not valid.");

	if (filename.size() < 5 || filename.compare(filename.size()-4,4,".dat") != 0)
		EXCEPTION("File path not valid. It must terminate with '.dat' extension.");

	OutputFileHandler handler(directory,false);
	out_stream p_file = handler.OpenOutputFile(filename);

	for (unsigned row=0; row<mStochasticMatrix.size();row++)
	{
		unsigned col=0;
        std::map<unsigned,double>::const_iterator iterator;
        for (iterator=mStochasticMatrix.at(row).begin(); iterator!=mStochasticMatrix.at(row).end(); ++iterator)
        {
			for (;col<iterator->first;col++)
			{
				if (col!=0) *p_file << ",";
				*p_file << 0.0;
			}
			if (col!=0) *p_file << ",";
			*p_file << iterator->second;
			col++;
		}
		for (;col<mStochasticMatrix.size();col++)
		{
			if (col!=0) *p_file << ",";
			*p_file << 0.0;
		}
		*p_file << "\n";
	}
	*p_file << "\n";
	for (unsigned i=0; i<mAttractorLength.size(); i++)
	{
		if (i!=0) *p_file << ",";
		*p_file << mAttractorLength.at(i);
	}
	*p_file << "\n";
	p_file->close();
}
