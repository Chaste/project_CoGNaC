#include "DifferentiationTree.hpp"
#include <cassert>
#include <math.h>

DifferentiationTree::DifferentiationTree()
{
    mVertices.push_back(new DifferentiationTreeNode());
    mLevelNodes.push_back(std::set<unsigned>());
    mLevelNodes.at(0).insert(0);
    mLeaves.insert(0);
    mNodesLevel.push_back(0);
}

DifferentiationTree::DifferentiationTree(double cell_cycle_length)
{
    assert(cell_cycle_length > 0.0);

    mVertices.push_back(
            new DifferentiationTreeNode(cell_cycle_length));
    mLevelNodes.push_back(std::set<unsigned>());
    mLevelNodes.at(0).insert(0);
    mLeaves.insert(0);
    mNodesLevel.push_back(0);

}

DifferentiationTree::DifferentiationTree(
        double cell_cycle_length,
        std::set<unsigned> component_states,
        std::vector<double> stationary_distribution
        )
{
    assert(cell_cycle_length > 0.0);

    mVertices.push_back(
            new DifferentiationTreeNode(
                    cell_cycle_length,
                    component_states,
                    stationary_distribution));
    mLevelNodes.push_back(std::set<unsigned>());
    mLevelNodes.at(0).insert(0);
    mLeaves.insert(0);
    mNodesLevel.push_back(0);
}

DifferentiationTree::~DifferentiationTree()
{
    if (mVertices.size() > 0)
    {
        for (unsigned i=mVertices.size() -1;
                mVertices.size() > 0;
                i = mVertices.size() -1)
        {
            if (mVertices.at(i))
            {
                delete mVertices[i];
            }
            mVertices.pop_back();
        }
    }
}

bool DifferentiationTree::checkIsomorphism(
            unsigned step,
            const DifferentiationTree* tree1,
            const DifferentiationTree* tree2,
            const std::vector<std::set<unsigned> >& levelNodes1,
            const std::vector<std::set<unsigned> >& levelNodes2,
            std::map<unsigned, unsigned>& sub_tree_type1,
            std::map<unsigned, unsigned>& sub_tree_type2,
            std::vector<std::list<unsigned> >& tuples1,
            std::vector<std::list<unsigned> >& tuples2
            ) const
{
    /*
     * Taken from The design and analysis of computer algorithms,
     * Aho, Hopcroft, Ullman, Pagg. 84 - 85 - 86.
     */
    assert(tree1 && tree2);
    std::map<unsigned, unsigned>::iterator map_it1 = sub_tree_type1.find(0);
    std::map<unsigned, unsigned>::iterator map_it2 = sub_tree_type2.find(0);
    if (step == 0)
    {
        if (map_it1 != sub_tree_type1.end() && map_it2 != sub_tree_type1.end())
        {
            if (*map_it1 == *map_it2) return true;
            else return false;
        }
        else return false;
    }
    else
    {
        if (levelNodes1.at(step).size() == levelNodes2.at(step).size())
        {
            std::list<std::list<unsigned> > list_of_tuples1;
            std::list<std::list<unsigned> > list_of_tuples2;
            std::map<std::list<unsigned>, std::set<unsigned> > map_list1;
            std::map<std::list<unsigned>, std::set<unsigned> > map_list2;

            // T1
            for (std::set<unsigned>::iterator set_it = levelNodes1.at(step).begin();
                    set_it != levelNodes1.at(step).end(); ++set_it)
            {
                unsigned parent = tree1->getNode(*set_it)->getParent();
                tuples1.at(parent).push_back(sub_tree_type1.at(*set_it));
            }
            for (std::set<unsigned>::iterator set_it = levelNodes1.at(step-1).begin();
                    set_it != levelNodes1.at(step-1).end(); ++set_it)
            {
                if (!tree1->getNode(*set_it)->isLeaf())
                {
                    tuples1.at(*set_it).sort();
                    list_of_tuples1.push_back(tuples1.at(*set_it));

                    std::set<unsigned> temp_set;
                    temp_set.insert(*set_it);
                    std::pair<std::map<std::list<unsigned>,
                        std::set<unsigned> >::iterator,bool> ret;
                    ret = map_list1.insert (std::pair<std::list<unsigned>,
                            std::set<unsigned> >
                            (tuples1.at(*set_it), temp_set));
                    if (ret.second==false) {
                        ret.first->second.insert(*set_it);
                    }

                }
            }
            list_of_tuples1.sort();
            //Same for T2
            for (std::set<unsigned>::iterator set_it = levelNodes2.at(step).begin();
                    set_it != levelNodes2.at(step).end(); ++set_it)
            {
                unsigned parent = tree2->getNode(*set_it)->getParent();
                tuples2.at(parent).push_back(sub_tree_type2.at(*set_it));
            }
            for (std::set<unsigned>::iterator set_it = levelNodes2.at(step-1).begin();
                    set_it != levelNodes2.at(step-1).end(); ++set_it)
            {
                if (!tree2->getNode(*set_it)->isLeaf())
                {
                    tuples2.at(*set_it).sort();
                    list_of_tuples2.push_back(tuples2.at(*set_it));

                    std::set<unsigned> temp_set;
                    temp_set.insert(*set_it);
                    std::pair<std::map<std::list<unsigned>,
                        std::set<unsigned> >::iterator,bool> ret;
                    ret = map_list2.insert (std::pair<std::list<unsigned>,
                            std::set<unsigned> >
                            (tuples2.at(*set_it), temp_set));
                    if (ret.second==false) {
                        ret.first->second.insert(*set_it);
                    }
                }
            }
            list_of_tuples2.sort();

            if (list_of_tuples1 == list_of_tuples2)
            {
                unsigned position = 1;
                unsigned i=0;
                unsigned nodes_equals=0;
                for (std::list<std::list<unsigned> >::iterator list_it = list_of_tuples1.begin();
                        list_it != list_of_tuples1.end();++list_it)
                {
                    if (i >= nodes_equals)
                    {
                        std::set<unsigned> list_nodes1 = map_list1.at(*list_it);
                        for (std::set<unsigned>::iterator set_it = list_nodes1.begin();
                                set_it != list_nodes1.end(); ++set_it)
                        {
                            sub_tree_type1.insert(std::pair<unsigned, unsigned>(*set_it,position));
                        }

                        std::set<unsigned> list_nodes2 = map_list2.at(*list_it);
                        for (std::set<unsigned>::iterator set_it = list_nodes2.begin();
                                set_it != list_nodes2.end(); ++set_it)
                        {
                            sub_tree_type2.insert(std::pair<unsigned, unsigned>(*set_it,position));
                        }
                        position++;
                    }
                    else
                    {
                        i++;
                    }
                }
                step--;
                return checkIsomorphism(step, tree1, tree2, levelNodes1, levelNodes2, sub_tree_type1, sub_tree_type2, tuples1, tuples2);
            }
            else return false;
        }
        else return false;
    }
}

unsigned DifferentiationTree::searchParentNode(
        std::set<unsigned> child_component_states,
        unsigned starting_node) const
{
    assert(starting_node < mVertices.size());
    /* Check the set inclusion of child component states in parent component state */
    std::set<unsigned> parent_component_states = mVertices.at(starting_node)->getComponentStates();
    if(std::includes(parent_component_states.begin(), parent_component_states.end(),
            child_component_states.begin(), child_component_states.end()))
    {
        std::vector<unsigned> children = mVertices.at(starting_node)->getChildren();
        std::vector<unsigned>::iterator iterator = children.begin();
        unsigned nodeFinded;
        for(nodeFinded = mVertices.size(); iterator!=children.end() && (nodeFinded==mVertices.size()); ++iterator)
        {
            nodeFinded = searchParentNode(child_component_states, *iterator);
        }
        if(nodeFinded < mVertices.size())
            return nodeFinded;
        else
            return starting_node;
    }
    else
    {
        return mVertices.size();
    }
}

void DifferentiationTree::addChild(unsigned parent, DifferentiationTreeNode* node)
{
    assert(parent < mVertices.size());
    assert(node);
    std::set<unsigned> parent_component_states = mVertices.at(parent)->getComponentStates();
    std::set<unsigned> child_component_states = node->getComponentStates();
    if(std::includes(parent_component_states.begin(), parent_component_states.end(),
                child_component_states.begin(), child_component_states.end()))
    {
        mVertices.push_back(node);
        if (mVertices.at(parent)->getNumberOfChildren() == 0)
        {
            std::set<unsigned>::iterator it = mLeaves.find(parent);
            if (it != mLeaves.end())
            {
                mLeaves.erase(it);
            }
        }
        mVertices.at(parent)->addChild(mVertices.size()-1);
        mVertices.at(mVertices.size() -1)->setParent(parent);
        mLeaves.insert(mVertices.size() - 1);
        unsigned parent_level = mNodesLevel.at(parent);
        mNodesLevel.push_back(parent_level + 1);
        // parent_level must not be > mLevelNodes.size()
        if (parent_level + 1 == mLevelNodes.size())
            mLevelNodes.push_back(std::set<unsigned>());

        mLevelNodes.at(parent_level + 1).insert(mVertices.size()-1);
        std::vector<double> new_distribution(
                mVertices.at(parent)->getNumberOfChildren(),
                1.0 / (double) mVertices.at(parent)->getNumberOfChildren());
        mVertices.at(parent)->setStationaryDistribution(new_distribution);
    }
}

void DifferentiationTree::addNode(DifferentiationTreeNode* node)
{
    assert(node);
    if (mVertices.size() > 0)
    {
        unsigned parent = searchParentNode(node->getComponentStates());
        if (parent < mVertices.size())
        {
            addChild(parent, node);
        }
        //Else no add.
    }
}

void DifferentiationTree::addNewChild(unsigned parent)
{
    assert(parent < mVertices.size());

    DifferentiationTreeNode* node = new DifferentiationTreeNode();
    addChild(parent, node);

}

void DifferentiationTree::addNewChild(unsigned parent, double cell_cycle_length)
{
    assert(parent < mVertices.size());
    assert(cell_cycle_length > 0.0);

    DifferentiationTreeNode* node = new DifferentiationTreeNode(cell_cycle_length);
    addChild(parent, node);
}

void DifferentiationTree::addNewChild(
        unsigned parent,
        double cell_cycle_length,
        std::set<unsigned> component_states,
        std::vector<double> stationary_distribution
        )
{
    assert(parent < mVertices.size());
    assert(cell_cycle_length > 0.0);

    DifferentiationTreeNode* node =
            new DifferentiationTreeNode(
                    cell_cycle_length,
                    component_states,
                    stationary_distribution);

	std::set<unsigned> parent_component_states = mVertices.at(parent)->getComponentStates();
	std::set<unsigned> child_component_states = node->getComponentStates();
	if(std::includes(parent_component_states.begin(), parent_component_states.end(),
				child_component_states.begin(), child_component_states.end()))

	mVertices.push_back(node);
	if (mVertices.at(parent)->getNumberOfChildren() == 0)
	{
		std::set<unsigned>::iterator it = mLeaves.find(parent);
		if (it != mLeaves.end())
		{
			mLeaves.erase(it);
		}
	}
	mVertices.at(parent)->addChild(mVertices.size()-1);
	mVertices.at(mVertices.size() -1)->setParent(parent);
	mLeaves.insert(mVertices.size() - 1);

	unsigned parent_level = mNodesLevel.at(parent);
	mNodesLevel.push_back(parent_level + 1);
	// parent_level must not be > mLevelNodes.size()
	if (parent_level + 1 == mLevelNodes.size())
		mLevelNodes.push_back(std::set<unsigned>());

	mLevelNodes.at(parent_level + 1).insert(mVertices.size()-1);
}

void DifferentiationTree::normaliseLength(double hours_mean)
{
    assert(hours_mean > 0.0);

    if (mVertices.size() > 0)
    {
        double length_mean = 0.0;
        for (unsigned i=0;i<mVertices.size();i++)
        {
            length_mean += mVertices.at(i)->getCellCycleLength();
        }
        length_mean /= (double) mVertices.size();
        for (unsigned i=0;i<mVertices.size();i++)
        {
            double normalized_length = mVertices.at(i)->getCellCycleLength() * hours_mean / length_mean;
            mVertices.at(i)->setCellCycleLength(normalized_length);
        }
    }
}

void DifferentiationTree::colorise()
{
    if (mVertices.size() > 1)
    {
        unsigned number_of_differentiated = mLeaves.size();
        double transit_step_colour = (.6 - .3) / 2.0;
        if (mVertices.size() > 2 + number_of_differentiated)
        {
            transit_step_colour = ((double) .6 - .3) / (double) (mVertices.size()
                    -2 - number_of_differentiated);
        }
        double differentiated_step_colour = (1.0 - .9) / 2.0;
        if (number_of_differentiated > 1)
        {
            differentiated_step_colour = ((double) 1.0 - .9) /
                    ((double) number_of_differentiated - 1.0) ;
        }
        //BFS
        std::vector<unsigned> queue;
        queue.push_back(0);
        unsigned differentiated_count = 0;
        unsigned transit_count = 0;
        while (!queue.empty())
        {
            unsigned node = queue.at(queue.size()-1);
            if (node == 0)
            {
                mVertices.at(node)->setColour(0.0);
            }
            else
            {
                if(mVertices.at(node)->isLeaf())
                {
                    double colour = ((double) differentiated_count) * differentiated_step_colour;
                    differentiated_count++;
                    mVertices.at(node)->setColour(.9 + colour);
                }
                else
                {
                    double colour = ((double) transit_count) * transit_step_colour;
                    transit_count++;
                    mVertices.at(node)->setColour(.3 + colour);
                }
            }
            std::vector<unsigned> children = mVertices.at(node)->getChildren();
            queue.pop_back();

            for (unsigned i=0; i<children.size(); i++)
            {
                queue.push_back(children.at(i));
            }
        }
    }
}

void DifferentiationTree::setColour(unsigned node_id, double colour)
{
    if (node_id >= mVertices.size())
        EXCEPTION("Id of the node is greater or equal to the number of vertices");
    mVertices.at(node_id)->setColour(colour);
}

void DifferentiationTree::printDifferentiationTree() const
{
    for (unsigned i=0; i<mVertices.size(); i++)
    {
        std::set<unsigned>::iterator iterator;
        std::cout << "- ID: " << i << ", parent: " <<
                mVertices.at(i)->getParent() << ", with " <<
                mVertices.at(i)->getNumberOfChildren() << " children.\n";
        std::cout << "---states (set): {";
        std::set<unsigned> states = mVertices.at(i)->getComponentStates();
        for (iterator=states.begin(); iterator!=states.end(); ++iterator)
        {
            std::cout << " " << *iterator << " ";
        }
        std::cout << "}\n";
    }
}

void DifferentiationTree::printDifferentiationTreeToGmlFile(std::string directory, std::string filename) const
{
	if (directory.empty())
		EXCEPTION("Directory name not valid.");

	if (filename.size() < 5 || filename.compare(filename.size()-4,4,".gml") != 0)
		EXCEPTION("File path not valid. It must terminate with '.gml' extension.");

	OutputFileHandler handler(directory,false);
	out_stream p_file = handler.OpenOutputFile(filename);

	*p_file << "graph [\n";
	*p_file << " comment \"Graph generated by CoGNaC\"\n";
	*p_file << " directed 1\n";
	unsigned square_position = 0;
	float x = 0.0, y = 0.0, step = 80.0, square_root = 0.0;
	square_root = sqrt((float) mVertices.size());
	if (ceilf(square_root) == square_root) square_position = (unsigned) square_root;
	else square_position = (unsigned) square_root + 1;

	/* for each node */
    for (unsigned i=0; i<mVertices.size(); i++)
    {
    	if (i > 0 && i % square_position == 0)
		{
			y += step;
			x = 0.0;
		}
		*p_file << " node [\n";
		*p_file << "  id " << i+1 << "\n";
		*p_file << "  label \"" << i+1 << "\"\n";
		*p_file << "  cell_cycle_length " << mVertices.at(i)->getCellCycleLength() << "\n";
		*p_file << "  graphics [\n   x " << x << "\n   y " << y << "\n";
		*p_file << "   w 40.0\n   h 40.0\n   type \"ellipse\"\n  ]\n";
		*p_file << "  LabelGraphics [\n   text \"" << i+1 << "\"\n  ]\n";
		*p_file << " ]\n";
		x += step;
    }

    for (unsigned i=0; i<mVertices.size(); i++)
	{
    	if (mVertices.at(i)->getChildren().size() > 0)
    	{
    	    std::vector<unsigned> children = mVertices.at(i)->getChildren();
    	    std::vector<double> probabilities = mVertices.at(i)->getStationaryDistribution();
    	    for (unsigned j=0; j<children.size(); j++)
    	    {
    		*p_file << " edge [\n";
    		*p_file << "  source " << i + 1 << "\n";
    		*p_file << "  target " << children.at(j) + 1 << "\n";
    		*p_file << "  label \" \"\n";
    		*p_file << "  probability " << probabilities.at(j) << "\n";
    		*p_file << "  graphics [\n   arrow \"last\"\n  ]\n";
    		*p_file << " ]\n";
    	    }
    	}
	}
    *p_file << "]";
    p_file->close();
}

bool DifferentiationTree::topologyTreeCompare(const DifferentiationTree* user_defined_tree) const
{
    if (this->size() == user_defined_tree->size())
    {
        std::set<unsigned> leaves1 = mLeaves;
        std::set<unsigned> leaves2 = user_defined_tree->getLeaves();

        if (leaves1.size() == leaves2.size())
        {
            if (mLevelNodes.size() == user_defined_tree->getLevelNodes().size())
            {
                if (mLevelNodes.size() == 1) return true;
                std::vector<std::set<unsigned> > levelNodes1 = mLevelNodes;
                std::vector<std::set<unsigned> > levelNodes2 = user_defined_tree->getLevelNodes();

                std::map<unsigned, unsigned> sub_tree_type1;
                std::map<unsigned, unsigned> sub_tree_type2;
                std::set<unsigned>::iterator set_it;
                for (set_it = leaves1.begin(); set_it != leaves1.end(); ++set_it)
                {
                    sub_tree_type1.insert(std::pair<unsigned, unsigned>(*set_it, 0));
                }
                for (set_it = leaves2.begin(); set_it != leaves2.end(); ++set_it)
                {
                    sub_tree_type2.insert(std::pair<unsigned, unsigned>(*set_it, 0));
                }
                std::vector<std::list<unsigned> > tuples1(this->size(), std::list<unsigned>());
                std::vector<std::list<unsigned> > tuples2(this->size(), std::list<unsigned>());

                return checkIsomorphism(mLevelNodes.size()-1, this, user_defined_tree, levelNodes1,levelNodes2,sub_tree_type1,sub_tree_type2, tuples1, tuples2);
            }
            else return false;
        }
        else return false;
    }
    else return false;
}

std::set<unsigned> DifferentiationTree::getLeaves() const
{
    return mLeaves;
}

std::vector<std::set<unsigned> > DifferentiationTree::getLevelNodes() const
{
    return mLevelNodes;
}

std::vector<unsigned> DifferentiationTree::getNodesLevel() const
{
    return mNodesLevel;
}

const DifferentiationTreeNode* DifferentiationTree::getRoot() const
{
    assert(mVertices.size()>0);
    return mVertices.at(0);
}

const DifferentiationTreeNode* DifferentiationTree::getNode(unsigned node_id) const
{
    assert(node_id < mVertices.size());
    return mVertices.at(node_id);
}

unsigned DifferentiationTree::size() const{
    return mVertices.size();
}
