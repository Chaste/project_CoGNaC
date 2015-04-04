#ifndef TESTSEARCHINGGENEACTIVATIONPATTERNSINTHELPERNETWORK_HPP_
#define TESTSEARCHINGGENEACTIVATIONPATTERNSINTHELPERNETWORK_HPP_

#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "SmartPointers.hpp"

#include <bdd.h>
#include "RandomBooleanNetwork.hpp"
#include "DifferentiationTree.hpp"
#include "DifferentiationTreeNode.hpp"

#include "ThresholdErgodicSetDifferentiationTree.hpp"
#include "FakePetscSetup.hpp"

class TestSearchingGeneActivationPatternsInThelperNetwork : public CxxTest::TestSuite
{
public:
    void testThelper() throw (Exception)
    {
        bdd_init(10000,1000);
        try
        {
        	RandomBooleanNetwork *rbn1 = new RandomBooleanNetwork("projects/SimoneRubinacci/networks_samples/thelper.cnet");
            rbn1->findAttractors();

            std::vector<std::map<unsigned,double> >  atn = rbn1->getAttractorMatrix();
            for (unsigned row=0; row< atn.size(); row++)
			{
				std::map<unsigned,double>::iterator iterator;
				for (iterator=atn.at(row).begin(); iterator!=atn.at(row).end(); ++iterator)
				{
					std::cout << "(" << row << "," << iterator->first << ") -> " << iterator->second << "\n";
				}
			}

        	ThresholdErgodicSetDifferentiationTree TES_tree ("projects/SimoneRubinacci/networks_samples/thelper.net");

            TS_ASSERT_EQUALS(TES_tree.getBooleanNetwork()->getAttractorsNumber(),3u);
            std::vector<unsigned> attractors_lengths = TES_tree.getBooleanNetwork()->getAttractorLength();
            DifferentiationTree* diff_tree = TES_tree.getDifferentiationTree();
            diff_tree->printDifferentiationTree();

        }
        catch (Exception& e)
        {
            TS_FAIL(e.GetMessage());
            bdd_done();
        }
        bdd_done();
    }
    void testpippolo()
    {
    	std::vector<std::map<unsigned,double> > stoc_matrix(3);
        stoc_matrix.at(0).insert(std::pair<unsigned,double>(2,1));
        stoc_matrix.at(1).insert(std::pair<unsigned,double>(1,1));
        stoc_matrix.at(2).insert(std::pair<unsigned,double>(0,1));

        ThresholdErgodicSetDifferentiationTree TES_tree(stoc_matrix,
                        std::vector<unsigned>(3,1));

        DifferentiationTree* diff_tree = TES_tree.getDifferentiationTree();
        diff_tree->printDifferentiationTree();

        delete diff_tree;
    }
};



#endif /* TESTSEARCHINGGENEACTIVATIONPATTERNSINTHELPERNETWORK_HPP_ */
