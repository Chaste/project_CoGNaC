#ifndef TESTSEARCHINGGENEACTIVATIONPATTERNSINTHELPERNETWORK_HPP_
#define TESTSEARCHINGGENEACTIVATIONPATTERNSINTHELPERNETWORK_HPP_

/*
 * = Searching gene activation patterns in T-helper cell differentiation =
 * In this class we show how to search gene activation patterns starting from a network,
 * calculate its Attractor Transition Network (Figure 3) and then calculate its {{{DifferentiationTree}}}.
 *
 * === Including header files ===
 * We begin by including the necessary header files.
*/

#include <cxxtest/TestSuite.h>

#include "RandomBooleanNetwork.hpp"
#include "DifferentiationTree.hpp"
#include "DifferentiationTreeNode.hpp"
#include "ThresholdErgodicSetDifferentiationTree.hpp"
//This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

/* Next, we define the test class. */
class TestSearchingGeneActivationPatternsInThelperNetwork : public CxxTest::TestSuite
{
public:
	/*
	 * == Finding attractors and calculating the ATN from 'thelper' network ==
	 * In this example we search the attractors of the thelper network. We assert that the network has three single-point attractors and
	 * we calculate the ATN.
	 */

    void testThelper() throw (Exception)
    {
    	/* First of all we initialise Buddy. */
        bdd_init(10000,1000);
        try
        {
        	/* We instantiate a {{{RandomBooleanNetwork}}} object using a {{{ThresholdErgodicSetDifferentiationTree}}}
        	 * which is usually used for generate a {{{DifferentiationTree}}} object. In the constructor, the
        	 * {{{ThresholdErgodicSetDifferentiationTree}}} object initialise a {{{RandomBooleanNetwork}}} from the
        	 * 'thelper.net' network, and then it search its attractors.
        	 */
        	ThresholdErgodicSetDifferentiationTree TES_tree ("projects/SimoneRubinacci/networks_samples/thelper.net");

        	/* We assert that the number of attractors found is three. */
            TS_ASSERT_EQUALS(TES_tree.getBooleanNetwork()->getAttractorsNumber(),3u);

        	/* We assert that the attractors found are all single-point. */
            std::vector<unsigned> attractors_lengths = TES_tree.getBooleanNetwork()->getAttractorLength();
            for (unsigned i=0; i<attractors_lengths.size(); i++){
            	TS_ASSERT_EQUALS(attractors_lengths.at(i),1);
            }

            /* We export the stochastic matrix in a file. */
            //TES_tree.printStochasticMatrixAndAttractorLengthsToDatFile("networks_generated","matrix.dat");

            /* We also export the differentiation tree into a .gml file. */
            DifferentiationTree* diff_tree = TES_tree.getDifferentiationTree();
            //diff_tree->printDifferentiationTreeToGmlFile("networks_generated", "Thelper_diff_tree.gml");
            diff_tree->printDifferentiationTree();
        }
        catch (Exception& e)
        {
            TS_FAIL(e.GetMessage());
            bdd_done();
        }
        /* We release Buddy. */
        bdd_done();
    }
};



#endif /* TESTSEARCHINGGENEACTIVATIONPATTERNSINTHELPERNETWORK_HPP_ */
