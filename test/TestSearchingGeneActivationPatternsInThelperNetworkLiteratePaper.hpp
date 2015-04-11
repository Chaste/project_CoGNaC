#ifndef TESTSEARCHINGGENEACTIVATIONPATTERNSINTHELPERNETWORKLITERATEPAPER_HPP_
#define TESTSEARCHINGGENEACTIVATIONPATTERNSINTHELPERNETWORKLITERATEPAPER_HPP_

/*
 * = Searching gene activation patterns in T-helper cell differentiation =
 *
 * == Introduction ==
 *
 * EMPTYLINE
 *
 * In this class we show how to search gene activation patterns starting from a network,
 * calculate its Attractor Transition Network (Figure 3).
 *
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
	 *
	 * EMPTYLINE
	 *
	 * In this example we search the attractors of the thelper network. We test that
	 * the network has three single-point attractors and calculate the ATN.
	 */

    void testThelper() throw (Exception)
    {
    	/* First of all we initialise Buddy. */
        bdd_init(10000,1000);
        try
        {
        	/* We instantiate a {{{RandomBooleanNetwork}}} object using a {{{ThresholdErgodicSetDifferentiationTree}}}
        	 * which is used for generate a {{{DifferentiationTree}}} object. In the constructor, the
        	 * {{{ThresholdErgodicSetDifferentiationTree}}} object initialise a {{{RandomBooleanNetwork}}} from the
        	 * 'thelper.net' network, and then it search the attractors of the network.
        	 */
        	ThresholdErgodicSetDifferentiationTree TES_tree ("projects/CoGNaC/networks_samples/thelper.net");

        	/* We test that the number of attractors found is three. */
            TS_ASSERT_EQUALS(TES_tree.getBooleanNetwork()->getAttractorsNumber(),3u);

        	/* We test that the attractors found are all single-point. */
            std::vector<unsigned> attractors_lengths = TES_tree.getBooleanNetwork()->getAttractorLength();
            for (unsigned i=0; i<attractors_lengths.size(); i++){
            	TS_ASSERT_EQUALS(attractors_lengths.at(i),1);
            }

            /* We export the stochastic matrix in a file, where we can
             * visualise data shown in Figure 3 (Attractor Transition Network). */
            TES_tree.printStochasticMatrixAndAttractorLengthsToDatFile("networks_generated","stochastic_matrix_thelper.dat");
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



#endif /* TESTSEARCHINGGENEACTIVATIONPATTERNSINTHELPERNETWORKLITERATEPAPER_HPP_ */
