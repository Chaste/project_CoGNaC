#ifndef TESTRANDOMBOOLEANNETWORKGENERATOR_HPP_
#define TESTRANDOMBOOLEANNETWORK_HPP_

/*
 * = Testing the classes {{{RandomBooleanNetwork}}} and {{{ThresholdErgodicSetDifferentiationTree}}} =
 * This class is used to test that the classes {{{RandomBooleanNetwork}}} and {{{ThresholdErgodicSetDifferentiationTree}}} are implemented correctly.
 * It is also used to check that Buddy bdd tool is imported correctly.
 *
 * === Including header files ===
 * We begin by including the necessary header files.
*/

#include <cxxtest/TestSuite.h>

#include "iostream"
#include <map>
#include <vector>
#include <set>

#include "RandomBooleanNetwork.hpp"
#include "ThresholdErgodicSetDifferentiationTree.hpp"
#include "DifferentiationTree.hpp"

/* The next header includes Buddy library. */
#include <bdd.h>

/* The next header includes a number random generator.
 * It is possible to generate random numbers reseeding the generator.
 * Two ways are:
 * RandomNumberGenerator::Instance()->Reseed(time(NULL));
 * where time() returns the number of seconds since January 1970.
 * RandomNumberGenerator::Instance()->Reseed(getpid());
 * where getpid() returns the system's process ID for the current program.
 */
#include "RandomNumberGenerator.hpp"
//This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

class TestRandomBooleanNetwork : public CxxTest::TestSuite
{
public:

	/*
	 * == Finding attractors from 'fission yeast' network ==
	 * The network shows 13 fixed-point attractors.
	 *
	 */

    void testFindAttractorsFromFile() throw (Exception)
    {
    	/* First of all we initialise Buddy. */
        bdd_init(10000,1000);
        try
        {
        	/* We instantiate a {{{RandomBooleanNetwork}}} object using a {{{ThresholdErgodicSetDifferentiationTree}}}
        	 * which is usually used for generate a {{{DifferentiationTree}}} object. In the constructor, the
        	 * {{{ThresholdErgodicSetDifferentiationTree}}} object initialise a {{{RandomBooleanNetwork}}} from the
        	 * 'fission_yeast.cnet' file, and then it search the attractors.
        	 */
            ThresholdErgodicSetDifferentiationTree TES_tree ("projects/CoGNaC/networks_samples/fission_yeast.net");

            /*
             * We get the check if the attractor number is 13 as expected.
             */
            TS_ASSERT_EQUALS(TES_tree.getBooleanNetwork()->getAttractorsNumber(),13u);

            /*
             * We also test if the attractors are all fixed-point (set of one state).
             */
            std::vector<unsigned> attractors_lengths = TES_tree.getBooleanNetwork()->getAttractorLength();
            for (unsigned i=0; i< TES_tree.getBooleanNetwork()->getAttractorsNumber(); i++)
            {
                TS_ASSERT_EQUALS(attractors_lengths.at(i), 1u);
            }
        }
        catch (Exception& e)
        {
            TS_FAIL(e.GetMessage());
            bdd_done();
        }
        /* We release Buddy. */
        bdd_done();
    }

	/*
	 * == Testing input and output in {{{RandomBooleanNetwork}}} object ==
	 * We generate a {{{RandomBooleanNetwork}}} from a file named 'mammalian.cnet', and then we search the attractors of the network.
	 * We finally export the network and its directed graph in several file formats.
	 */

    void testRandomBooleanNetwork()
    {
    	/* First of all we initialise Buddy.
    	 */
        bdd_init(10000,1000);

        /* We generate a network from 'mammalian.cnet' file.
         */
        RandomBooleanNetwork *rbn1 = new RandomBooleanNetwork("projects/CoGNaC/networks_samples/mammalian.cnet");
        /* We search the attractors of the network.
         */
        rbn1->findAttractors();
        /* We save the network in several distinct file formats:
         * .net, BooleanNet and BoolNet. We also save its graph
         * in .gml, .dot and .sif file format.
         */
        rbn1->printNetworkToNetFile("networks_generated", "mammalian_GEN_ATTRACTORS.net");
        rbn1->printNetworkToBooleanNetFile("networks_generated", "mammalian_booleannet_TXT_GEN.txt");
        rbn1->printNetworkToBoolNetFile("networks_generated", "mammalian_TXT_GEN.txt");
        rbn1->printGraphToGmlFile("networks_generated", "mammalian_graph_GEN.gml");
        rbn1->printGraphToDotFile("networks_generated", "mammalian_graph_GEN.dot");
        rbn1->printGraphToSifFile("networks_generated", "mammalian_graph_GEN.sif");
        /* We release Buddy. */
        bdd_done();
    }

	/*
	 * == Testing properties of a fixed {{{RandomBooleanNetwork}}} ==
	 * We generate a {{{RandomBooleanNetwork}}} using a fixed seed of the {{{RandomNumberGenerator}}} and
	 * then we test the properties of the network.
	 */

    void testRbnGenerator()
    {
    	/* First of all we initialise Buddy. */
        bdd_init(10000,1000);
        /* We instantiate the {{{RandomNumberGenerator}}} into a specific seed.
         * Usually it is chosen at random, using time(NULL) or getpid() functions.
         */
        RandomNumberGenerator::Instance()->Reseed(0);
        /* We generate a random network with 11 nodes, K=2, scale-free functions and completely canalyzing functions.
		 */
        ThresholdErgodicSetDifferentiationTree TES_tree(11,2,true,1);

        /* We test that the generated  network has four attractors. */
        TS_ASSERT_EQUALS(TES_tree.getBooleanNetwork()->getAttractorsNumber(), 4u);

        /* We also test that the generated differentiation tree has five nodes
         * with a root and four leaves.
         */
        DifferentiationTree* Diff_tree = TES_tree.getDifferentiationTree();
        Diff_tree->printDifferentiationTree();
        TS_ASSERT_EQUALS(Diff_tree->size(), 5u);
        TS_ASSERT_EQUALS(Diff_tree->getLeaves().size(), 4u);

        /* We release Buddy. */
        bdd_done();
    }
    /*
	 * == Searching a network which shows a particular {{{DifferentiationTree}}} ==
	 * We randomly generate {{{RandomBooleanNetwork}}} objects until we obtain a fixed
	 * {{{Differentiation}}} using the function {{{topologyTreeCompare}}}.
	 */
    void testSearchDifferentiationTree()
    {
    	/* First of all we initialise Buddy. */
        bdd_init(200000,10000);

        bool tree_isomorphism = false;
        unsigned attempt_number = 0;
        RandomNumberGenerator::Instance()->Reseed(0);

        DifferentiationTree* diff_tree = new DifferentiationTree();
        diff_tree->addNewChild(0);
        diff_tree->addNewChild(0);

        ThresholdErgodicSetDifferentiationTree* tes_tree;
        DifferentiationTree* differentiation_tree;
        do{
        	if (attempt_number > 0)
            {
        		delete tes_tree;
        		delete differentiation_tree;
            }
        	tes_tree = new ThresholdErgodicSetDifferentiationTree(10,2,true,.5);
            differentiation_tree = tes_tree->getDifferentiationTree();
            differentiation_tree->printDifferentiationTree();
            tree_isomorphism = differentiation_tree->topologyTreeCompare(diff_tree);
            attempt_number++;
        } while (!tree_isomorphism);

        delete diff_tree;

        TS_ASSERT(tree_isomorphism);
        tes_tree->getBooleanNetwork()->printNetworkToBoolNetFile("networks_generated", "RandomNetwork_Match_GEN.txt");
        delete tes_tree;
        delete differentiation_tree;
        /* We release Buddy. */
        bdd_done();
    }

    /*
	 * == Testing the diffentiation tree from a fixed network ==
	 * We generate a {{{RandomBooleanNetwork}}} from a file, and we generate its
	 * {{{DifferentiationTree}}}. Then, we check its topological properties.
	 */
    void testThresholdErgodicSetDifferentiationTree()
    {
        /* Matrix taken from Graudenzi et al. http://dx.doi.org/10.1101/000927 */
    	ThresholdErgodicSetDifferentiationTree TES_tree("projects/CoGNaC/networks_samples/graudenzi_matrix.dat");
    	/* We obtain its differentiation tree. */
    	DifferentiationTree* diff_tree = TES_tree.getDifferentiationTree();

    	/* We assert properties on the number of children for each node. */
        TS_ASSERT_EQUALS(diff_tree->getRoot()->getNumberOfChildren(), 2u);
        TS_ASSERT_EQUALS(diff_tree->getLeaves().size(), 4u);
        TS_ASSERT_EQUALS(diff_tree->size(), 7u);
        std::vector<std::set<unsigned> > level_nodes = diff_tree->getLevelNodes();

        /* We export the differentiation tree into a .gml file. */
        //diff_tree->printDifferentiationTreeToGmlFile("networks_generated","graudenzi_cell_tree.gml");
        diff_tree->printDifferentiationTree();
        delete diff_tree;

    	/* We assert properties for each level in the tree. */
        TS_ASSERT_EQUALS(level_nodes.size(), 4u);
        TS_ASSERT_EQUALS(level_nodes.at(0).size(), 1u);
        TS_ASSERT_EQUALS(level_nodes.at(1).size(), 2u);
        TS_ASSERT_EQUALS(level_nodes.at(2).size(), 2u);
        TS_ASSERT_EQUALS(level_nodes.at(3).size(), 2u);
    }

    /*
	 * == Generating {{{RandomBooleanNetwork}}} with fixed topology ==
	 * We randomly generate {{{RandomBooleanNetwork}}} objects starting from a
	 * fixed graph and then we save it to a .gml file in order to test the
	 * equality.
	 */
    void testFixedTopologyRandomBooleanNetwork() throw (Exception)
    {
    	/* First of all we initialise Buddy. */
    	bdd_init(10000,1000);
        try
        {
        	/* We generate a {{{RandomBooleanNetwork}}} starting from a graph. */
            RandomBooleanNetwork rbn1("projects/CoGNaC/networks_samples/mammalian_graph.gml", 0.5);
            /* And then we save it into a .gml file, in order to check its topological properties using Cytoscape */
            rbn1.printGraphToGmlFile("networks_generated", "mammalian_graph_FixedTopologyRBN_GEN.gml");
        } catch (Exception& e)
        {
            TS_FAIL(e.GetMessage());
            bdd_done();
        }
        /* We release Buddy. */
        bdd_done();
    }

    /*
	 * == Testing Buddy ==
	 * We test that BuDDy works properly, using an example contained in BuDDy's guide.
	 */

    void testBuddy()
    {
    	/* First of all we initialise Buddy. */
        bdd_init(1000,100);

        /* We use three variables. */
        bdd x,y,z;

        bdd_setvarnum(5);

        /* We impose z be x AND y. */
        x = bdd_ithvar(0);
        y = bdd_ithvar(1);
        z = x & y;

        /* We release Buddy. */
        bdd_done();
    }

};
#endif /* TESTRANDOMBOOLEANNETWORK_HPP_ */
