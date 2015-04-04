#ifndef TESTRANDOMBOOLEANNETWORKGENERATOR_HPP_
#define TESTRANDOMBOOLEANNETWORK_HPP_

#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "SmartPointers.hpp"

#include <bdd.h>
#include "RandomBooleanNetwork.hpp"
#include "ThresholdErgodicSetDifferentiationTree.hpp"
#include "DifferentiationTree.hpp"
#include "RandomNumberGenerator.hpp"

#include "iostream"
#include <map>
#include <vector>
#include <set>
#include "FakePetscSetup.hpp"


/*
 * About randomness:
 * RandomNumberGenerator::Instance()->Reseed(0);
 * RandomNumberGenerator::Instance()->Reseed(time(NULL)); //time() returns the number of seconds since January 1970
 * RandomNumberGenerator::Instance()->Reseed(getpid()); //getpid() returns the system's process ID for the current program
 */

class TestRandomBooleanNetwork : public CxxTest::TestSuite
{
public:

    void testFindAttractorsFromFile() throw (Exception)
    {
        bdd_init(10000,1000);
        try
        {
            ThresholdErgodicSetDifferentiationTree TES_tree ("projects/SimoneRubinacci/networks_samples/fission_yeast.net");

            TS_ASSERT_EQUALS(TES_tree.getBooleanNetwork()->getAttractorsNumber(),13u);
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
        bdd_done();
    }

    void testRandomBooleanNetwork()
    {
        bdd_init(10000,1000);
        RandomBooleanNetwork *rbn1 = new RandomBooleanNetwork("projects/SimoneRubinacci/networks_samples/mammalian.cnet");

        delete rbn1;
        rbn1 = new RandomBooleanNetwork("projects/SimoneRubinacci/networks_samples/mammalian.cnet");
        rbn1->findAttractors();
        rbn1->printNetworkToNetFile("networks_generated", "mammalian_GEN_ATTRACTORS.net");
        rbn1->printNetworkToBooleanNetFile("networks_generated", "mammalian_booleannet_TXT_GEN.txt");
        rbn1->printNetworkToBoolNetFile("networks_generated", "mammalian_TXT_GEN.txt");
        rbn1->printGraphToGmlFile("networks_generated", "mammalian_graph_GEN.gml");
        rbn1->printGraphToDotFile("networks_generated", "mammalian_graph_GEN.dot");
        rbn1->printGraphToSifFile("networks_generated", "mammalian_graph_GEN.sif");
        bdd_done();
    }
    void testRbnGenerator()
    {
        bdd_init(10000,1000);
        RandomNumberGenerator::Instance()->Reseed(0);

        ThresholdErgodicSetDifferentiationTree TES_tree(11,2,true,1);

        TS_ASSERT_EQUALS(TES_tree.getBooleanNetwork()->getAttractorsNumber(), 4u);
        DifferentiationTree* Diff_tree = TES_tree.getDifferentiationTree();
        TS_ASSERT_EQUALS(Diff_tree->size(), 5u);
        TS_ASSERT_EQUALS(Diff_tree->getLeaves().size(), 4u);

        bdd_done();
    }

    void testSearchDifferentiationTree()
    {
        bool tree_isomorphism = false;
        unsigned attempt_number = 0;
        RandomNumberGenerator::Instance()->Reseed(0);

        bdd_init(200000,10000);

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
            tree_isomorphism = differentiation_tree->topologyTreeCompare(diff_tree);
            attempt_number++;
        } while (!tree_isomorphism);

        delete diff_tree;

        TS_ASSERT(tree_isomorphism);
        tes_tree->getBooleanNetwork()->printNetworkToBoolNetFile("networks_generated", "RandomNetwork_Match_GEN.txt");
        delete tes_tree;
        delete differentiation_tree;
        bdd_done();
    }

    void testThresholdErgodicSetDifferentiationTree()
    {
        /* Matrix taken from Graudenzi et al. http://dx.doi.org/10.1101/000927 */
    	std::vector<std::map<unsigned,double> > stoc_matrix(4);
        stoc_matrix.at(0).insert(std::pair<unsigned,double>(0,0.886));
        stoc_matrix.at(0).insert(std::pair<unsigned,double>(1,0.064));
        stoc_matrix.at(0).insert(std::pair<unsigned,double>(2,0.001));
        stoc_matrix.at(0).insert(std::pair<unsigned,double>(3,0.049));
        stoc_matrix.at(1).insert(std::pair<unsigned,double>(0,0.068));
        stoc_matrix.at(1).insert(std::pair<unsigned,double>(1,0.877));
        stoc_matrix.at(1).insert(std::pair<unsigned,double>(2,0.054));
        stoc_matrix.at(1).insert(std::pair<unsigned,double>(3,0.001));
        stoc_matrix.at(2).insert(std::pair<unsigned,double>(1,0.06));
        stoc_matrix.at(2).insert(std::pair<unsigned,double>(2,0.896));
        stoc_matrix.at(2).insert(std::pair<unsigned,double>(3,0.044));
        stoc_matrix.at(3).insert(std::pair<unsigned,double>(0,0.053));
        stoc_matrix.at(3).insert(std::pair<unsigned,double>(1,0.004));
        stoc_matrix.at(3).insert(std::pair<unsigned,double>(2,0.052));
        stoc_matrix.at(3).insert(std::pair<unsigned,double>(3,0.891));

        ThresholdErgodicSetDifferentiationTree TES_tree(stoc_matrix,
                        std::vector<unsigned>(4,1));

        DifferentiationTree* diff_tree = TES_tree.getDifferentiationTree();

        TS_ASSERT_EQUALS(diff_tree->getRoot()->getNumberOfChildren(), 2u);
        TS_ASSERT_EQUALS(diff_tree->getLeaves().size(), 4u);
        TS_ASSERT_EQUALS(diff_tree->size(), 10u);
        std::vector<std::set<unsigned> > level_nodes = diff_tree->getLevelNodes();

        delete diff_tree;

        TS_ASSERT_EQUALS(level_nodes.size(), 4u);
        TS_ASSERT_EQUALS(level_nodes.at(0).size(), 1u);
        TS_ASSERT_EQUALS(level_nodes.at(1).size(), 2u);
        TS_ASSERT_EQUALS(level_nodes.at(2).size(), 3u);
        TS_ASSERT_EQUALS(level_nodes.at(3).size(), 4u);
    }

    void testFixedTopologyRandomBooleanNetwork() throw (Exception)
    {
    	bdd_init(10000,1000);
        try
        {
            RandomBooleanNetwork rbn1("projects/SimoneRubinacci/networks_samples/mammalian_graph.gml", 0.5);
            rbn1.printGraphToGmlFile("networks_generated", "mammalian_graph_FixedTopologyRBN_GEN.gml");
            rbn1.printNetworkToNetFile("networks_generated", "mammalian_graph_FixedTopologyRBN_GEN.net");
        } catch (Exception& e)
        {
            TS_FAIL(e.GetMessage());
            bdd_done();
        }
        bdd_done();
    }

    void testBuddy()
    {

    	bdd x_0, x_1, x_2, x_3, x_4, nx_0, nx_1, nx_2, nx_3, nx_4, T=bddtrue, RT, I=bddtrue;
        bdd* replaceBdd;
        replaceBdd = new bdd[5];

        bddPair* renamepair;
        bddPair* replacePair;

        bdd_init(1000,100);
        bdd_setvarnum(2*5);
        x_0 = bdd_ithvar(0);
        x_1 = bdd_ithvar(2);
        x_2 = bdd_ithvar(4);
        x_3 = bdd_ithvar(6);
        x_4 = bdd_ithvar(8);
        nx_0 = bdd_ithvar(1);
        nx_1 = bdd_ithvar(3);
        nx_2 = bdd_ithvar(5);
        nx_3 = bdd_ithvar(7);
        nx_4 = bdd_ithvar(9);

        replaceBdd[0] = !(x_0 & x_2 & x_3 & x_4);
        replaceBdd[1] = !(x_0 & x_2 & x_3 & x_4);
        replaceBdd[2] = !(x_0 & x_2 & x_3 & x_4);
        replaceBdd[3] = !(x_0 & x_2 & x_3 & x_4);
        replaceBdd[4] = !(x_0 & x_2 & x_3 & x_4);

        int nvar[] = {0,2,4,6,8};
        int pvar[] = {1,3,5,7,9};

        renamepair = bdd_newpair();
        replacePair = bdd_newpair();
        bdd_setpairs(renamepair, pvar, nvar, 5); //used to replace operate
        bdd_setbddpairs(replacePair, nvar, replaceBdd, 5);


        replaceBdd[0] = bdd_veccompose(replaceBdd[0], replacePair);
        replaceBdd[1] = bdd_veccompose(replaceBdd[1], replacePair);
        replaceBdd[2] = bdd_veccompose(replaceBdd[2], replacePair);
        replaceBdd[3] = bdd_veccompose(replaceBdd[3], replacePair);
        replaceBdd[4] = bdd_veccompose(replaceBdd[4], replacePair);



        T &= bdd_apply(nx_0, replaceBdd[0], bddop_biimp);
        T &= bdd_apply(nx_1, replaceBdd[1], bddop_biimp);
        T &= bdd_apply(nx_2, replaceBdd[2], bddop_biimp);
        T &= bdd_apply(nx_3, replaceBdd[3], bddop_biimp);
        T &= bdd_apply(nx_4, replaceBdd[4], bddop_biimp);

        replaceBdd[0] = bdd_veccompose(replaceBdd[0], replacePair);
        replaceBdd[1] = bdd_veccompose(replaceBdd[1], replacePair);
        replaceBdd[2] = bdd_veccompose(replaceBdd[2], replacePair);
        replaceBdd[3] = bdd_veccompose(replaceBdd[3], replacePair);
        replaceBdd[4] = bdd_veccompose(replaceBdd[4], replacePair);

        replaceBdd[0] = bdd_veccompose(replaceBdd[0], replacePair);
        replaceBdd[1] = bdd_veccompose(replaceBdd[1], replacePair);
        replaceBdd[2] = bdd_veccompose(replaceBdd[2], replacePair);
        replaceBdd[3] = bdd_veccompose(replaceBdd[3], replacePair);
        replaceBdd[4] = bdd_veccompose(replaceBdd[4], replacePair);

        replaceBdd[0] = bdd_veccompose(replaceBdd[0], replacePair);
        replaceBdd[1] = bdd_veccompose(replaceBdd[1], replacePair);
        replaceBdd[2] = bdd_veccompose(replaceBdd[2], replacePair);
        replaceBdd[3] = bdd_veccompose(replaceBdd[3], replacePair);
        replaceBdd[4] = bdd_veccompose(replaceBdd[4], replacePair);

        bdd var_set = bdd_makeset(nvar, 5);
        bdd assignment = (x_0 & x_1 & x_2 & x_3 & x_4);
        assignment = bdd_appex(assignment, T, bddop_and, var_set);
        assignment = bdd_replace(assignment, renamepair);

        bdd_done();
    }

};
#endif /* TESTRANDOMBOOLEANNETWORK_HPP_ */
