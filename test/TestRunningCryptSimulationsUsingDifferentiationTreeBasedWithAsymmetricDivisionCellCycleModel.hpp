#ifndef TESTRUNNINGCRYPTSIMULATIONSUSINGDIFFERENTIATIONTREEBASEDCELLCYCLEMODEL_HPP_
#define TESTRUNNINGCRYPTSIMULATIONSUSINGDIFFERENTIATIONTREEBASEDCELLCYCLEMODEL_HPP_

#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "SmartPointers.hpp"
#include "AbstractCellBasedTestSuite.hpp"

#include <vector>
#include <set>
#include "float.h"

#include "ThresholdErgodicSetDifferentiationTree.hpp"
#include "DifferentiationTree.hpp"
#include "DifferentiationTreeNode.hpp"
#include "DifferentiationTreeBasedWithAsymmetricDivisionCellCycleModel.hpp"
#include "CellDifferentiationTypeWriter.hpp"

#include "CellsGenerator.hpp"
#include "StemCellProliferativeType.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "OffLatticeSimulation.hpp"
#include "VoronoiDataWriter.hpp"
#include "SloughingCellKiller.hpp"
#include "GeneralisedLinearSpringForce.hpp"

#include "FakePetscSetup.hpp"


class BoundaryConditionWidthAndBottom : public AbstractCellPopulationBoundaryCondition<2>
{
private:

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellPopulationBoundaryCondition<2> >(*this);
    }

public:
    BoundaryConditionWidthAndBottom(AbstractCellPopulation<2>* pCellPopulation)
        : AbstractCellPopulationBoundaryCondition<2>(pCellPopulation)
    {
    }

    void ImposeBoundaryCondition(const std::map<Node<2>*, c_vector<double, 2> >& rOldLocations)
    {
        for (AbstractCellPopulation<2>::Iterator cell_iter = this->mpCellPopulation->Begin();
             cell_iter != this->mpCellPopulation->End();
             ++cell_iter)
        {
            unsigned node_index = this->mpCellPopulation->GetLocationIndexUsingCell(*cell_iter);
            Node<2>* p_node = this->mpCellPopulation->GetNode(node_index);
            double x_coordinate = p_node->rGetLocation()[0];
            double y_coordinate = p_node->rGetLocation()[1];

            if (x_coordinate > 20.0)
            {
                p_node->rGetModifiableLocation()[0] = 20.0;
            }
            else if (x_coordinate < 0.0)
            {
                p_node->rGetModifiableLocation()[0] = 0.0;
            }

            if (y_coordinate < 0.0)
            {
            	p_node->rGetModifiableLocation()[1] = 0.0;
            }
        }
    }

    bool VerifyBoundaryCondition()
    {
        bool condition_satisfied = true;

        for (AbstractCellPopulation<2>::Iterator cell_iter = this->mpCellPopulation->Begin();
             cell_iter != this->mpCellPopulation->End();
             ++cell_iter)
        {
            c_vector<double, 2> cell_location = this->mpCellPopulation->GetLocationOfCellCentre(*cell_iter);
            double x_coordinate = cell_location(0);
            double y_coordinate = cell_location(1);

            if ((x_coordinate < 0.0) || (x_coordinate > 20.0))
            {
                condition_satisfied = false;
                break;
            }
            if ((y_coordinate < 0.0))
            {
                condition_satisfied = false;
                break;
            }
        }
        return condition_satisfied;
    }

    void OutputCellPopulationBoundaryConditionParameters(out_stream& rParamsFile)
    {
        AbstractCellPopulationBoundaryCondition<2>::OutputCellPopulationBoundaryConditionParameters(rParamsFile);
    }
};

#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(BoundaryConditionWidthAndBottom)
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(BoundaryConditionWidthAndBottom)

namespace boost
{
    namespace serialization
    {
        template<class Archive>
        inline void save_construct_data(
            Archive & ar, const BoundaryConditionWidthAndBottom * t, const BOOST_PFTO unsigned int file_version)
        {
            const AbstractCellPopulation<2>* const p_cell_population = t->GetCellPopulation();
            ar << p_cell_population;
        }

        template<class Archive>
        inline void load_construct_data(
            Archive & ar, BoundaryConditionWidthAndBottom * t, const unsigned int file_version)
        {
            AbstractCellPopulation<2>* p_cell_population;
            ar >> p_cell_population;

            ::new(t)BoundaryConditionWidthAndBottom(p_cell_population);
        }
    }
}

class TestRunningCryptSimulationsUsingDifferentiationTreeBasedWithAsymmetricDivisionCellCycleModel : public AbstractCellBasedTestSuite
{
public:

    void TestThresholdErgodicSetDifferentiationTreeFromAnExample()
    {
    	/*
    	std::vector<std::map<unsigned,double> > stoc_matrix(3);
        stoc_matrix.at(0).insert(std::pair<unsigned,double>(0,0.9746));
        stoc_matrix.at(0).insert(std::pair<unsigned,double>(1,0.0246));
        stoc_matrix.at(0).insert(std::pair<unsigned,double>(2,0.0007));
        stoc_matrix.at(1).insert(std::pair<unsigned,double>(0,0.4464));
        stoc_matrix.at(1).insert(std::pair<unsigned,double>(1,0.5457));
        stoc_matrix.at(1).insert(std::pair<unsigned,double>(2,0.0079));
        stoc_matrix.at(2).insert(std::pair<unsigned,double>(0,0.4050));
        stoc_matrix.at(2).insert(std::pair<unsigned,double>(1,0.1350));
        stoc_matrix.at(2).insert(std::pair<unsigned,double>(2,0.4600));

        std::vector<unsigned> attractor_lengths;
        attractor_lengths.push_back(113);
        attractor_lengths.push_back(29);
        attractor_lengths.push_back(4);

        ThresholdErgodicSetDifferentiationTree TES_tree(stoc_matrix,
        		attractor_lengths);
        */
    	ThresholdErgodicSetDifferentiationTree TES_tree("projects/CoGNaC/networks_samples/fig4_atn.dat");

        DifferentiationTree* diff_tree = TES_tree.getDifferentiationTree();

        //Assert topological properties of the tree.
        TS_ASSERT_EQUALS(diff_tree->getRoot()->getNumberOfChildren(), 3u);
		TS_ASSERT_EQUALS(diff_tree->getLeaves().size(), 3u);
		TS_ASSERT_EQUALS(diff_tree->size(), 4u);
		std::vector<std::set<unsigned> > level_nodes = diff_tree->getLevelNodes();

		TS_ASSERT_EQUALS(level_nodes.size(), 2u);
		TS_ASSERT_EQUALS(level_nodes.at(0).size(), 1u);
		TS_ASSERT_EQUALS(level_nodes.at(1).size(), 3u);


		std::vector<double> stochatic_probabilities;
		stochatic_probabilities.push_back(0.9459740);
		stochatic_probabilities.push_back(0.0518657);
		stochatic_probabilities.push_back(0.00216022);
		std::vector<double> stochatic_probabilities_root = diff_tree->getRoot()->getStationaryDistribution();
		std::vector<double>::iterator iterator;

		TS_ASSERT_EQUALS(stochatic_probabilities.size(), stochatic_probabilities_root.size());
		for (unsigned i = 0; i<stochatic_probabilities.size(); ++i)
		{
			TS_ASSERT_DELTA(stochatic_probabilities.at(i), stochatic_probabilities_root.at(i), 1e-6);
		}

        delete diff_tree;
    }

    void TestBoundaryConditionWidthAndBottom() throw(Exception)
    {
        HoneycombMeshGenerator generator(25, 4);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumNodes());

        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        BoundaryConditionWidthAndBottom bc(&cell_population);

        bool population_satisfies_bc = bc.VerifyBoundaryCondition();
        TS_ASSERT_EQUALS(population_satisfies_bc, false);

        std::map<Node<2>*, c_vector<double, 2> > old_node_locations;
        for (AbstractMesh<2,2>::NodeIterator node_iter = p_mesh->GetNodeIteratorBegin();
                node_iter != p_mesh->GetNodeIteratorEnd();
                ++node_iter)
        {
            old_node_locations[&(*node_iter)] = node_iter->rGetLocation();
        }

        bc.ImposeBoundaryCondition(old_node_locations);

        population_satisfies_bc = bc.VerifyBoundaryCondition();
        TS_ASSERT_EQUALS(population_satisfies_bc, true);
    }

    void noTestSimulationWithBoundaries()
	{
    	//RandomNumberGenerator::Instance()->Reseed(time(NULL));
    	/*
		std::vector<std::map<unsigned,double> > stoc_matrix(3);
		stoc_matrix.at(0).insert(std::pair<unsigned,double>(0,0.9746));
		stoc_matrix.at(0).insert(std::pair<unsigned,double>(1,0.0246));
		stoc_matrix.at(0).insert(std::pair<unsigned,double>(2,0.0007));
		stoc_matrix.at(1).insert(std::pair<unsigned,double>(0,0.4464));
		stoc_matrix.at(1).insert(std::pair<unsigned,double>(1,0.5457));
		stoc_matrix.at(1).insert(std::pair<unsigned,double>(2,0.0079));
		stoc_matrix.at(2).insert(std::pair<unsigned,double>(0,0.4050));
		stoc_matrix.at(2).insert(std::pair<unsigned,double>(1,0.1350));
		stoc_matrix.at(2).insert(std::pair<unsigned,double>(2,0.4600));

		std::vector<unsigned> attractor_lengths;
		attractor_lengths.push_back(113);
		attractor_lengths.push_back(29);
		attractor_lengths.push_back(4);

		ThresholdErgodicSetDifferentiationTree TES_tree(stoc_matrix,
						attractor_lengths);
		*/
    	ThresholdErgodicSetDifferentiationTree TES_tree("projects/CoGNaC/networks_samples/fig4_atn.dat");

		DifferentiationTree* diff_tree = TES_tree.getDifferentiationTree();
		diff_tree->printDifferentiationTreeToGmlFile("networks_generated","fig4_differentiation_tree.gml");

		markLessProbableWithBlackColour(diff_tree);
		diff_tree->normaliseLength(8.0);
		HoneycombMeshGenerator generator(20, 20, 4);
		MutableMesh<2,2>* p_mesh = generator.GetMesh();
		std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

		std::vector<CellPtr> cells;
		MAKE_PTR(StemCellProliferativeType, p_stem_type);
		CellsGenerator<DifferentiationTreeBasedWithAsymmetricDivisionCellCycleModel, 2> cells_generator;
		cells_generator.GenerateBasicRandom(cells, location_indices.size(), p_stem_type);

		for (unsigned i=0; i<cells.size(); i++)
		{
			DifferentiationTreeBasedWithAsymmetricDivisionCellCycleModel* p_cell_cycle_model =
								new DifferentiationTreeBasedWithAsymmetricDivisionCellCycleModel(diff_tree,0);
			CellPtr p_cell = cells.at(i);
			p_cell->SetCellCycleModel(p_cell_cycle_model);
			p_cell->SetBirthTime(-diff_tree->getRoot()->getCellCycleLength()*RandomNumberGenerator::Instance()->ranf());
			p_cell->InitialiseCellCycleModel();
		}

		MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, location_indices);
		cell_population.AddCellWriter<CellDifferentiationTypeWriter>();
		cell_population.SetWriteVtkAsPoints(false);
		cell_population.AddPopulationWriter<VoronoiDataWriter>();

        MAKE_PTR_ARGS(BoundaryConditionWidthAndBottom, p_bc, (&cell_population));

		OffLatticeSimulation<2> simulator(cell_population);
		simulator.SetOutputDirectory("BoundariesMeshsDifferentiationTreeBasedWithAsymmetricDivisionCellCycleModel");
		simulator.SetEndTime(16.0);

		MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
		p_linear_force->SetMeinekeSpringStiffness(30.0);
		p_linear_force->SetMeinekeSpringGrowthDuration(0.0);
		//p_linear_force->SetCutOffLength(30);

		simulator.AddForce(p_linear_force);

        simulator.AddCellPopulationBoundaryCondition(p_bc);

        double tissue_height = 17.00;
        MAKE_PTR_ARGS(SloughingCellKiller<2>, p_killer, (&cell_population, tissue_height));
        simulator.AddCellKiller(p_killer);

		simulator.Solve();
		delete diff_tree;

	}

    void markLessProbableWithBlackColour(DifferentiationTree* diff_tree)
    {
    	diff_tree->setColour(0, 0.0);
    	double minimum = DBL_MAX;
    	unsigned minimum_index = 0;
    	double maximum = 0.0;
    	unsigned maximum_index = 0;

    	for (unsigned i=1;i<diff_tree->size(); i++)
    	{
			if (diff_tree->getNode(i)->getCellCycleLength() < minimum)
			{
				minimum = diff_tree->getNode(i)->getCellCycleLength();
				minimum_index = i;
			}
			if (diff_tree->getNode(i)->getCellCycleLength() > maximum)
			{
				maximum = diff_tree->getNode(i)->getCellCycleLength();
				maximum_index = i;
			}
    	}
    	if (minimum_index > 0)
    	{
    		for (unsigned i=1;i<diff_tree->size(); i++)
			{
    			if (i == minimum_index)
    			{
    				diff_tree->setColour(i, 4.0);
    			}
    			else if (i == maximum_index)
    			{
    				diff_tree->setColour(i, 2.0);
    			}
    			else diff_tree->setColour(i, 1.0);
			}
    	}
    }

    void noTestGenerateFig2()
    {
        DifferentiationTree* tree = new DifferentiationTree(4);
        tree->addNewChild(0,4);
        tree->addNewChild(0,4);

        tree->colorise();


        HoneycombMeshGenerator generator(1, 1, 4);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        std::vector<CellPtr> cells;
        MAKE_PTR(StemCellProliferativeType, p_stem_type);
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        CellsGenerator<DifferentiationTreeBasedWithAsymmetricDivisionCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, location_indices.size(), p_stem_type);

        for (unsigned i=0; i<cells.size(); i++)
        {
        	DifferentiationTreeBasedWithAsymmetricDivisionCellCycleModel* p_cell_cycle_model =
                                new DifferentiationTreeBasedWithAsymmetricDivisionCellCycleModel(tree,0);
            CellPtr p_cell = cells.at(i);
            p_cell->SetCellCycleModel(p_cell_cycle_model);
            p_cell->InitialiseCellCycleModel();
        }

        MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, location_indices);
        cell_population.AddCellWriter<CellDifferentiationTypeWriter>();
        cell_population.SetWriteVtkAsPoints(false);
        cell_population.AddPopulationWriter<VoronoiDataWriter>();

        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestFig2");
        simulator.SetEndTime(16.0);

        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        simulator.AddForce(p_linear_force);

        simulator.Solve();
        delete tree;
    }
};
#endif /* TESTRUNNINGCRYPTSIMULATIONSUSINGDIFFERENTIATIONTREEBASEDCELLCYCLEMODEL_HPP_ */
