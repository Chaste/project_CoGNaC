#ifndef TESTCRYPTBASEDSIMULATIONSUSINGDIFFERENTIATIONTREECELLCYCLEMODEL_HPP_
#define TESTCRYPTBASEDSIMULATIONSUSINGDIFFERENTIATIONTREECELLCYCLEMODEL_HPP_

#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "SmartPointers.hpp"
#include "AbstractCellBasedTestSuite.hpp"

#include "DifferentiationTree.hpp"
#include "DifferentiationTreeNode.hpp"
#include "DifferentiationTreeBasedCellCycleModel.hpp"

#include "CryptCellsGenerator.hpp"
#include "CylindricalHoneycombMeshGenerator.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "CryptSimulation2d.hpp"
#include "SloughingCellKiller.hpp"
#include "VoronoiDataWriter.hpp"
#include "CellDifferentiationTypeWriter.hpp"

#include "HoneycombVertexMeshGenerator.hpp"
#include "CylindricalHoneycombVertexMeshGenerator.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "NagaiHondaForce.hpp"
#include "SimpleTargetAreaModifier.hpp"

#include "CellMutationStatesCountWriter.hpp"

#include "FakePetscSetup.hpp"
#include "RandomNumberGenerator.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"

class TestCryptBasedSimulationsUsingDifferentiationTreeCellCycleModel : public AbstractCellBasedTestSuite
{
public:
    void Test2DCryptMeshBasedSimulationWithMutationUsingDifferentiationTreeCellCycleModel()
    {
        /*** TREE ***/
        DifferentiationTree* tree = new DifferentiationTree(4);
        tree->addNewChild(0,3);
        tree->addNewChild(0,3);
        tree->addNewChild(2,3);
        tree->addNewChild(2,3);
        tree->addNewChild(3,2);
        tree->addNewChild(3,1);
        tree->addNewChild(4,2);


        //tree->normaliseLength(.8);
        tree->colorise();
        /*** END TREE ***/

        CylindricalHoneycombMeshGenerator generator(6, 9, 2);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();

        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        std::vector<CellPtr> cells;
        CryptCellsGenerator<DifferentiationTreeBasedCellCycleModel> cells_generator;
        cells_generator.Generate(cells, p_mesh, location_indices, true);

        boost::shared_ptr<AbstractCellProperty> p_state(CellPropertyRegistry::Instance()->Get<ApcTwoHitCellMutationState>());

        for (unsigned i=0; i<cells.size(); i++)
        {
            DifferentiationTreeBasedCellCycleModel* p_cell_cycle_model =
                    new DifferentiationTreeBasedCellCycleModel(tree,0);
            CellPtr p_cell = cells.at(i);
            p_cell->SetCellCycleModel(p_cell_cycle_model);
            p_cell->InitialiseCellCycleModel();
        }

        MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, location_indices);
        cell_population.AddCellWriter<CellDifferentiationTypeWriter>();
        cell_population.SetWriteVtkAsPoints(false);
        cell_population.AddPopulationWriter<VoronoiDataWriter>();
        cell_population.AddCellPopulationCountWriter<CellMutationStatesCountWriter>();

        double crypt_height = 8.0;
        CryptSimulation2d simulator(cell_population);
        simulator.SetOutputDirectory("CryptMeshWithMutationsDifferentiationTreeCellCycleModel");
        simulator.SetSamplingTimestepMultiple(12);
        simulator.SetEndTime(10);

        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        simulator.AddForce(p_linear_force);
        MAKE_PTR_ARGS(SloughingCellKiller<2>, p_killer, (&cell_population, crypt_height));
        simulator.AddCellKiller(p_killer);

        simulator.Solve();

        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            unsigned node_index = cell_population.GetLocationIndexUsingCell(*cell_iter);

            if (node_index == 74) // Chosen from looking at the results from steady state
            {
                cell_iter->SetMutationState(p_state);
            }
        }

       double normal_damping_constant = cell_population.GetDampingConstantNormal();
       cell_population.SetDampingConstantMutant(10*normal_damping_constant);

       simulator.SetEndTime(20);

       simulator.Solve();
       delete tree;
    }

    void Test2DCryptMeshBasedSimulationWithMorphologyUsingDifferentiationTreeCellCycleModel()
    {
        /*** TREE ***/
        DifferentiationTree* tree = new DifferentiationTree(4);
        tree->addNewChild(0,3);
        tree->addNewChild(0,3);
        tree->addNewChild(2,2);
        tree->addNewChild(2,2);
        tree->addNewChild(3,1);
        tree->addNewChild(3,1);
        tree->addNewChild(4,1);

        //tree->normaliseLength(.8);
        tree->colorise();
        /*** END TREE ***/

        CylindricalHoneycombMeshGenerator generator(20, 22, 2);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();

        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        std::vector<CellPtr> cells;
        CryptCellsGenerator<DifferentiationTreeBasedCellCycleModel> cells_generator;
        cells_generator.Generate(cells, p_mesh, location_indices, true);

        boost::shared_ptr<AbstractCellProperty> p_state(CellPropertyRegistry::Instance()->Get<ApcTwoHitCellMutationState>());

        for (unsigned i=0; i<cells.size(); i++)
        {
            DifferentiationTreeBasedCellCycleModel* p_cell_cycle_model;
            if (i < 60)
            {
                p_cell_cycle_model = new DifferentiationTreeBasedCellCycleModel(tree,1);
            }
            else if(i <120)
            {
                p_cell_cycle_model = new DifferentiationTreeBasedCellCycleModel(tree,0);
            }
            else if(i < 360)
            {
                p_cell_cycle_model = new DifferentiationTreeBasedCellCycleModel(tree,2);
            }
            else if(i < 400)
            {
                double random_number = RandomNumberGenerator::Instance()->ranf();
                if (random_number > .5)
                {
                    p_cell_cycle_model = new DifferentiationTreeBasedCellCycleModel(tree,3);
                }
                else
                {
                    p_cell_cycle_model = new DifferentiationTreeBasedCellCycleModel(tree,4);
                }
            }
            else
            {
                double random_number = RandomNumberGenerator::Instance()->ranf();
                if (random_number<.333333334)
                {
                    p_cell_cycle_model = new DifferentiationTreeBasedCellCycleModel(tree,5);
                }
                else if (random_number < .666666667)
                {
                    p_cell_cycle_model = new DifferentiationTreeBasedCellCycleModel(tree,6);
                }
                else
                {
                    p_cell_cycle_model = new DifferentiationTreeBasedCellCycleModel(tree,7);
                }
            }
            CellPtr p_cell = cells.at(i);
            p_cell->SetCellCycleModel(p_cell_cycle_model);
            p_cell->InitialiseCellCycleModel();
        }

        MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, location_indices);
        cell_population.AddCellWriter<CellDifferentiationTypeWriter>();
        cell_population.SetWriteVtkAsPoints(false);
        cell_population.AddPopulationWriter<VoronoiDataWriter>();
        cell_population.AddCellPopulationCountWriter<CellMutationStatesCountWriter>();

        double crypt_height = 19;
        CryptSimulation2d simulator(cell_population);
        simulator.SetOutputDirectory("CryptMeshMorphologyDifferentiationTreeCellCycleModel");
        simulator.SetSamplingTimestepMultiple(12);
        simulator.SetEndTime(10);

        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        simulator.AddForce(p_linear_force);
        MAKE_PTR_ARGS(SloughingCellKiller<2>, p_killer, (&cell_population, crypt_height));
        simulator.AddCellKiller(p_killer);

        simulator.Solve();
        delete tree;
    }

    void Test2DCryptMeshBasedSimulationUsingDifferentiationTreeCellCycleModel()
    {
        /*** TREE ***/
        DifferentiationTree* tree = new DifferentiationTree(5);
        tree->addNewChild(0,2);
        tree->addNewChild(0,2);
        tree->addNewChild(2,2);
        tree->addNewChild(2,2);
        tree->addNewChild(3,2);
        tree->addNewChild(3,2);
        tree->addNewChild(4,2);

        //tree->normaliseLength(8);
        tree->colorise();
        /*** END TREE ***/

        CylindricalHoneycombMeshGenerator generator(6, 9, 2);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        std::vector<CellPtr> cells;
        CryptCellsGenerator<DifferentiationTreeBasedCellCycleModel> cells_generator;
        cells_generator.Generate(cells, p_mesh, location_indices, true);

        for (unsigned i=0; i<cells.size(); i++)
        {
            DifferentiationTreeBasedCellCycleModel* p_cell_cycle_model =
                    new DifferentiationTreeBasedCellCycleModel(tree,0);
            CellPtr p_cell = cells.at(i);
            p_cell->SetCellCycleModel(p_cell_cycle_model);
            p_cell->InitialiseCellCycleModel();
        }
        MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, location_indices);

        cell_population.AddCellWriter<CellDifferentiationTypeWriter>();
        cell_population.SetWriteVtkAsPoints(false);
        cell_population.AddPopulationWriter<VoronoiDataWriter>();

        CryptSimulation2d simulator(cell_population);
        simulator.SetOutputDirectory("CryptMeshDifferentiationTreeCellCycleModel");
        simulator.SetEndTime(8);
        simulator.SetSamplingTimestepMultiple(12);

        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        simulator.AddForce(p_linear_force);

        double crypt_height = 8.0;
        MAKE_PTR_ARGS(SloughingCellKiller<2>, p_killer, (&cell_population, crypt_height));
        simulator.AddCellKiller(p_killer);

        simulator.Solve();
        delete tree;
    }

    void Test2DCryptNodeBasedSimulationUsingDifferentiationTreeCellCycleModel()
    {
		/*** TREE ***/
		DifferentiationTree* tree = new DifferentiationTree(2);
		tree->addNewChild(0,1);
		tree->addNewChild(0,2);
		tree->addNewChild(2,2);
		tree->addNewChild(2,1);
		tree->addNewChild(3,2);
		tree->addNewChild(3,1);
		tree->addNewChild(4,1);
		//tree->normaliseLength(.8);
		tree->colorise();
		/*** END TREE ***/

		CylindricalHoneycombVertexMeshGenerator generator(6, 9);
		Cylindrical2dVertexMesh* p_mesh = generator.GetCylindricalMesh();
		std::vector<CellPtr> cells;
		CryptCellsGenerator<DifferentiationTreeBasedCellCycleModel> cells_generator;
		cells_generator.Generate(cells, p_mesh, std::vector<unsigned>(), true, 1.0, 2.0, 3.0, 4.0);
		for (unsigned i=0; i<cells.size(); i++)
		{
			DifferentiationTreeBasedCellCycleModel* p_cell_cycle_model =
			new DifferentiationTreeBasedCellCycleModel(tree,0);
			CellPtr p_cell = cells.at(i);
			p_cell->SetCellCycleModel(p_cell_cycle_model);
			p_cell->InitialiseCellCycleModel();
		}
		VertexBasedCellPopulation<2> crypt(*p_mesh, cells);

		CryptSimulation2d simulator(crypt);
		simulator.SetOutputDirectory("CryptVertexDifferentiationTreeCellCycleModel");
		simulator.SetEndTime(.5);

		MAKE_PTR(NagaiHondaForce<2>, p_force);
		simulator.AddForce(p_force);
		MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
		simulator.AddSimulationModifier(p_growth_modifier);
		double crypt_length = 6.0;
		MAKE_PTR_ARGS(SloughingCellKiller<2>, p_killer, (&crypt, crypt_length));
		simulator.AddCellKiller(p_killer);
		simulator.Solve();
		delete tree;
    }

};
#endif /* TESTCRYPTBASEDSIMULATIONSUSINGDIFFERENTIATIONTREECELLCYCLEMODEL_HPP_ */
