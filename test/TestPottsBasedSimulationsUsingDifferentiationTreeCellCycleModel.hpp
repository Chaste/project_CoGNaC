#ifndef TESTPOTTSBASEDSIMULATIONSUSINGDIFFERENTIATIONTREECELLCYCLEMODEL_HPP_
#define TESTPOTTSBASEDSIMULATIONSUSINGDIFFERENTIATIONTREECELLCYCLEMODEL_HPP_

#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "SmartPointers.hpp"
#include "AbstractCellBasedTestSuite.hpp"

#include "DifferentiationTree.hpp"
#include "DifferentiationTreeNode.hpp"
#include "DifferentiationTreeBasedCellCycleModel.hpp"

#include "CellsGenerator.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "PottsMeshGenerator.hpp"
#include "OnLatticeSimulation.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "VolumeConstraintPottsUpdateRule.hpp"
#include "AdhesionPottsUpdateRule.hpp"
#include "DifferentialAdhesionPottsUpdateRule.hpp"
#include "StemCellProliferativeType.hpp"
#include "CellMutationStatesCountWriter.hpp"

#include "CryptCellsGenerator.hpp"

#include "FakePetscSetup.hpp"

class TestPottsBasedSimulationsUsingDifferentiationTreeBasedCellCycleModel : public AbstractCellBasedTestSuite
{
public:
    void TestRunningPottsBasedSimulations()
    {
        /*** TREE ***/
        DifferentiationTree* tree = new DifferentiationTree(5);
        tree->addNewChild(0, 5);
        tree->addNewChild(0, 3);
        tree->addNewChild(2, 3);
        tree->addNewChild(2, 3);
        tree->addNewChild(3, 3);
        tree->addNewChild(3, 3);
        tree->addNewChild(4, 3);
        //tree->normaliseLength(.8);
        tree->colorise();
        /*** END TREE ***/

        EXIT_IF_PARALLEL;

        PottsMeshGenerator<2> generator(50, 2, 4, 50, 2, 4);  // Parameters are: lattice sites across; num elements across; element width; lattice sites up; num elements up; and element height
        PottsMesh<2>* p_mesh = generator.GetMesh();

        std::vector<CellPtr> cells;
        MAKE_PTR(StemCellProliferativeType, p_stem_type);
        CellsGenerator<DifferentiationTreeBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_stem_type);

        for (unsigned i=0; i<cells.size(); i++)
        {
            DifferentiationTreeBasedCellCycleModel* p_cell_cycle_model =
                                            new DifferentiationTreeBasedCellCycleModel(tree,0);
            CellPtr p_cell = cells.at(i);
            p_cell->SetCellCycleModel(p_cell_cycle_model);
            p_cell->InitialiseCellCycleModel();
        }
        PottsBasedCellPopulation<2> cell_population(*p_mesh, cells);

        cell_population.SetTemperature(0.1);
        cell_population.SetNumSweepsPerTimestep(1);

        OnLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("PottsDifferentiationTreeCellCycleModel");
        simulator.SetEndTime(10.0);

        MAKE_PTR(VolumeConstraintPottsUpdateRule<2>, p_volume_constraint_update_rule);
        p_volume_constraint_update_rule->SetMatureCellTargetVolume(16);
        p_volume_constraint_update_rule->SetDeformationEnergyParameter(0.2);
        simulator.AddPottsUpdateRule(p_volume_constraint_update_rule);
        MAKE_PTR(AdhesionPottsUpdateRule<2>, p_adhesion_update_rule);
        simulator.AddPottsUpdateRule(p_adhesion_update_rule);

        simulator.Solve();
        delete tree;
    }
/* next is not working. fix it! */
    void TestPottsCryptSimulationUsingDifferentiationTreeCellCycleModel()
    {
        /*** TREE ***/
        DifferentiationTree* tree = new DifferentiationTree(4);
        tree->addNewChild(0, 3);
        tree->addNewChild(0, 3);
        tree->addNewChild(2, 2);
        tree->addNewChild(2, 2);
        tree->addNewChild(3, 1);
        tree->addNewChild(3, 1);
        tree->addNewChild(4, 1);

        //tree->normaliseLength(.8);
        tree->colorise();
        /*** END TREE ***/

        PottsMeshGenerator<2> generator(100, 20, 5, 150, 29, 5, 1u, 1u, 1u, false, true, false, false);  // Parameters are: lattice sites across; num elements across; element width; lattice sites up; num elements up; and element height
        PottsMesh<2>* p_mesh = generator.GetMesh();

        std::vector<CellPtr> cells;
        MAKE_PTR(StemCellProliferativeType, p_stem_type);
        CellsGenerator<DifferentiationTreeBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_stem_type);
        boost::shared_ptr<AbstractCellProperty> p_state(CellPropertyRegistry::Instance()->Get<ApcTwoHitCellMutationState>());
        MAKE_PTR(CellLabel, p_label);

        for (unsigned i=0; i<cells.size(); i++)
        {
            DifferentiationTreeBasedCellCycleModel* p_cell_cycle_model;
            if (i < 60)
            {
                p_cell_cycle_model = new DifferentiationTreeBasedCellCycleModel(tree,1);
                cells[i]->AddCellProperty(p_label);
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

        PottsBasedCellPopulation<2> cell_population(*p_mesh, cells);

        //double crypt_height = 19;
        cell_population.SetNumSweepsPerTimestep(4);

        OnLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("PottsCryptTreeCellCycleModel");
        simulator.SetEndTime(10.0);

        MAKE_PTR(VolumeConstraintPottsUpdateRule<2>, p_volume_constraint_update_rule);
        p_volume_constraint_update_rule->SetMatureCellTargetVolume(16);
        p_volume_constraint_update_rule->SetDeformationEnergyParameter(0.2);
        simulator.AddPottsUpdateRule(p_volume_constraint_update_rule);
        MAKE_PTR(AdhesionPottsUpdateRule<2>, p_adhesion_update_rule);
        simulator.AddPottsUpdateRule(p_adhesion_update_rule);

        //MAKE_PTR_ARGS(SloughingCellKiller<2>, p_killer, (&cell_population, crypt_height));
        //simulator.AddCellKiller(p_killer);

        simulator.Solve();
        delete tree;
    }
};

#endif /* TESTPOTTSBASEDSIMULATIONSUSINGDIFFERENTIATIONTREECELLCYCLEMODEL_HPP_ */
