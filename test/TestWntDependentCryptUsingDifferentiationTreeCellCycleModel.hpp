#ifndef TESTWNTDEPENDENTCRYPTUSINGDIFFERENTIATIONTREECELLCYCLEMODEL_HPP_
#define TESTWNTDEPENDENTCRYPTUSINGDIFFERENTIATIONTREECELLCYCLEMODEL_HPP_

#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "SmartPointers.hpp"
#include "AbstractCellBasedTestSuite.hpp"

#include "DifferentiationTree.hpp"
#include "DifferentiationTreeNode.hpp"
#include "DifferentiationTreeBasedWntCellCycleModel.hpp"
#include "CellDifferentiationTypeWriter.hpp"

#include "DifferentiationTreeCryptCellsGenerator.hpp"
#include "CryptCellsGenerator.hpp"
#include "CylindricalHoneycombMeshGenerator.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "CryptSimulation2d.hpp"
#include "WntConcentration.hpp"
#include "SloughingCellKiller.hpp"
#include "CellMutationStatesCountWriter.hpp"
#include "VoronoiDataWriter.hpp"
#include "FakePetscSetup.hpp"


class TestWntDependentCryptUsingDifferentiationTreeCellCycleModel : public AbstractCellBasedTestSuite
{
public:
    void TestDifferentiationTreeCryptWithMutations() throw(Exception)
    {
        /*** TREE ***/
        DifferentiationTree* tree = new DifferentiationTree(9);
        tree->addNewChild(0,10);
        tree->addNewChild(0,10);
        tree->addNewChild(2,11);
        tree->addNewChild(2,11);
        tree->addNewChild(3,10);
        tree->addNewChild(3,11);
        tree->addNewChild(4,10);

        tree->normaliseLength(12);
        tree->colorise();
        /*** END TREE ***/

        CylindricalHoneycombMeshGenerator generator(20, 20, 3);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();

        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        std::vector<CellPtr> cells;
        CryptCellsGenerator<DifferentiationTreeBasedWntCellCycleModel> cells_generator;
        cells_generator.Generate(cells, p_mesh, location_indices, true);

        for (unsigned i=0; i<cells.size(); i++)
        {
            DifferentiationTreeBasedWntCellCycleModel* p_cell_cycle_model;
            if (i<24)
            {
                p_cell_cycle_model = new DifferentiationTreeBasedWntCellCycleModel(tree,1);

            }
            else if(i<48)
            {
                p_cell_cycle_model = new DifferentiationTreeBasedWntCellCycleModel(tree,0);
            }
            else if(i<112)
            {
                p_cell_cycle_model = new DifferentiationTreeBasedWntCellCycleModel(tree,2);
            }
            else if(i<144)
            {
                p_cell_cycle_model = new DifferentiationTreeBasedWntCellCycleModel(tree,3);
            }
            else
            {
                p_cell_cycle_model = new DifferentiationTreeBasedWntCellCycleModel(tree,7);
            }
            CellPtr p_cell = cells.at(i);
            double birth_time = p_cell->GetBirthTime();
            p_cell_cycle_model->SetDimension(2);
            p_cell->SetCellCycleModel(p_cell_cycle_model);
            p_cell->SetBirthTime(birth_time);
            p_cell->InitialiseCellCycleModel();
        }

        MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, location_indices);
        cell_population.AddCellWriter<CellDifferentiationTypeWriter>();
        cell_population.SetWriteVtkAsPoints(false);
        cell_population.AddCellPopulationCountWriter<CellMutationStatesCountWriter>();
        cell_population.AddPopulationWriter<VoronoiDataWriter>();

        double crypt_height = 18.0;

        WntConcentration<2>::Instance()->SetType(LINEAR);
        WntConcentration<2>::Instance()->SetCellPopulation(cell_population);
        WntConcentration<2>::Instance()->SetCryptLength(crypt_height);

        CryptSimulation2d simulator(cell_population);
        simulator.SetOutputDirectory("DifferentiationTreeWntCryptWithMutations");
        simulator.SetSamplingTimestepMultiple(12);
        simulator.SetEndTime(20);

        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        simulator.AddForce(p_linear_force);
        MAKE_PTR_ARGS(SloughingCellKiller<2>, p_killer, (&cell_population, crypt_height));
        simulator.AddCellKiller(p_killer);

        simulator.Solve();

        WntConcentration<2>::Destroy();

    }
    void notactiveTestMeshBasedCryptWithMutations() throw(Exception)
    {
        CylindricalHoneycombMeshGenerator generator(6, 9, 2);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();

        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        std::vector<CellPtr> cells;
        CryptCellsGenerator<SimpleWntCellCycleModel> cells_generator;
        cells_generator.Generate(cells, p_mesh, location_indices, true);

        boost::shared_ptr<AbstractCellProperty> p_state(CellPropertyRegistry::Instance()->Get<ApcTwoHitCellMutationState>());

        MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, location_indices);
        cell_population.AddCellPopulationCountWriter<CellMutationStatesCountWriter>();

        double crypt_height = 8.0;

        WntConcentration<2>::Instance()->SetType(LINEAR);
        WntConcentration<2>::Instance()->SetCellPopulation(cell_population);
        WntConcentration<2>::Instance()->SetCryptLength(crypt_height);

        CryptSimulation2d simulator(cell_population);
        simulator.SetOutputDirectory("MeshBasedCryptWithMutations");
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

       WntConcentration<2>::Destroy();
    }

    void testCryptGenerator()
    {
        /*** TREE ***/
        DifferentiationTree* tree = new DifferentiationTree(9);
        tree->addNewChild(0,10);
        tree->addNewChild(0,10);
        tree->addNewChild(2,11);
        tree->addNewChild(2,11);
        tree->addNewChild(3,10);
        tree->addNewChild(3,11);
        tree->addNewChild(4,10);

        tree->normaliseLength(12);
        tree->colorise();
        /*** END TREE ***/

        CylindricalHoneycombMeshGenerator generator(6, 11, 2);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();

        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        double crypt_height = 10.0;
        std::vector<CellPtr> cells;
        DifferentiationTreeCryptCellsGenerator<DifferentiationTreeBasedWntCellCycleModel> cells_generator;
        cells_generator.Generate(cells, p_mesh, tree, location_indices, true, crypt_height);

        boost::shared_ptr<AbstractCellProperty> p_state(CellPropertyRegistry::Instance()->Get<ApcTwoHitCellMutationState>());

        MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, location_indices);
        cell_population.AddCellWriter<CellDifferentiationTypeWriter>();
        cell_population.SetWriteVtkAsPoints(false);
        cell_population.AddCellPopulationCountWriter<CellMutationStatesCountWriter>();
        cell_population.AddPopulationWriter<VoronoiDataWriter>();

        WntConcentration<2>::Instance()->SetType(LINEAR);
        WntConcentration<2>::Instance()->SetCellPopulation(cell_population);
        WntConcentration<2>::Instance()->SetCryptLength(crypt_height);

        CryptSimulation2d simulator(cell_population);
        simulator.SetOutputDirectory("GeneratorCCM");
        simulator.SetSamplingTimestepMultiple(12);
        simulator.SetEndTime(10);

        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        simulator.AddForce(p_linear_force);
        MAKE_PTR_ARGS(SloughingCellKiller<2>, p_killer, (&cell_population, crypt_height));
        simulator.AddCellKiller(p_killer);

        simulator.Solve();
        WntConcentration<2>::Destroy();

    }
};

#endif /* TESTWNTDEPENDENTCRYPTUSINGDIFFERENTIATIONTREECELLCYCLEMODEL_HPP_ */
