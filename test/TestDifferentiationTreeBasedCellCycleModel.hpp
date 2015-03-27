#ifndef TESTDIFFERENTIATIONTREEBASEDCELLCYCLEMODEL_HPP_
#define TESTDIFFERENTIATIONTREEBASEDCELLCYCLEMODEL_HPP_

#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "SmartPointers.hpp"
#include "AbstractCellBasedTestSuite.hpp"

#include "DifferentiationTree.hpp"
#include "DifferentiationTreeNode.hpp"
#include "DifferentiationTreeBasedCellCycleModel.hpp"

#include "DifferentiatedCellProliferativeType.hpp"
#include "CheckReadyToDivideAndPhaseIsUpdated.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "WildTypeCellMutationState.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "StemCellProliferativeType.hpp"
//#include "TransitCellProliferativeType.hpp"
#include "OffLatticeSimulation.hpp"
#include "VoronoiDataWriter.hpp"
#include "CellDifferentiationTypeWriter.hpp"
#include "CellsGenerator.hpp"

#include "HoneycombVertexMeshGenerator.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "NagaiHondaForce.hpp"
#include "SimpleTargetAreaModifier.hpp"

#include "FakePetscSetup.hpp"
#include <vector>
#include <set>

class TestDifferentiationTreeBasedCellCycleModel : public AbstractCellBasedTestSuite
{
public:
    void TestCorrectnessDifferentiationTreeBasedCellCycleModel()
    {
        TS_ASSERT_THROWS_NOTHING(DifferentiationTreeBasedCellCycleModel cell_model3);

/***BUILD THE TREE***/
        DifferentiationTree* tree = new DifferentiationTree(7);
        tree->addNewChild(0,5);
        tree->addNewChild(0,4);
        tree->addNewChild(2,5);
        tree->addNewChild(2,5);
        tree->addNewChild(3,5);
        tree->addNewChild(3,5);
        tree->addNewChild(4,5);
/***END***/

        unsigned num_cells = 10;
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        for (unsigned i=0; i<num_cells; i++)
        {
            DifferentiationTreeBasedCellCycleModel* p_cell_cycle_model =
                    new DifferentiationTreeBasedCellCycleModel(tree,0);
            CellPtr p_cell(new Cell(p_state, p_cell_cycle_model));
            p_cell->InitialiseCellCycleModel();
            cells.push_back(p_cell);
        }

        DifferentiationTreeBasedCellCycleModel* p_my_model =
                new DifferentiationTreeBasedCellCycleModel(tree,0);
        CellPtr p_my_cell(new Cell(p_state, p_my_model));
        p_my_cell->InitialiseCellCycleModel();

        unsigned num_steps = 100;
        double cell_cycle_time = cells[0]->GetCellCycleModel()->GetG1Duration()
                                        + cells[0]->GetCellCycleModel()->GetSG2MDuration();

        TS_ASSERT_DELTA(cell_cycle_time, tree->getNode(0)->getCellCycleLength(), 0.0001);
        TS_ASSERT_DELTA(cell_cycle_time, 7.0, 0.0001);

        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(cell_cycle_time, num_steps);

        for (unsigned i=0; i<num_steps; i++)
        {
            SimulationTime::Instance()->IncrementTimeOneStep();
            CheckReadyToDivideAndPhaseIsUpdated(p_my_model, tree->getNode(0)->getCellCycleLength()/4);
        }

        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "diff_tree_based_cell_cycle_model.arch";
        {
            SimulationTime::Destroy();
            SimulationTime::Instance()->SetStartTime(0.0);
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(4.0, 10);

            DifferentiationTreeBasedCellCycleModel* p_model =
                    new DifferentiationTreeBasedCellCycleModel(tree,2);
            CellPtr p_cell(new Cell(p_state, p_model));
            p_cell->InitialiseCellCycleModel();

            p_simulation_time->IncrementTimeOneStep();
            p_simulation_time->IncrementTimeOneStep();
            p_simulation_time->IncrementTimeOneStep();

            p_model->SetBirthTime(-1.0);
            p_model->ReadyToDivide();
            TS_ASSERT_EQUALS(p_model->GetCurrentCellCyclePhase(), S_PHASE);

            CellPtr const p_const_cell = p_cell;

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);
            output_arch << p_const_cell;
            }

            {
            SimulationTime::Destroy();
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetStartTime(0.0);
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

            CellPtr p_cell;

            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            input_arch >> p_cell;

            AbstractCellCycleModel* p_model = p_cell->GetCellCycleModel();

            TS_ASSERT_DELTA(p_model->GetBirthTime(), -1.0, 1e-12);
            TS_ASSERT_DELTA(p_model->GetAge(), 2.2, 1e-12);
            TS_ASSERT_EQUALS(p_model->GetCurrentCellCyclePhase(), S_PHASE);
        }
        delete tree;
    }

    void Test2DParaviewMeshBasedSimulationUsingDifferentiationTreeBasedCellCycleModel()
    {
        DifferentiationTree* tree = new DifferentiationTree(.5);
        tree->addNewChild(0,.2);
        tree->addNewChild(0,.2);
        tree->addNewChild(2,.2);
        tree->addNewChild(2,.2);
        tree->addNewChild(3,.2);
        tree->addNewChild(3,.2);
        tree->addNewChild(4,.2);

        //tree->normaliseLength(.8);
        tree->colorise();


        HoneycombMeshGenerator generator(10, 10, 2);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        std::vector<CellPtr> cells;
        MAKE_PTR(StemCellProliferativeType, p_stem_type);
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        CellsGenerator<DifferentiationTreeBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, location_indices.size(), p_stem_type);

        for (unsigned i=0; i<cells.size(); i++)
        {
            DifferentiationTreeBasedCellCycleModel* p_cell_cycle_model =
                                new DifferentiationTreeBasedCellCycleModel(tree,0);
            CellPtr p_cell = cells.at(i);
            if (i==50)
            {
                delete p_cell_cycle_model;
                p_cell_cycle_model = new DifferentiationTreeBasedCellCycleModel(tree,0);
                p_cell->SetCellProliferativeType(p_transit_type);
            }
            //DifferentiationTreeBasedCellCycleModel* p_cell_cycle_model =
            //                    new DifferentiationTreeBasedCellCycleModel(tree,&node0);
            //CellPtr p_cell = cells.at(i);
            p_cell->SetCellCycleModel(p_cell_cycle_model);
            p_cell->InitialiseCellCycleModel();
        }

        MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, location_indices);
        cell_population.AddCellWriter<CellDifferentiationTypeWriter>();
        cell_population.SetWriteVtkAsPoints(false);
        cell_population.AddPopulationWriter<VoronoiDataWriter>();

        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("MeshBasedSimulationWithDifferentiationTreeBasedCellCycleModel");
        simulator.SetEndTime(1.0);

        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        p_linear_force->SetCutOffLength(1.5);
        simulator.AddForce(p_linear_force);

        simulator.Solve();
        delete tree;
    }

    void Test2DParaviewNodeBasedSimulationUsingDifferentiationTreeBasedCellCycleModel() throw(Exception)
    {
        HoneycombMeshGenerator generator(10, 10, 0);
        TetrahedralMesh<2,2>* p_generating_mesh = generator.GetMesh();

        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(*p_generating_mesh, 1.5);

        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);

        DifferentiationTree* tree = new DifferentiationTree(.1);
        tree->addNewChild(0,.05);
        tree->addNewChild(0,.05);
        tree->addNewChild(2,.05);
        tree->addNewChild(2,.05);
        tree->addNewChild(3,.05);
        tree->addNewChild(3,.05);
        tree->addNewChild(4,.05);

        tree->colorise();

        CellsGenerator<DifferentiationTreeBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes());

        for (unsigned i=0; i<cells.size(); i++)
        {
           DifferentiationTreeBasedCellCycleModel* p_cell_cycle_model =
                   new DifferentiationTreeBasedCellCycleModel(tree,0);
           CellPtr p_cell = cells.at(i);
           p_cell->SetCellCycleModel(p_cell_cycle_model);
           p_cell->InitialiseCellCycleModel();
        }
        NodeBasedCellPopulation<2> cell_population(mesh, cells);

        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("NodeBasedSimulationWithDifferentiationTreeBasedCellCycleModel");
        simulator.SetEndTime(0.25);

        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        p_linear_force->SetCutOffLength(1.5);
        simulator.AddForce(p_linear_force);

        simulator.Solve();
        delete tree;

   }

   void Test2DParaviewVertexBasedSimulationsUsingDifferentiationTreeBasedCellCycleModel()
   {
       DifferentiationTree* tree = new DifferentiationTree(.05);
       tree->addNewChild(0,.05);
       tree->addNewChild(0,.05);
       tree->addNewChild(2,.05);
       tree->addNewChild(2,.05);
       tree->addNewChild(3,.05);
       tree->addNewChild(3,.05);
       tree->addNewChild(4,.05);

       tree->colorise();

       HoneycombVertexMeshGenerator generator(2, 2);
       MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

       std::vector<CellPtr> cells;
       //CellsGenerator<DifferentiationTreeBasedCellCycleModel, 2> cells_generator;
       //cells_generator.GenerateBasic(cells, p_mesh->GetNumElements());
       MAKE_PTR(WildTypeCellMutationState, p_state);
       cells.reserve(p_mesh->GetNumElements());

       for (unsigned i=0; i<p_mesh->GetNumElements(); i++)
       {
          DifferentiationTreeBasedCellCycleModel* p_cell_cycle_model =
                  new DifferentiationTreeBasedCellCycleModel(tree,0);
          p_cell_cycle_model->SetDimension(2);
          CellPtr p_cell(new Cell(p_state, p_cell_cycle_model));
          p_cell->InitialiseCellCycleModel();
          cells.push_back(p_cell);
       }

       VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

       OffLatticeSimulation<2> simulator(cell_population);
       simulator.SetOutputDirectory("VertexSimulationUsingDifferentiationTreeBasedCellCycleModel");
       simulator.SetEndTime(0.1);

       MAKE_PTR(NagaiHondaForce<2>, p_nagai_honda_force);
       simulator.AddForce(p_nagai_honda_force);

       MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
       simulator.AddSimulationModifier(p_growth_modifier);

       simulator.Solve();
       delete tree;
   }

};
#endif /* TESTDIFFERENTIATIONTREEBASEDCELLCYCLEMODEL_HPP_ */
