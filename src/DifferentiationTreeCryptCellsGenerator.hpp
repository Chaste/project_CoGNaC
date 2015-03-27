#ifndef DIFFERENTIATIONTREECRYPTCELLSGENERATOR_HPP_
#define DIFFERENTIATIONTREECRYPTCELLSGENERATOR_HPP_

#include <boost/mpl/integral_c.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/mpl/if.hpp>

#include "CellsGenerator.hpp"
#include "CryptCellsGenerator.hpp"

#include "CellPropertyRegistry.hpp"
#include "TetrahedralMesh.hpp"
#include "VertexMesh.hpp"
#include "PottsMesh.hpp"

#include "DifferentiationTree.hpp"
#include "DifferentiationTreeBasedCellCycleModel.hpp"
#include "DifferentiationTreeBasedWntCellCycleModel.hpp"
#include "StochasticDurationGenerationBasedCellCycleModel.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "TysonNovakCellCycleModel.hpp"
#include "WntCellCycleModel.hpp"
#include "SimpleWntCellCycleModel.hpp"
#include "StochasticWntCellCycleModel.hpp"
#include "VanLeeuwen2009WntSwatCellCycleModelHypothesisOne.hpp"
#include "VanLeeuwen2009WntSwatCellCycleModelHypothesisTwo.hpp"
#include "Exception.hpp"
#include "StemCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"


/**
 * A subclass of CellsGenerator that generates cells for crypt simulations.
 *
 * It is templated over types of cell-cycle model.
 */
template<class CELL_CYCLE_MODEL>
class DifferentiationTreeCryptCellsGenerator : public CellsGenerator<CELL_CYCLE_MODEL,2>
{
public:

    /**
     * Generates cells of a specified cell cycle type under the correct
     * crypt conditions and gives random ages if required,
     * or gives them an age of 0.0 - creates least work for solver startup.
     *
     * @param rCells  An empty cells vector for this function to fill up
     * @param pMesh  The crypt mesh (can be cylindrical)
     * @param locationIndices the node indices corresponding to real cells
     * @param randomBirthTimes  Whether to assign the cells random birth times
     *    (this can be expensive computationally with ODE models)
     * @param y0
     * @param y1
     * @param y2
     * @param y3
     * @param initialiseCells  whether to initialise the cell-cycle models as each
     *   cell is created
     */
    void Generate(
                  std::vector<CellPtr>& rCells,
                  AbstractMesh<2,2>* pMesh,
                  const DifferentiationTree* crypt_tree,
                  const std::vector<unsigned> locationIndices,
                  bool randomBirthTimes,
                  double crypt_height = 10,
                  double y0 = 0.1,
                  double y1 = 0.2,
                  double y2 = 0.6,
                  double y3 = 0.8,
                  bool initialiseCells = true
                  );
};

template<class CELL_CYCLE_MODEL>
void DifferentiationTreeCryptCellsGenerator<CELL_CYCLE_MODEL>::Generate(
                                      std::vector<CellPtr>& rCells,
                                      AbstractMesh<2,2>* pMesh,
                                      const DifferentiationTree* crypt_tree,
                                      const std::vector<unsigned> locationIndices,
                                      bool randomBirthTimes,
                                      double crypt_height,
                                      double y0,
                                      double y1,
                                      double y2,
                                      double y3,
                                      bool initialiseCells)
{
    y0 *= crypt_height;
    y1 *= crypt_height;
    y2 *= crypt_height;
    y3 *= crypt_height;

    rCells.clear();

    RandomNumberGenerator* p_random_num_gen = RandomNumberGenerator::Instance();

    /* Verify if crypt tree is a correct tree*/
    DifferentiationTree* tree = new DifferentiationTree(4);
    tree->addNewChild(0);
    tree->addNewChild(0);
    tree->addNewChild(2);
    tree->addNewChild(2);
    tree->addNewChild(3);
    tree->addNewChild(3);
    tree->addNewChild(4);
    /* */
    assert(tree->topologyTreeCompare(crypt_tree));
    delete tree;

    unsigned child = crypt_tree->getRoot()->getChildren().at(0);
    unsigned paneth = 0, ta1 = 0, ta2_a = 0, ta2_b = 0, enterocyte = 0, enteroendocrine = 0, goblet = 0;
    if (crypt_tree->getNode(child)->getNumberOfChildren() == 0)
    {
        paneth = child;
        ta1 = crypt_tree->getRoot()->getChildren().at(1);
    }else
    {
        paneth = crypt_tree->getRoot()->getChildren().at(1);
        ta1 = child;
    }
    child = crypt_tree->getNode(ta1)->getChildren().at(0);
    if (crypt_tree->getNode(child)->getNumberOfChildren()==2)
    {
        ta2_a = child;
        ta2_b = crypt_tree->getNode(ta1)->getChildren().at(1);
    }
    else
    {
        ta2_a = crypt_tree->getNode(ta1)->getChildren().at(1);
        ta2_b = child;
    }
    enterocyte = crypt_tree->getNode(ta2_a)->getChildren().at(0);
    enteroendocrine = crypt_tree->getNode(ta2_a)->getChildren().at(1);
    goblet = crypt_tree->getNode(ta2_b)->getChildren().at(0);

    unsigned mesh_size;
    if (dynamic_cast<TetrahedralMesh<2,2>*>(pMesh))
    {
        mesh_size = pMesh->GetNumNodes();
        unsigned num_cells = locationIndices.empty() ? pMesh->GetNumNodes() : locationIndices.size();
        rCells.reserve(num_cells);
    }
    else if (dynamic_cast<PottsMesh<2>*>(pMesh))
    {
        mesh_size = static_cast<PottsMesh<2>*>(pMesh)->GetNumElements();
        rCells.reserve(mesh_size);
    }
    else
    {
        // Note the double brackets, to stop the macro thinking is has two arguments.
        assert((dynamic_cast<VertexMesh<2,2>*>(pMesh)));
        mesh_size = static_cast<VertexMesh<2,2>*>(pMesh)->GetNumElements();
        rCells.reserve(mesh_size);
    }

    boost::shared_ptr<AbstractCellProperty> p_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());

    // Loop over the mesh and populate rCells
    for (unsigned i=0; i<mesh_size; i++)
    {
        // Find the location of this cell
        double y = 0.0;
        if (dynamic_cast<TetrahedralMesh<2,2>*>(pMesh))
        {
            if (locationIndices.empty())
            {
                y = pMesh->GetNode(i)->GetPoint().rGetLocation()[1];

            }
            else if (std::find(locationIndices.begin(), locationIndices.end(), i) != locationIndices.end())
            {
                y = pMesh->GetNode(i)->GetPoint().rGetLocation()[1];
            }
        }
        else if (dynamic_cast<PottsMesh<2>*>(pMesh))
        {
            y = static_cast<PottsMesh<2>*>(pMesh)->GetCentroidOfElement(i)[1];
        }
        else
        {
            // Note the double brackets, to stop the macro thinking is has two arguments.
            assert((dynamic_cast<VertexMesh<2,2>*>(pMesh)));
            y = static_cast<VertexMesh<2,2>*>(pMesh)->GetCentroidOfElement(i)[1];
        }

        // Create a cell-cycle model and set the spatial dimension
        CELL_CYCLE_MODEL* p_cell_cycle_model = new CELL_CYCLE_MODEL;
        /*
        if (dynamic_cast<DifferentiationTreeBasedCellCycleModel*>(p_cell_cycle_model))
        {
            p_cell_cycle_model = new DifferentiationTreeBasedCellCycleModel();
        } else if (dynamic_cast<DifferentiationTreeBasedWntCellCycleModel*>(p_cell_cycle_model))
        {
            p_cell_cycle_model = new DifferentiationTreeBasedWntCellCycleModel();
        } else
        {
            p_cell_cycle_model =
        }
        */
        p_cell_cycle_model->SetDimension(2);

        double typical_transit_cycle_time = p_cell_cycle_model->GetAverageTransitCellCycleTime();
        double typical_stem_cycle_time = p_cell_cycle_model->GetAverageStemCellCycleTime();

        double birth_time = 0.0;
        if (randomBirthTimes)
        {
            birth_time = -p_random_num_gen->ranf();
        }

        // Set the cell-cycle model's generation if required
        unsigned generation = 4;
        if (y <= y0)
        {
            if (dynamic_cast<DifferentiationTreeBasedCellCycleModel*>(p_cell_cycle_model))
            {
                dynamic_cast<DifferentiationTreeBasedCellCycleModel*>(p_cell_cycle_model)->SetTreeAndType(crypt_tree, paneth);
            } else if (dynamic_cast<DifferentiationTreeBasedWntCellCycleModel*>(p_cell_cycle_model))
            {
                dynamic_cast<DifferentiationTreeBasedWntCellCycleModel*>(p_cell_cycle_model)->SetTreeAndType(crypt_tree, paneth);
            }
            birth_time *= crypt_tree->getNode(paneth)->getCellCycleLength();
            generation = 0;
        }
        else if (y < y1)
        {
            if (dynamic_cast<DifferentiationTreeBasedCellCycleModel*>(p_cell_cycle_model))
            {
                dynamic_cast<DifferentiationTreeBasedCellCycleModel*>(p_cell_cycle_model)->SetTreeAndType(crypt_tree, 0);
            } else if (dynamic_cast<DifferentiationTreeBasedWntCellCycleModel*>(p_cell_cycle_model))
            {
                dynamic_cast<DifferentiationTreeBasedWntCellCycleModel*>(p_cell_cycle_model)->SetTreeAndType(crypt_tree, 0);
            }
            birth_time *= crypt_tree->getRoot()->getCellCycleLength();
            generation = 1;
        }
        else if (y < y2)
        {
            if (dynamic_cast<DifferentiationTreeBasedCellCycleModel*>(p_cell_cycle_model))
            {
                dynamic_cast<DifferentiationTreeBasedCellCycleModel*>(p_cell_cycle_model)->SetTreeAndType(crypt_tree, ta1);
            } else if (dynamic_cast<DifferentiationTreeBasedWntCellCycleModel*>(p_cell_cycle_model))
            {
                dynamic_cast<DifferentiationTreeBasedWntCellCycleModel*>(p_cell_cycle_model)->SetTreeAndType(crypt_tree, ta1);
            }
            birth_time *= crypt_tree->getNode(ta1)->getCellCycleLength();
            generation = 2;
        }
        else if (y < y3)
        {
            if (dynamic_cast<DifferentiationTreeBasedCellCycleModel*>(p_cell_cycle_model))
            {
                if (RandomNumberGenerator::Instance()->ranf() <= .5)
                {
                    dynamic_cast<DifferentiationTreeBasedCellCycleModel*>(p_cell_cycle_model)->SetTreeAndType(crypt_tree, ta2_a);
                    birth_time *= crypt_tree->getNode(ta2_a)->getCellCycleLength(); // hours
                }
                else
                {
                    dynamic_cast<DifferentiationTreeBasedCellCycleModel*>(p_cell_cycle_model)->SetTreeAndType(crypt_tree, ta2_b);
                    birth_time *= crypt_tree->getNode(ta2_b)->getCellCycleLength(); // hours
                }
            }
            else if (dynamic_cast<DifferentiationTreeBasedWntCellCycleModel*>(p_cell_cycle_model))
            {
                if (RandomNumberGenerator::Instance()->ranf() <= .5)
                {
                    dynamic_cast<DifferentiationTreeBasedWntCellCycleModel*>(p_cell_cycle_model)->SetTreeAndType(crypt_tree, ta2_a);
                    birth_time *= crypt_tree->getNode(ta2_a)->getCellCycleLength(); // hours
                }
                else
                {
                    dynamic_cast<DifferentiationTreeBasedWntCellCycleModel*>(p_cell_cycle_model)->SetTreeAndType(crypt_tree, ta2_b);
                    birth_time *= crypt_tree->getNode(ta2_a)->getCellCycleLength(); // hours
                }
            } else  generation = 3;
        }
        else
        {
            if (dynamic_cast<DifferentiationTreeBasedCellCycleModel*>(p_cell_cycle_model))
            {
                double random_number = RandomNumberGenerator::Instance()->ranf();
                if (random_number <= 1/3)
                {
                    dynamic_cast<DifferentiationTreeBasedCellCycleModel*>(p_cell_cycle_model)->SetTreeAndType(crypt_tree, enterocyte);
                    birth_time *= crypt_tree->getNode(enterocyte)->getCellCycleLength(); // hours
                }
                else if (random_number <= 2/3)
                {
                    dynamic_cast<DifferentiationTreeBasedCellCycleModel*>(p_cell_cycle_model)->SetTreeAndType(crypt_tree, enteroendocrine);
                    birth_time *= crypt_tree->getNode(enteroendocrine)->getCellCycleLength(); // hours
                }
                else
                {
                    dynamic_cast<DifferentiationTreeBasedCellCycleModel*>(p_cell_cycle_model)->SetTreeAndType(crypt_tree, goblet);
                    birth_time *= crypt_tree->getNode(goblet)->getCellCycleLength(); // hours
                }
            } else if (dynamic_cast<DifferentiationTreeBasedWntCellCycleModel*>(p_cell_cycle_model))
            {
                double random_number = RandomNumberGenerator::Instance()->ranf();
                if (random_number <= 1/3)
                {
                    dynamic_cast<DifferentiationTreeBasedWntCellCycleModel*>(p_cell_cycle_model)->SetTreeAndType(crypt_tree, enterocyte);
                    birth_time *= crypt_tree->getNode(enterocyte)->getCellCycleLength(); // hours
                }
                else if (random_number <= 2/3)
                {
                    dynamic_cast<DifferentiationTreeBasedWntCellCycleModel*>(p_cell_cycle_model)->SetTreeAndType(crypt_tree, enteroendocrine);
                    birth_time *= crypt_tree->getNode(enteroendocrine)->getCellCycleLength(); // hours
                }
                else
                {
                    dynamic_cast<DifferentiationTreeBasedWntCellCycleModel*>(p_cell_cycle_model)->SetTreeAndType(crypt_tree, goblet);
                    birth_time *= crypt_tree->getNode(goblet)->getCellCycleLength(); // hours
                }
            }
        }
        if (dynamic_cast<AbstractSimpleGenerationBasedCellCycleModel*>(p_cell_cycle_model))
        {
            dynamic_cast<AbstractSimpleGenerationBasedCellCycleModel*>(p_cell_cycle_model)->SetGeneration(generation);
        }

        // Create a cell
        CellPtr p_cell(new Cell(p_state, p_cell_cycle_model));

        /* Set the cell's proliferative type, dependent on its height up the crypt and whether it can terminally differentiate
        if (y <= y0)
        {
            p_cell->SetCellProliferativeType(CellPropertyRegistry::Instance()->Get<StemCellProliferativeType>());
        }
        else
        {
            p_cell->SetCellProliferativeType(CellPropertyRegistry::Instance()->Get<TransitCellProliferativeType>());
            if (y >= y3 && p_cell_cycle_model->CanCellTerminallyDifferentiate())
            {
                p_cell->SetCellProliferativeType(CellPropertyRegistry::Instance()->Get<DifferentiatedCellProliferativeType>());
            }
        }
        */

        // Initialise the cell-cycle model if this is required
        if (initialiseCells)
        {
            p_cell->InitialiseCellCycleModel();
        }

        // Set the cell's birth time
        if (!dynamic_cast<DifferentiationTreeBasedCellCycleModel*>(p_cell_cycle_model)
            && !dynamic_cast<DifferentiationTreeBasedWntCellCycleModel*>(p_cell_cycle_model))
        {
            if (y <= y1 && y > y0)
            {
                birth_time *= typical_stem_cycle_time; // hours
            }
            else
            {
                birth_time *= typical_transit_cycle_time; // hours
            }
        }
        p_cell->SetBirthTime(birth_time);

        if (locationIndices.empty())
        {
            rCells.push_back(p_cell);
        }
        else if (std::find(locationIndices.begin(), locationIndices.end(), i) != locationIndices.end())
        {
            rCells.push_back(p_cell);
        }
    }
}

#endif /* DIFFERENTIATIONTREECRYPTCELLSGENERATOR_HPP_ */


