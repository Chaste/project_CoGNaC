#include "CellDifferentiationTypeWriter.hpp"
#include "AbstractCellPopulation.hpp"

#include "Exception.hpp"
#include <iostream>

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
CellDifferentiationTypeWriter<ELEMENT_DIM, SPACE_DIM>::CellDifferentiationTypeWriter()
    : AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>("results.vizdiffcolour")
{
    this->mVtkCellDataName = "Differentiation Colour";
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double CellDifferentiationTypeWriter<ELEMENT_DIM, SPACE_DIM>::GetCellDataForVtkOutput(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    double colour = (double) pCell->GetCellProliferativeType()->GetColour();
    // Set colour dependent on cell differentiation node
    try{
        double diff_colour = pCell->GetCellData()->GetItem("Colour");
        if (diff_colour >= 0.0)
        {
            colour = diff_colour;
        }
    } catch (Exception& error)
    {

    }
    // Set colour dependent on cell mutation state
    if (!pCell->GetMutationState()->IsType<WildTypeCellMutationState>())
    {
        colour = (double) pCell->GetMutationState()->GetColour();
    }
    if (pCell->HasCellProperty<CellLabel>())
    {
        CellPropertyCollection collection = pCell->rGetCellPropertyCollection().GetProperties<CellLabel>();
        boost::shared_ptr<CellLabel> p_label = boost::static_pointer_cast<CellLabel>(collection.GetProperty());
        colour = (double) p_label->GetColour() /* + 1.0 */;
    }
    if (pCell->HasCellProperty<ApoptoticCellProperty>() || pCell->HasApoptosisBegun())
    {
        ///\todo: replace this hard-coded 6 with the ApoptoticCellProperty member mColour?
        colour = 6.0;
    }
    return colour;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellDifferentiationTypeWriter<ELEMENT_DIM, SPACE_DIM>::VisitCell(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    double colour = (double) pCell->GetCellProliferativeType()->GetColour();
    // Set colour dependent on cell differentiation node
    try{
        double diff_colour = pCell->GetCellData()->GetItem("Colour");
        if (diff_colour >= 0.0)
        {
            colour = diff_colour;
        }
    } catch (Exception& error)
    {

    }
    // Set colour dependent on cell mutation state
    if (!pCell->GetMutationState()->IsType<WildTypeCellMutationState>())
    {
        colour = (double) pCell->GetMutationState()->GetColour();
    }
    if (pCell->HasCellProperty<CellLabel>())
    {
        CellPropertyCollection collection = pCell->rGetCellPropertyCollection().GetProperties<CellLabel>();
        boost::shared_ptr<CellLabel> p_label = boost::static_pointer_cast<CellLabel>(collection.GetProperty());
        colour = (double) p_label->GetColour()/*+1.0*/;
    }
    if (pCell->HasCellProperty<ApoptoticCellProperty>() || pCell->HasApoptosisBegun())
    {
        ///\todo: replace this hard-coded 6 with the ApoptoticCellProperty member mColour?
        colour = 6.0;
    }

    *this->mpOutStream << colour << " ";
}

// Explicit instantiation
template class CellDifferentiationTypeWriter<1,1>;
template class CellDifferentiationTypeWriter<1,2>;
template class CellDifferentiationTypeWriter<2,2>;
template class CellDifferentiationTypeWriter<1,3>;
template class CellDifferentiationTypeWriter<2,3>;
template class CellDifferentiationTypeWriter<3,3>;

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
EXPORT_TEMPLATE_CLASS_ALL_DIMS(CellDifferentiationTypeWriter)
