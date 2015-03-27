#include "DifferentiationTreeBasedCellCycleModel.hpp"
#include "StemCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include <iostream>

DifferentiationTreeBasedCellCycleModel::DifferentiationTreeBasedCellCycleModel()
     : AbstractSimpleCellCycleModel(),
       mpDifferentiationTree(NULL),
       mDifferentiationType(0)
{
    mG1Duration = 2.0; //? Taken from Meineke model
    mSDuration = 8.0; //? Taken from Meineke model
    mG2Duration = 1.0; //? Taken from Meineke model
    mMDuration = 1.0; //? Taken from Meineke model
    mStemCellG1Duration = 2.0;
    mTransitCellG1Duration = 2.0;
}

DifferentiationTreeBasedCellCycleModel::DifferentiationTreeBasedCellCycleModel
(
        const DifferentiationTree* differentiation_tree,
        unsigned differentiation_type
)
:       AbstractSimpleCellCycleModel(),
        mpDifferentiationTree(differentiation_tree),
        mDifferentiationType(differentiation_type)
{
    assert(mpDifferentiationTree);
    assert(differentiation_type < mpDifferentiationTree->size());
}

DifferentiationTreeBasedCellCycleModel::~DifferentiationTreeBasedCellCycleModel()
{
    mpDifferentiationTree = NULL;
}

void DifferentiationTreeBasedCellCycleModel::SetTreeAndType
(
        const DifferentiationTree* differentiation_tree,
        unsigned differentiation_type
)
{
    assert(mpDifferentiationTree);
    assert(differentiation_type < mpDifferentiationTree->size());

    mpDifferentiationTree = differentiation_tree;
    mDifferentiationType = differentiation_type;
}

void DifferentiationTreeBasedCellCycleModel::Initialise()
{
    assert(mpCell != NULL);

    if (mpDifferentiationTree)
    {
        mSDuration = mpDifferentiationTree->getNode(mDifferentiationType)->getCellCycleLength() / 4.0;
        mG2Duration = mpDifferentiationTree->getNode(mDifferentiationType)->getCellCycleLength() / 4.0;
        mMDuration = mpDifferentiationTree->getNode(mDifferentiationType)->getCellCycleLength() / 4.0;
        mStemCellG1Duration = mpDifferentiationTree->getNode(mDifferentiationType)->getCellCycleLength() / 4.0;
        mTransitCellG1Duration = mpDifferentiationTree->getNode(mDifferentiationType)->getCellCycleLength() / 4.0;

        if (mDifferentiationType == 0)
        {
            boost::shared_ptr<AbstractCellProperty> p_diff_type =
                        mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<StemCellProliferativeType>();
            mpCell->SetCellProliferativeType(p_diff_type);
            mG1Duration = mpDifferentiationTree->getNode(mDifferentiationType)->getCellCycleLength() / 4.0;
            mpCell->GetCellData()->SetItem("Colour", 0.0);
            //mpCell->SetCellProliferativeType(STEM);
        }
        else if (mpDifferentiationTree->getNode(mDifferentiationType)->getNumberOfChildren() == 0)
        {
            boost::shared_ptr<AbstractCellProperty> p_diff_type =
                        mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<DifferentiatedCellProliferativeType>();
            mpCell->SetCellProliferativeType(p_diff_type);
            mG1Duration = DBL_MAX;
            mpCell->GetCellData()->SetItem("Colour", mpDifferentiationTree->getNode(mDifferentiationType)->getColour());
            //mpCell->SetCellProliferativeType(DIFFERENTIATED);
        }
        else
        {
            boost::shared_ptr<AbstractCellProperty> p_diff_type =
                        mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<TransitCellProliferativeType>();
            mpCell->SetCellProliferativeType(p_diff_type);
            mG1Duration = mpDifferentiationTree->getNode(mDifferentiationType)->getCellCycleLength() / 4.0;
            mpCell->GetCellData()->SetItem("Colour", mpDifferentiationTree->getNode(mDifferentiationType)->getColour());
            //mpCell->SetCellProliferativeType(TRANSIT);
        }
    }
    else
    {
        mpCell->GetCellData()->SetItem("Colour", 0.0);
        mG1Duration = 2.0; //? Taken from Meineke model
        mSDuration = 8.0; //? Taken from Meineke model
        mG2Duration = 1.0; //? Taken from Meineke model
        mMDuration = 1.0; //? Taken from Meineke model
        mStemCellG1Duration = 2.0;
        mTransitCellG1Duration = 2.0;
        boost::shared_ptr<AbstractCellProperty> p_diff_type =
                    mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<StemCellProliferativeType>();
        mpCell->SetCellProliferativeType(p_diff_type);
    }
}

void DifferentiationTreeBasedCellCycleModel::InitialiseAfterDivision()
{
    assert(mpCell != NULL);

    if (mpDifferentiationTree)
    {
        mSDuration = mpDifferentiationTree->getNode(mDifferentiationType)->getCellCycleLength() / 4.0;
        mG2Duration = mpDifferentiationTree->getNode(mDifferentiationType)->getCellCycleLength() / 4.0;
        mMDuration = mpDifferentiationTree->getNode(mDifferentiationType)->getCellCycleLength() / 4.0;
        mStemCellG1Duration = mpDifferentiationTree->getNode(mDifferentiationType)->getCellCycleLength() / 4.0;
        mTransitCellG1Duration = mpDifferentiationTree->getNode(mDifferentiationType)->getCellCycleLength() / 4.0;

        if (mDifferentiationType == 0)
        {
            boost::shared_ptr<AbstractCellProperty> p_diff_type =
                        mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<StemCellProliferativeType>();
            mpCell->SetCellProliferativeType(p_diff_type);
            mG1Duration = mpDifferentiationTree->getNode(mDifferentiationType)->getCellCycleLength() / 4.0;
            mpCell->GetCellData()->SetItem("Colour", 0.0);
            //mpCell->SetCellProliferativeType(STEM);
        }
        else if (mpDifferentiationTree->getNode(mDifferentiationType)->getNumberOfChildren() == 0)
        {
            boost::shared_ptr<AbstractCellProperty> p_diff_type =
                        mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<DifferentiatedCellProliferativeType>();
            mpCell->SetCellProliferativeType(p_diff_type);
            mG1Duration = DBL_MAX;
            mpCell->GetCellData()->SetItem("Colour", mpDifferentiationTree->getNode(mDifferentiationType)->getColour());
            //mpCell->SetCellProliferativeType(DIFFERENTIATED);
        }
        else
        {
            boost::shared_ptr<AbstractCellProperty> p_diff_type =
                        mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<TransitCellProliferativeType>();
            mpCell->SetCellProliferativeType(p_diff_type);
            mG1Duration = mpDifferentiationTree->getNode(mDifferentiationType)->getCellCycleLength() / 4.0;
            mpCell->GetCellData()->SetItem("Colour", mpDifferentiationTree->getNode(mDifferentiationType)->getColour());
            //mpCell->SetCellProliferativeType(TRANSIT);
        }
    }
    else
    {
        mpCell->GetCellData()->SetItem("Colour", 0.0);
        mG1Duration = 2.0; //? Taken from Meineke model
        mSDuration = 8.0; //? Taken from Meineke model
        mG2Duration = 1.0; //? Taken from Meineke model
        mMDuration = 1.0; //? Taken from Meineke model
        mStemCellG1Duration = 2.0;
        mTransitCellG1Duration = 2.0;
        boost::shared_ptr<AbstractCellProperty> p_diff_type =
                    mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<StemCellProliferativeType>();
        mpCell->SetCellProliferativeType(p_diff_type);
    }
}

void DifferentiationTreeBasedCellCycleModel::SetG1Duration()
{
    assert(mpCell != NULL);
    if (mpCell->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>())
    {
        mG1Duration = DBL_MAX;
    }
    else
    {
        if (mpDifferentiationTree)
        {
            mG1Duration = mpDifferentiationTree->getNode(mDifferentiationType)->getCellCycleLength() / 4.0;
        }
        else
        {
            mG1Duration = 2.0;
        }
    }
}

AbstractCellCycleModel*  DifferentiationTreeBasedCellCycleModel::CreateCellCycleModel()
{
    DifferentiationTreeBasedCellCycleModel* p_model;
    if (mpDifferentiationTree)
    {
        p_model = new DifferentiationTreeBasedCellCycleModel(
                        mpDifferentiationTree,
                        mDifferentiationType
                    );
    }
    else
    {
        p_model = new DifferentiationTreeBasedCellCycleModel();
    }
    p_model->SetBirthTime(mBirthTime);
    //Differentiation
    InitialiseDaughterCell();

    return p_model;
}

void DifferentiationTreeBasedCellCycleModel::InitialiseDaughterCell()
{
    AbstractSimpleCellCycleModel::InitialiseDaughterCell();
    if (mpDifferentiationTree)
    {
        if (mpDifferentiationTree->getNode(mDifferentiationType)->getNumberOfChildren() != 0)
        {
            //Divide!
            //Important: stationary distribution MUST BE stochastic!
            unsigned child_position = 0;
            std::vector<double> probabilities = mpDifferentiationTree->getNode(mDifferentiationType)->getStationaryDistribution();
            double cumulative_probability = 0.0;
            double random_number = RandomNumberGenerator::Instance()->ranf();
            for (child_position = 0; child_position<probabilities.size()-1; child_position++)
            {
                cumulative_probability += probabilities.at(child_position);
                if (random_number <= cumulative_probability)
                {
                    break;
                }
            }
            mDifferentiationType = mpDifferentiationTree->getNode(mDifferentiationType)->getChildren().at(child_position);
        }
    }
    InitialiseAfterDivision();
}

void DifferentiationTreeBasedCellCycleModel::OutputCellCycleModelParameters(out_stream& rParamsFile)
{
    // No new parameters to output, so just call method on direct parent class
    AbstractSimpleCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(DifferentiationTreeBasedCellCycleModel)
