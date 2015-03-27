#include "DifferentiationTreeBasedWntCellCycleModel.hpp"
#include "StemCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "VoronoiDataWriter.hpp"
#include <iostream>

DifferentiationTreeBasedWntCellCycleModel::DifferentiationTreeBasedWntCellCycleModel()
     : mpDifferentiationTree(NULL),
       mDifferentiationType(0),
       mStemType(0),
       mPanethType(0),
       mTa1Type(0),
       mTa2AType(0),
       mTa2BType(0),
       mEnterocyteType(0),
       mEnteroendocrineType(0),
       mGobletType(0),
       mUseCellProliferativeTypeDependentG1Duration(false),
       mWntPanethThreshold(.9),//.9
       mWntTa1Threshold(.4),//.6
       mWntTa2Threshold(.2)//.4
{
    mG1Duration = 2.0; //? Taken from Meineke model
    mSDuration = 8.0; //? Taken from Meineke model
    mG2Duration = 1.0; //? Taken from Meineke model
    mMDuration = 1.0; //? Taken from Meineke model
    mStemCellG1Duration = 2.0;
    mTransitCellG1Duration = 2.0;
    mWntStemThreshold = .8;//.8
    mWntTransitThreshold = .2;//.5
    mWntLabelledThreshold = .3;//.5
}

DifferentiationTreeBasedWntCellCycleModel::DifferentiationTreeBasedWntCellCycleModel
(
        const DifferentiationTree* differentiation_tree,
        unsigned differentiation_type
)
:       mpDifferentiationTree(differentiation_tree),
        mDifferentiationType(differentiation_type),
        mStemType(0),
        mPanethType(0),
        mTa1Type(0),
        mTa2AType(0),
        mTa2BType(0),
        mEnterocyteType(0),
        mEnteroendocrineType(0),
        mGobletType(0),
        mUseCellProliferativeTypeDependentG1Duration(false),
        mWntPanethThreshold(.9),//.9
        mWntTa1Threshold(.4),//.6
        mWntTa2Threshold(.2)//.4
{
    assert(differentiation_tree);
    assert(differentiation_type < differentiation_tree->size());

    mWntStemThreshold = .8;//.8
    mWntTransitThreshold = .2;//.5
    mWntLabelledThreshold = .3;//.5
}

DifferentiationTreeBasedWntCellCycleModel::~DifferentiationTreeBasedWntCellCycleModel()
{
    mpDifferentiationTree = NULL;
}

void DifferentiationTreeBasedWntCellCycleModel::SetTreeAndType
(
        const DifferentiationTree* differentiation_tree,
        unsigned differentiation_type
)
{
    assert(differentiation_tree);
    assert(differentiation_type < differentiation_tree->size());

    mpDifferentiationTree = differentiation_tree;
    mDifferentiationType = differentiation_type;
}

void DifferentiationTreeBasedWntCellCycleModel::Initialise()
{
    assert(mpCell != NULL);

    if (mpDifferentiationTree)
    {
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
        assert(tree->topologyTreeCompare(mpDifferentiationTree));
        delete tree;

        unsigned child = mpDifferentiationTree->getRoot()->getChildren().at(0);

        if (mpDifferentiationTree->getNode(child)->getNumberOfChildren() == 0)
        {
            mPanethType = child;
            mTa1Type = mpDifferentiationTree->getRoot()->getChildren().at(1);
        }else
        {
            mPanethType = mpDifferentiationTree->getRoot()->getChildren().at(1);
            mTa1Type = child;
        }
        child = mpDifferentiationTree->getNode(mTa1Type)->getChildren().at(0);
        if (mpDifferentiationTree->getNode(child)->getNumberOfChildren()==2)
        {
            mTa2AType = child;
            mTa2BType = mpDifferentiationTree->getNode(mTa1Type)->getChildren().at(1);
        }
        else
        {
            mTa2AType = mpDifferentiationTree->getNode(mTa1Type)->getChildren().at(1);
            mTa2BType = child;
        }
        mEnterocyteType = mpDifferentiationTree->getNode(mTa2AType)->getChildren().at(0);
        mEnteroendocrineType = mpDifferentiationTree->getNode(mTa2AType)->getChildren().at(1);
        mGobletType = mpDifferentiationTree->getNode(mTa2BType)->getChildren().at(0);

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

void DifferentiationTreeBasedWntCellCycleModel::InitialiseAfterDivision()
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

AbstractCellCycleModel*  DifferentiationTreeBasedWntCellCycleModel::CreateCellCycleModel()
{
    DifferentiationTreeBasedWntCellCycleModel* p_model;
    if (mpDifferentiationTree)
    {
        p_model = new DifferentiationTreeBasedWntCellCycleModel(
                        mpDifferentiationTree,
                        mDifferentiationType
                    );
    }
    else
    {
        p_model = new DifferentiationTreeBasedWntCellCycleModel();
    }
    p_model->SetBirthTime(mBirthTime);
    p_model->SetDimension(mDimension);
    p_model->SetMinimumGapDuration(mMinimumGapDuration);
    p_model->SetStemCellG1Duration(mStemCellG1Duration);
    p_model->SetTransitCellG1Duration(mTransitCellG1Duration);
    p_model->SetSDuration(mSDuration);
    p_model->SetG2Duration(mG2Duration);
    p_model->SetMDuration(mMDuration);
    p_model->SetUseCellProliferativeTypeDependentG1Duration(mUseCellProliferativeTypeDependentG1Duration);
    p_model->SetWntStemThreshold(mWntStemThreshold);
    p_model->SetWntTransitThreshold(mWntTransitThreshold);
    p_model->SetWntLabelledThreshold(mWntLabelledThreshold);
    p_model->SetWntPanethThreshold(mWntPanethThreshold);
    p_model->SetWntTa1Threshold(mWntTa1Threshold);
    p_model->SetWntTa2Threshold(mWntTa2Threshold);
    p_model->SetPanethType(mPanethType);
    p_model->SetTa1Type(mTa1Type);
    p_model->SetTa2AType(mTa2AType);
    p_model->SetTa2BType(mTa2BType);
    p_model->SetEnterocyteType(mEnterocyteType);
    p_model->SetEnteroendocrineType(mEnteroendocrineType);
    p_model->SetGobletType(mGobletType);

    //Differentiation for child (the "new" cell)
    //InitialiseDaughterCell();

    return p_model;
}

void DifferentiationTreeBasedWntCellCycleModel::SetUseCellProliferativeTypeDependentG1Duration(bool useCellProliferativeTypeDependentG1Duration)
{
    mUseCellProliferativeTypeDependentG1Duration = useCellProliferativeTypeDependentG1Duration;
}

void DifferentiationTreeBasedWntCellCycleModel::SetG1Duration()
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
    if (mG1Duration < mMinimumGapDuration)
    {
        mG1Duration = mMinimumGapDuration;
    }
}

double DifferentiationTreeBasedWntCellCycleModel::GetWntLevel()
{
    assert(mpCell != NULL);
    double level = 0;

    switch (mDimension)
    {
        case 1:
        {
            const unsigned DIM = 1;
            level = WntConcentration<DIM>::Instance()->GetWntLevel(mpCell);
            break;
        }
        case 2:
        {
            const unsigned DIM = 2;
            level = WntConcentration<DIM>::Instance()->GetWntLevel(mpCell);
            break;
        }
        case 3:
        {
            const unsigned DIM = 3;
            level = WntConcentration<DIM>::Instance()->GetWntLevel(mpCell);
            break;
        }
        default:
            NEVER_REACHED;
    }
    return level;
}

WntConcentrationType DifferentiationTreeBasedWntCellCycleModel::GetWntType()
{
    WntConcentrationType wnt_type;
    switch (mDimension)
    {
        case 1:
        {
            const unsigned DIM = 1;
            wnt_type = WntConcentration<DIM>::Instance()->GetType();
            break;
        }
        case 2:
        {
            const unsigned DIM = 2;
            wnt_type = WntConcentration<DIM>::Instance()->GetType();
            break;
        }
        case 3:
        {
            const unsigned DIM = 3;
            wnt_type = WntConcentration<DIM>::Instance()->GetType();
            break;
        }
        default:
            NEVER_REACHED;
    }
    return wnt_type;
}

void DifferentiationTreeBasedWntCellCycleModel::UpdateCellCyclePhase()
{
    // The cell can divide if the Wnt concentration >= wnt_division_threshold
    double wnt_division_threshold = DBL_MAX;

    // Set up under what level of Wnt stimulus a cell will divide
    if (mpCell->GetMutationState()->IsType<WildTypeCellMutationState>())
    {
        wnt_division_threshold = mWntTa2Threshold;
    }
    else if (mpCell->GetMutationState()->IsType<ApcOneHitCellMutationState>())
    {
        // should be less than healthy values
        wnt_division_threshold = 0.77*mWntTa2Threshold;
    }
    else if (mpCell->GetMutationState()->IsType<BetaCateninOneHitCellMutationState>())
    {
        // less than above value
        wnt_division_threshold = 0.155*mWntTa2Threshold;
    }
    else if (mpCell->GetMutationState()->IsType<ApcTwoHitCellMutationState>())
    {
        // should be zero (no Wnt-dependence)
        wnt_division_threshold = 0.0;
    }
    else
    {
        NEVER_REACHED;
    }

    if (mpCell->HasCellProperty<CellLabel>())
    {
        wnt_division_threshold = mWntLabelledThreshold;
    }

    double wnt_level = GetWntLevel();
    //WntConcentrationType wnt_type = GetWntType();

    // Set the cell type to TransitCellProliferativeType if the Wnt stimulus exceeds wnt_division_threshold
    if (wnt_level >= wnt_division_threshold)
    {
        /*
         * This method is usually called within a CellBasedSimulation, after the CellPopulation
         * has called CellPropertyRegistry::TakeOwnership(). This means that were we to call
         * CellPropertyRegistry::Instance() here when setting the CellProliferativeType, we
         * would be creating a new CellPropertyRegistry. In this case the cell proliferative
         * type counts, as returned by AbstractCellPopulation::GetCellProliferativeTypeCount(),
         * would be incorrect. We must therefore access the CellProliferativeType via the cell's
         * CellPropertyCollection.
         */

        if (wnt_level >= mWntPanethThreshold)
        {
            if (mDifferentiationType != mPanethType)
            {
                mDifferentiationType = mPanethType;
                InitialiseAfterDivision();
            }
        } else if (wnt_level >= mWntStemThreshold)
        {
            if (mDifferentiationType != mStemType)
            {
                mDifferentiationType = mStemType;
                InitialiseAfterDivision();
            }
        } else if (wnt_level >= mWntTa1Threshold)
        {
            if (mDifferentiationType != mTa1Type)
            {
                mDifferentiationType = mTa1Type;
                InitialiseAfterDivision();
            }
        } else if (wnt_level >= mWntTa2Threshold)
        {
            if (mDifferentiationType != mTa2AType &&
                mDifferentiationType != mTa2BType)
            {
                if (mDifferentiationType != mTa1Type)
                {
                    mDifferentiationType = mTa1Type;
                    //InitialiseAfterDivision();
                }
                bool is_ta2a_first_child = true;
                if (mpDifferentiationTree->getNode(mDifferentiationType)->getChildren().at(0) == mTa2BType)
                {
                    is_ta2a_first_child = false;
                }
                double probability_first_child = mpDifferentiationTree->getNode(mDifferentiationType)->getStationaryDistribution().at(0);
                if (RandomNumberGenerator::Instance()->ranf() <= probability_first_child)
                {
                    if (is_ta2a_first_child)
                        mDifferentiationType = mTa2AType;
                    else
                        mDifferentiationType = mTa2BType;
                }
                else
                {
                    if (is_ta2a_first_child)
                         mDifferentiationType = mTa2BType;
                    else
                         mDifferentiationType = mTa2AType;
                }
                InitialiseAfterDivision();
            }
        }
    }
    else
    {
        if (mDifferentiationType != mEnterocyteType &&
            mDifferentiationType != mEnteroendocrineType &&
            mDifferentiationType != mGobletType)
        {
            if (mDifferentiationType != mTa2AType &&
                mDifferentiationType != mTa2BType)
            {
                if (mDifferentiationType != mTa1Type)
                {
                    mDifferentiationType = mTa1Type;
                    //InitialiseAfterDivision();
                }
                bool is_ta2a_first_child = true;
                if (mpDifferentiationTree->getNode(mDifferentiationType)->getChildren().at(0) == mTa2BType)
                {
                    is_ta2a_first_child = false;
                }
                double probability_first_child = mpDifferentiationTree->getNode(mDifferentiationType)->getStationaryDistribution().at(0);
                if (RandomNumberGenerator::Instance()->ranf() <= probability_first_child)
                {
                    if (is_ta2a_first_child)
                        mDifferentiationType = mTa2AType;
                    else
                        mDifferentiationType = mTa2BType;
                }
                else
                {
                    if (is_ta2a_first_child)
                         mDifferentiationType = mTa2BType;
                    else
                         mDifferentiationType = mTa2AType;
                }
            }
            if (mDifferentiationType == mTa2BType)
            {
                mDifferentiationType = mGobletType;
            }
            else
            {
                double probability_enterocyte = mpDifferentiationTree->getNode(mDifferentiationType)->getStationaryDistribution().at(0);
                if (RandomNumberGenerator::Instance()->ranf() <= probability_enterocyte)
                {
                    mDifferentiationType = mEnterocyteType;
                }
                else
                {
                    mDifferentiationType = mEnteroendocrineType;
                }
            }
            InitialiseAfterDivision();
        }
    }
    AbstractSimpleCellCycleModel::UpdateCellCyclePhase();
}

void DifferentiationTreeBasedWntCellCycleModel::InitialiseDaughterCell()
{
    AbstractSimpleCellCycleModel::InitialiseDaughterCell();
    //InitialiseAfterDivision();
}

bool DifferentiationTreeBasedWntCellCycleModel::CanCellTerminallyDifferentiate()
{
    return false;
}

double DifferentiationTreeBasedWntCellCycleModel::GetWntStemThreshold()
{
    return mWntStemThreshold;
}

void DifferentiationTreeBasedWntCellCycleModel::SetWntStemThreshold(double wntStemThreshold)
{
    assert(wntStemThreshold <= 1.0);
    assert(wntStemThreshold >= 0.0);
    mWntStemThreshold = wntStemThreshold;
}

double DifferentiationTreeBasedWntCellCycleModel::GetWntTransitThreshold()
{
    return mWntTransitThreshold;
}

void DifferentiationTreeBasedWntCellCycleModel::SetWntTransitThreshold(double wntTransitThreshold)
{
    //assert(wntTransitThreshold <= 1.0);
    //assert(wntTransitThreshold >= 0.0);
    mWntTransitThreshold = wntTransitThreshold;
}

double DifferentiationTreeBasedWntCellCycleModel::GetWntLabelledThreshold()
{
    return mWntLabelledThreshold;
}

void DifferentiationTreeBasedWntCellCycleModel::SetWntLabelledThreshold(double wntLabelledThreshold)
{
//    assert(wntLabelledThreshold <= 1.0);
//    assert(wntLabelledThreshold >= 0.0);
    mWntLabelledThreshold = wntLabelledThreshold;
}

double DifferentiationTreeBasedWntCellCycleModel::GetWntPanethThreshold()
{
    return mWntPanethThreshold;
}

void DifferentiationTreeBasedWntCellCycleModel::SetWntPanethThreshold(double wntPanethThreshold)
{
    mWntTransitThreshold = wntPanethThreshold;
}

double DifferentiationTreeBasedWntCellCycleModel::GetWntTa1Threshold()
{
    return mWntTa1Threshold;
}

void DifferentiationTreeBasedWntCellCycleModel::SetWntTa1Threshold(double wntTa1Threshold)
{
    mWntTa1Threshold = wntTa1Threshold;
}

double DifferentiationTreeBasedWntCellCycleModel::GetWntTa2Threshold()
{
    return mWntTa2Threshold;
}

void DifferentiationTreeBasedWntCellCycleModel::SetWntTa2Threshold(double wntTa2Threshold)
{
    mWntTa2Threshold = wntTa2Threshold;
}

void DifferentiationTreeBasedWntCellCycleModel::SetPanethType(unsigned paneth)
{
    assert(paneth > 0 && paneth < mpDifferentiationTree->size());
    mPanethType = paneth;
}

void DifferentiationTreeBasedWntCellCycleModel::SetTa1Type(unsigned ta1)
{
    assert(ta1 > 0 && ta1 < mpDifferentiationTree->size());
    mTa1Type = ta1;
}
void DifferentiationTreeBasedWntCellCycleModel::SetTa2AType(unsigned ta2a)
{
    assert(ta2a > 0 && ta2a < mpDifferentiationTree->size());
    mTa2AType = ta2a;
}

void DifferentiationTreeBasedWntCellCycleModel::SetTa2BType(unsigned ta2b)
{
    assert(ta2b > 0 && ta2b < mpDifferentiationTree->size());
    mTa2BType = ta2b;
}

void DifferentiationTreeBasedWntCellCycleModel::SetEnterocyteType(unsigned enterocyte)
{
    assert(enterocyte > 0 && enterocyte < mpDifferentiationTree->size());
    mEnterocyteType = enterocyte;
}

void DifferentiationTreeBasedWntCellCycleModel::SetEnteroendocrineType(unsigned enteroendocrine)
{
    assert(enteroendocrine > 0 && enteroendocrine < mpDifferentiationTree->size());
    mEnteroendocrineType = enteroendocrine;
}

void DifferentiationTreeBasedWntCellCycleModel::SetGobletType(unsigned goblet)
{
    assert(goblet > 0 && goblet < mpDifferentiationTree->size());
    mGobletType = goblet;
}

void DifferentiationTreeBasedWntCellCycleModel::OutputCellCycleModelParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<UseCellProliferativeTypeDependentG1Duration>" << mUseCellProliferativeTypeDependentG1Duration << "</UseCellProliferativeTypeDependentG1Duration>\n";
    *rParamsFile << "\t\t\t<WntStemThreshold>" << mWntStemThreshold << "</WntStemThreshold>\n";
    *rParamsFile << "\t\t\t<WntTransitThreshold>" << mWntTransitThreshold << "</WntTransitThreshold>\n";
    *rParamsFile << "\t\t\t<WntLabelledThreshold>" << mWntLabelledThreshold << "</WntLabelledThreshold>\n";
    *rParamsFile << "\t\t\t<WntPanethThreshold>" << mWntPanethThreshold << "</WntPanethThreshold>\n";
    *rParamsFile << "\t\t\t<WntTa1Threshold>" << mWntTa1Threshold << "</WntTa1Threshold>\n";
    *rParamsFile << "\t\t\t<WntTa2Threshold>" << mWntTa2Threshold << "</WntTa2Threshold>\n";
    *rParamsFile << "\t\t\t<StemType>" << mStemType << "</StemType>\n";
    *rParamsFile << "\t\t\t<PanethType>" << mPanethType << "</PanethType>\n";
    *rParamsFile << "\t\t\t<Ta1Type>" << mTa1Type << "</Ta1Type>\n";
    *rParamsFile << "\t\t\t<Ta2AType>" << mTa2AType << "</Ta2AType>\n";
    *rParamsFile << "\t\t\t<Ta2BType>" << mTa2BType << "</Ta2BType>\n";
    *rParamsFile << "\t\t\t<EnterocyteType>" << mEnterocyteType << "</EnterocyteType>\n";
    *rParamsFile << "\t\t\t<EnteroendocrineType>" << mEnteroendocrineType << "</EnteroendocrineType>\n";
    *rParamsFile << "\t\t\t<GobletType>" << mGobletType << "</GobletType>\n";

    // No new parameters to output, so just call method on direct parent class
    AbstractSimpleCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(DifferentiationTreeBasedWntCellCycleModel)
