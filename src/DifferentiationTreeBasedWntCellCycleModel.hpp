#ifndef DIFFERENTIATIONTREEBASEDWNTCELLCYCLEMODEL_HPP_
#define DIFFERENTIATIONTREEBASEDWNTCELLCYCLEMODEL_HPP_

#include "DifferentiationTreeNode.hpp"
#include "DifferentiationTree.hpp"
#include "AbstractSimpleCellCycleModel.hpp"
#include "RandomNumberGenerator.hpp"
#include "WntConcentration.hpp"

// Needed here to avoid serialization errors
#include "WildTypeCellMutationState.hpp"
#include "CellLabel.hpp"
#include "ApcOneHitCellMutationState.hpp"
#include "ApcTwoHitCellMutationState.hpp"
#include "BetaCateninOneHitCellMutationState.hpp"

class DifferentiationTreeBasedWntCellCycleModel : public AbstractSimpleCellCycleModel
{
private:
    /** Pointer to the DifferentiatioTree */
    const DifferentiationTree* mpDifferentiationTree;

    /** Index to the ID of current node of the tree */
    unsigned mDifferentiationType;

    unsigned mStemType;
    unsigned mPanethType;
    unsigned mTa1Type;
    unsigned mTa2AType;
    unsigned mTa2BType;
    unsigned mEnterocyteType;
    unsigned mEnteroendocrineType;
    unsigned mGobletType;

    /** Needed for serialization. */
    friend class boost::serialization::access;

    /**
     * Archive the cell-cycle model, never used directly - boost uses this.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractSimpleCellCycleModel>(*this);

        RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
        archive & *p_gen;

        archive & mUseCellProliferativeTypeDependentG1Duration;
        archive & mWntStemThreshold;
        archive & mWntTransitThreshold;
        archive & mWntPanethThreshold;
        archive & mWntTa1Threshold;
        archive & mWntTa2Threshold;
        archive & mWntLabelledThreshold;
        //archive & mpDifferentiationTree;
        archive & mDifferentiationType;
        archive & mStemType;
        archive & mPanethType;
        archive & mTa1Type;
        archive & mTa2AType;
        archive & mTa2BType;
        archive & mEnterocyteType;
        archive & mEnteroendocrineType;
        archive & mGobletType;

    }

protected:

    /**
     * Whether to use different mean G1 durations for different cell types.
     * For use in SetG1Duration().
     */
    bool mUseCellProliferativeTypeDependentG1Duration;

    /**
     * Non-dimensionalized Wnt threshold, above which cells behave as stem cells.
     */
    double mWntStemThreshold;

    /**
     * Non-dimensionalized Wnt threshold, above which cells progress through the cell cycle.
     */
    double mWntTransitThreshold;

    /**
     * Non-dimensionalized Wnt threshold, above which labelled cells progress through the cell cycle.
     */
    double mWntLabelledThreshold;

    /**
     * @return the Wnt level experienced by the cell.
     */
    double GetWntLevel();

    /**
     * @return the type of Wnt concentration (LINEAR, RADIAL, EXPONENTIAL or NONE).
     * This affects how the cell cycle phase is updated.
     */
    WntConcentrationType GetWntType();

    /**
     * Non-dimensionalized Wnt threshold, above which cells behave as Paneth cells.
     */
    double mWntPanethThreshold;

    /**
     * Non-dimensionalized Wnt threshold, above which cells behave as TA-1 cells.
     */
    double mWntTa1Threshold;

    /**
     * Non-dimensionalized Wnt threshold, above which cells behave as TA-2 cells.
     */
    double mWntTa2Threshold;

public:
    DifferentiationTreeBasedWntCellCycleModel();
    DifferentiationTreeBasedWntCellCycleModel(
            const DifferentiationTree* differentiation_tree,
            unsigned differentiation_type
            );

    /**
     * The destruptor doesn't delete the object. Just put they to NULL.
     */
    virtual ~DifferentiationTreeBasedWntCellCycleModel();

    /**
     * Overridden SetG1Duration Method
     */
    void SetG1Duration();

    /**
     * Setter for mpDifferentiationTree and mpDifferentiationType
     */
    void SetTreeAndType(const DifferentiationTree* differentiation_tree,
            unsigned diffreentiation_type);

    /**
     * Initialise the duration based on the DifferentaitionType
     */
    void Initialise();

    /**
     * Re-Initialise the duration based on the new DifferentaitionType
     */
    void InitialiseAfterDivision();

    /**
     * Set the new cell type based on the new DifferentiationTree once
     * it has been created after division.
     * The duration of the cell cycle length will be based on cell type.
     */
    virtual void InitialiseDaughterCell();

    /**
     * Overridden UpdateCellCyclePhase() method.
     */
    virtual void UpdateCellCyclePhase();
    /**
     * Overridden builder method to create new copies of
     * this cell-cycle model.
     *
     * @return new cell-cycle model
     */
    virtual AbstractCellCycleModel*  CreateCellCycleModel();

    /**
     * Set whether Whether the duration of the G1 phase is dependent on cell type
     * @param useCellProliferativeTypeDependentG1Duration - boolean, defaults to true.
     */
    void SetUseCellProliferativeTypeDependentG1Duration(bool useCellProliferativeTypeDependentG1Duration=true);

    /**
     * Overridden CanCellTerminallyDifferentiate() method.
     * @return whether cell can terminally differentiate
     */
    virtual bool CanCellTerminallyDifferentiate();

    /**
     * @return mWntStemThreshold
     */
    double GetWntStemThreshold();

    /**
     * Set mWntStemThreshold.
     *
     * @param wntStemThreshold the value of mWntStemThreshold
     */
    void SetWntStemThreshold(double wntStemThreshold);

    /**
     * @return mWntTransitThreshold
     */
    double GetWntTransitThreshold();

    /**
     * Set mWntTransitThreshold.
     *
     * @param wntTransitThreshold the value of mWntTransitThreshold
     */
    void SetWntTransitThreshold(double wntTransitThreshold);

    /**
     * @return mWntLabelledThreshold
     */
    double GetWntLabelledThreshold();

    /**
     * Set mWntLabelledThreshold.
     *
     * @param wntLabelledThreshold the value of mWntLabelledThreshold
     */
    void SetWntLabelledThreshold(double wntLabelledThreshold);

    /**
     * @return mWntPanethThreshold
     */
    double GetWntPanethThreshold();

    /**
     * Set mWntPanethThreshold.
     *
     * @param wntPanethThreshold the value of mWntPanethThreshold
     */
    void SetWntPanethThreshold(double wntPanethThreshold);

    /**
     * @return mWntTa1Threshold
     */
    double GetWntTa1Threshold();

    /**
     * Set mWntTa1Threshold.
     *
     * @param wntTa1Threshold the value of mWntTa1Threshold
     */
    void SetWntTa1Threshold(double wntTa1Threshold);

    /**
     * @return mWntTa2Threshold
     */
    double GetWntTa2Threshold();

    /**
     * Set mWntTa2Threshold.
     *
     * @param wntTa2Threshold the value of mWntTa2Threshold
     */
    void SetWntTa2Threshold(double wntTa2Threshold);

    void SetPanethType(unsigned paneth);
    void SetTa1Type(unsigned ta1);
    void SetTa2AType(unsigned ta2a);
    void SetTa2BType(unsigned ta2b);
    void SetEnterocyteType(unsigned enterocyte);
    void SetEnteroendocrineType(unsigned enteroendocrine);
    void SetGobletType(unsigned goblet);

    /**
     * Outputs cell cycle model parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputCellCycleModelParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(DifferentiationTreeBasedWntCellCycleModel)

#endif /* DIFFERENTIATIONTREEBASEDWNTCELLCYCLEMODEL_HPP_ */
