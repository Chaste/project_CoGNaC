#ifndef DifferentiationTreeBasedWithAsymmetricDivisionCellCycleModel_HPP_
#define DifferentiationTreeBasedWithAsymmetricDivisionCellCycleModel_HPP_

#include "DifferentiationTreeNode.hpp"
#include "DifferentiationTree.hpp"
#include "AbstractSimpleCellCycleModel.hpp"
#include "RandomNumberGenerator.hpp"
//friend class TestSimpleCellCycleModels;


class DifferentiationTreeBasedWithAsymmetricDivisionCellCycleModel : public AbstractSimpleCellCycleModel
{
private:
    /** Pointer to the DifferentiatioTree */
    const DifferentiationTree* mpDifferentiationTree;

    /** Index to the ID of current node of the tree */
    unsigned mDifferentiationType;

    /** Needed for serialization. */
    friend class boost::serialization::access;

    /**
     * Archive the cell-cycle model and random number generator, never used directly - boost uses this.
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
        archive & p_gen;
    }

public:
    DifferentiationTreeBasedWithAsymmetricDivisionCellCycleModel();
    DifferentiationTreeBasedWithAsymmetricDivisionCellCycleModel(
            const DifferentiationTree* differentiation_tree,
            unsigned differentiation_type
            );

    /**
     * The destruptor doesn't delete the object. Just put they to NULL.
     */
    virtual ~DifferentiationTreeBasedWithAsymmetricDivisionCellCycleModel();

    /**
     * Overridden SetG1Duration Method
     */
    void SetG1Duration();

    /**
     * Setter for mpDifferentiationTree and mpDifferentiationType
     */
    void SetTreeAndType(const DifferentiationTree* differentiation_tree,
            unsigned differentiation_type);

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
    void InitialiseDaughterCell();

    /**
     * Check and perform division in the parent cell.
     */
    void InitialiseParentCell();

    /**
     * Overridden builder method to create new copies of
     * this cell-cycle model.
     *
     * @return new cell-cycle model
     */
    AbstractCellCycleModel*  CreateCellCycleModel();

    /**
     * Outputs cell cycle model parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputCellCycleModelParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(DifferentiationTreeBasedWithAsymmetricDivisionCellCycleModel)

#endif /* DifferentiationTreeBasedWithAsymmetricDivisionCellCycleModel_HPP_ */
