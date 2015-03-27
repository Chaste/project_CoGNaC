#ifndef DIFFERENTIATIONTREEBASEDCELLCYCLEMODEL_HPP_
#define DIFFERENTIATIONTREEBASEDCELLCYCLEMODEL_HPP_

#include "DifferentiationTreeNode.hpp"
#include "DifferentiationTree.hpp"
#include "AbstractSimpleCellCycleModel.hpp"
#include "RandomNumberGenerator.hpp"
//friend class TestSimpleCellCycleModels;


class DifferentiationTreeBasedCellCycleModel : public AbstractSimpleCellCycleModel
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
    DifferentiationTreeBasedCellCycleModel();
    DifferentiationTreeBasedCellCycleModel(
            const DifferentiationTree* differentiation_tree,
            unsigned differentiation_type
            );

    /**
     * The destruptor doesn't delete the object. Just put they to NULL.
     */
    virtual ~DifferentiationTreeBasedCellCycleModel();

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
CHASTE_CLASS_EXPORT(DifferentiationTreeBasedCellCycleModel)

#endif /* DIFFERENTIATIONTREEBASEDCELLCYCLEMODEL_HPP_ */
