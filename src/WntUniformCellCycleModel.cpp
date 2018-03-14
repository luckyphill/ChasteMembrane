

#include "WntUniformCellCycleModel.hpp"
#include "Debug.hpp"

WntUniformCellCycleModel::WntUniformCellCycleModel()
    : UniformCellCycleModel()
{
};

WntUniformCellCycleModel::~WntUniformCellCycleModel()
{

}

WntUniformCellCycleModel::WntUniformCellCycleModel(const WntUniformCellCycleModel& rModel)
   : UniformCellCycleModel(rModel)
{
    /*
     * Initialize only those member variables defined in this class.
     *
     * The member variable mCellCycleDuration is initialized in the
     * AbstractSimpleCellCycleModel constructor.
     *
     * The member variables mBirthTime, mReadyToDivide and mDimension
     * are initialized in the AbstractCellCycleModel constructor.
     *
     * Note that mCellCycleDuration is (re)set as soon as
     * InitialiseDaughterCell() is called on the new cell-cycle model.
     */
};

AbstractCellCycleModel* WntUniformCellCycleModel::CreateCellCycleModel()
{
    return new WntUniformCellCycleModel(*this);
};

bool WntUniformCellCycleModel::IsAbovetWntThreshold()
{
    assert(mpCell != nullptr);
    double level = 0;
    bool AboveThreshold = false;
    mDimension = 2;
    TRACE("Manually set dimension in line 43 of WntUniformCellCycleModel.cpp. Need to work how to do this properly")
    switch (mDimension)
    {
        case 1:
        {
            const unsigned DIM = 1;
            if (WntConcentrationXSection<DIM>::Instance()->GetWntLevel(mpCell) > WntConcentrationXSection<DIM>::Instance()->GetWntThreshold())
            {
               AboveThreshold = true;
            }
            break;
        }
        case 2:
        {
            const unsigned DIM = 2;
            if (WntConcentrationXSection<DIM>::Instance()->GetWntLevel(mpCell) > WntConcentrationXSection<DIM>::Instance()->GetWntThreshold())
            {
               AboveThreshold = true;
            }
            break;
        }
        case 3:
        {
            const unsigned DIM = 3;
            if (WntConcentrationXSection<DIM>::Instance()->GetWntLevel(mpCell) > WntConcentrationXSection<DIM>::Instance()->GetWntThreshold())
            {
               AboveThreshold = true;
            }
            break;
        }
        default:
            NEVER_REACHED;
    }

    return AboveThreshold;
};

// Overloading ReadyToDivide to account for Wnt Concentration
bool WntUniformCellCycleModel::ReadyToDivide()
{
    assert(mpCell != nullptr);

    if (!mReadyToDivide)
    {
        if (GetAge() >= mCellCycleDuration && IsAbovetWntThreshold())
        {
            mReadyToDivide = true;
        }
    }
    return mReadyToDivide;
};

void WntUniformCellCycleModel::OutputCellCycleModelParameters(out_stream& rParamsFile)
{

    // Call method on direct parent class
    AbstractSimpleCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
};

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(WntUniformCellCycleModel)