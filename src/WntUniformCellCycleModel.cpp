

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

    if (WntConcentrationXSection<2>::Instance()->GetWntLevel(mpCell) > WntConcentrationXSection<2>::Instance()->GetWntThreshold())
    {
       AboveThreshold = true;
    }

    return AboveThreshold;
};

// Overloading ReadyToDivide to account for Wnt Concentration
bool WntUniformCellCycleModel::ReadyToDivide()
{
    assert(mpCell != nullptr);
    // Assume that the crypt can be broken into three sections:
    // The niche, where division happens slowly
    // The transient amplifying region where division is rapid
    // The top region where cells have terminally differentiated and stop dividing
    // The point where these regimes change can be controlled by changing the threshold
    if (!mReadyToDivide)
    {
        double wntLevel = WntConcentrationXSection<2>::Instance()->GetWntLevel(mpCell);
        
        if (wntLevel > mTransientRegimeThreshold){
            // Niche division rate
            if (wntLevel > mNicheDivisionRegimeThreshold){
            // Niche division rate
                if (GetAge() >= mNicheCellCycleTime)
                {
                    mReadyToDivide = true;
                }
            } else {
                if (GetAge() >= mTransientCellCycleTime)
                {
                    mReadyToDivide = true;
                }
            }
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