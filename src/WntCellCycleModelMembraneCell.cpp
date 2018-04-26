

#include "WntCellCycleModelMembraneCell.hpp"
#include "AnoikisCellTagged.hpp"
#include "Debug.hpp"

WntCellCycleModelMembraneCell::WntCellCycleModelMembraneCell()
{
};

WntCellCycleModelMembraneCell::~WntCellCycleModelMembraneCell()
{

}

WntCellCycleModelMembraneCell::WntCellCycleModelMembraneCell(const WntCellCycleModelMembraneCell& rModel)
   : mNicheCellCycleTime(rModel.mNicheCellCycleTime),
   mTransientCellCycleTime(rModel.mTransientCellCycleTime)
{
};

AbstractCellCycleModel* WntCellCycleModelMembraneCell::CreateCellCycleModel()
{
    return new WntCellCycleModelMembraneCell(*this);
};


bool WntCellCycleModelMembraneCell::IsAbovetWntThreshold()
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
bool WntCellCycleModelMembraneCell::ReadyToDivide()
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

// Since we are inheriting from UniformCellCycleModel, need to overload this again.
void WntCellCycleModelMembraneCell::ResetForDivision()
{
    assert(mReadyToDivide);
    mReadyToDivide = false;
    mBirthTime = SimulationTime::Instance()->GetTime();
    
}

void WntCellCycleModelMembraneCell::InitialiseDaughterCell()
{
    CellPtr this_cell = GetCell();
    if (this_cell->HasCellProperty<AnoikisCellTagged>())
    {
        TRACE("Removed tag from popped-up daughter")
        this_cell->RemoveCellProperty<AnoikisCellTagged>();
    }
    SetCellCycleTimesForDaughter();
}

void WntCellCycleModelMembraneCell::OutputCellCycleModelParameters(out_stream& rParamsFile)
{

    // Call method on direct parent class
    AbstractCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
};

void WntCellCycleModelMembraneCell::SetNicheCellCycleTime(double nicheCellCycleTime)
{
    mNicheCellCycleTime = nicheCellCycleTime;
    mStoredNicheCellCycleTime = nicheCellCycleTime;

};
void WntCellCycleModelMembraneCell::SetTransientCellCycleTime(double transientCellCycleTime)
{   
    mTransientCellCycleTime = transientCellCycleTime;
    mStoredTransientCellCycleTime = transientCellCycleTime;
};

double WntCellCycleModelMembraneCell::GetAverageTransitCellCycleTime()
{
    return mTransientCellCycleTime;
};

double WntCellCycleModelMembraneCell::GetAverageStemCellCycleTime()
{
    return mNicheCellCycleTime;
};

void WntCellCycleModelMembraneCell::SetCellCycleTimesForDaughter()
{
    //Adds an element of randomness to the cell cycle lengths for daughter cells to avoid synchronising behaviour
    //Can't allow CCT to be less than 2, ideally at least 5  to avoid ridiculously small CCTs
    RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
    double wiggle = 2*p_gen->ranf() - 1; //Number can be +ve or -ve
    
    mTransientCellCycleTime = mStoredTransientCellCycleTime + wiggle;
    mNicheCellCycleTime = mStoredNicheCellCycleTime + wiggle;

};

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(WntCellCycleModelMembraneCell)