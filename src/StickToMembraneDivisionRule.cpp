/*

A division rule to keep the cells both on the membrane

*/

#include "StickToMembraneDivisionRule.hpp"
#include "RandomNumberGenerator.hpp"
#include "AnoikisCellTagged.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::pair<c_vector<double, SPACE_DIM>, c_vector<double, SPACE_DIM> > 
    StickToMembraneDivisionRule<ELEMENT_DIM, SPACE_DIM>::CalculateCellDivisionVector(
        CellPtr pParentCell,
        AbstractCentreBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation)
{
    // Get separation parameter
    double separation = rCellPopulation.GetMeinekeDivisionSeparation();

    c_vector<double, 2> random_vector;
    
    //Still need to allow for division if cell has popped up, so use the anoikis tag to check
    if (pParentCell->HasCellProperty<AnoikisCellTagged>())
    {
        //If its popped up, divide in any direction
        double random_angle = 2.0 * M_PI*RandomNumberGenerator::Instance()->ranf();

        random_vector(0) = 0.5 * separation * cos(random_angle);
        random_vector(1) = 0.5 * separation * sin(random_angle);
    } else {
        //If normal division, split in the direction of membrane axis
        random_vector(0) = 0.5 * separation * mMembraneAxis(0);
        random_vector(1) = 0.5 * separation * mMembraneAxis(1);
        //Need to add in some wiggle to this so that it isn't perfectly in line each time
    }
    
    
    c_vector<double, 2> parent_position = rCellPopulation.GetLocationOfCellCentre(pParentCell) - random_vector;
    c_vector<double, 2> daughter_position = parent_position + random_vector;

    std::pair<c_vector<double, 2>, c_vector<double, 2> > positions(parent_position, daughter_position);

    return positions;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void StickToMembraneDivisionRule<ELEMENT_DIM, SPACE_DIM>::SetMembraneAxis(c_vector<double, 2> membraneAxis)
{
    mMembraneAxis = membraneAxis;
}

// Explicit instantiation
template class StickToMembraneDivisionRule<1,1>;
template class StickToMembraneDivisionRule<1,2>;
template class StickToMembraneDivisionRule<2,2>;
template class StickToMembraneDivisionRule<1,3>;
template class StickToMembraneDivisionRule<2,3>;
template class StickToMembraneDivisionRule<3,3>;

// Serialization for Boost >= 1.36
// #include "SerializationExportWrapperForCpp.hpp"
// EXPORT_TEMPLATE_CLASS_ALL_DIMS(StickToMembraneDivisionRule)
