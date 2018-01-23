/*
Initial strucure borrows heavily from EpithelialLayerBasementMembraneForce by Axel Almet
- Added mutation that turns a differentiated cell into a "membrane cell" in 
in order to test a method of introducing a membrane
- The modifications here only change the way the "mutant" cells interact with each
other. Otherwise they are still considered "differentiated" cells for other interactions
*/

#include "IsNan.hpp"
#include "AbstractCellProperty.hpp"

#include "MembraneCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "StemCellProliferativeType.hpp"

#include "Debug.hpp"

#include "LinearSpringForceMembraneCellNodeBased.hpp"


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
LinearSpringForceMembraneCellNodeBased<ELEMENT_DIM,SPACE_DIM>::LinearSpringForceMembraneCellNodeBased()
   : AbstractTwoBodyInteractionForce<ELEMENT_DIM,SPACE_DIM>(),
    mEpithelialSpringStiffness(15.0), // Epithelial covers stem and transit
    mMembraneSpringStiffness(15.0),
    mStromalSpringStiffness(15.0), // Stromal is the differentiated "filler" cells
    mEpithelialMembraneSpringStiffness(15.0),
    mMembraneStromalSpringStiffness(15.0),
    mStromalEpithelialSpringStiffness(15.0),
    mEpithelialRestLength(1.0),
    mMembraneRestLength(1.0),
    mStromalRestLength(1.0),
    mEpithelialMembraneRestLength(1.0),
    mMembraneStromalRestLength(0.6),
    mStromalEpithelialRestLength(1.0),
    mEpithelialCutOffRadius(1.5), // Epithelial covers stem and transit
    mMembraneCutOffRadius(1.5),
    mStromalCutOffRadius(1.5), // Stromal is the differentiated "filler" cells
    mEpithelialMembraneCutOffRadius(1.5),
    mMembraneStromalCutOffRadius(1.5),
    mStromalEpithelialCutOffRadius(1.5)
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
LinearSpringForceMembraneCellNodeBased<ELEMENT_DIM,SPACE_DIM>::~LinearSpringForceMembraneCellNodeBased()
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> LinearSpringForceMembraneCellNodeBased<ELEMENT_DIM,SPACE_DIM>::CalculateForceBetweenNodes(unsigned nodeAGlobalIndex,
                                                                                    unsigned nodeBGlobalIndex,
                                                                                    AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation)
{
    // We should only ever calculate the force between two distinct nodes
    assert(nodeAGlobalIndex != nodeBGlobalIndex);

    Node<SPACE_DIM>* p_node_a = rCellPopulation.GetNode(nodeAGlobalIndex);
    Node<SPACE_DIM>* p_node_b = rCellPopulation.GetNode(nodeBGlobalIndex);

    // Get the node locations
    c_vector<double, SPACE_DIM> node_a_location = p_node_a->rGetLocation();
    c_vector<double, SPACE_DIM> node_b_location = p_node_b->rGetLocation();

    double node_a_radius = 0.0;
    double node_b_radius = 0.0;

    // Get the unit vector parallel to the line joining the two nodes
    c_vector<double, SPACE_DIM> unitForceDirection;

    unitForceDirection = rCellPopulation.rGetMesh().GetVectorFromAtoB(node_a_location, node_b_location);

    // Calculate the distance between the two nodes
    double distance_between_nodes = norm_2(unitForceDirection);
    assert(distance_between_nodes > 0);
    assert(!std::isnan(distance_between_nodes));

    unitForceDirection /= distance_between_nodes;

    /*
     * Calculate the rest length of the spring connecting the two nodes with a default
     * value of 1.0.
     */

    // We have three types of cells, with 6 different possible pairings as demarked by the 6 different spring stiffnesses
    // Need to check which types we have and set spring_constant accordingly

    CellPtr p_cell_A = rCellPopulation.GetCellUsingLocationIndex(nodeAGlobalIndex);
    CellPtr p_cell_B = rCellPopulation.GetCellUsingLocationIndex(nodeBGlobalIndex);

    // First, determine what we've got
    bool membraneA = p_cell_A->GetCellProliferativeType()->IsType<MembraneCellProliferativeType>();
    bool membraneB = p_cell_B->GetCellProliferativeType()->IsType<MembraneCellProliferativeType>();

    bool stromalA = p_cell_A->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>();
    bool stromalB = p_cell_B->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>();

    bool epiA = ( p_cell_A->GetCellProliferativeType()->IsType<TransitCellProliferativeType>() || p_cell_A->GetCellProliferativeType()->IsType<StemCellProliferativeType>() );
    bool epiB = ( p_cell_B->GetCellProliferativeType()->IsType<TransitCellProliferativeType>() || p_cell_B->GetCellProliferativeType()->IsType<StemCellProliferativeType>() );


    double rest_length_final = 1.0;
    //rest_length_final = static_cast<MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>*>(&rCellPopulation)->GetRestLength(nodeAGlobalIndex, nodeBGlobalIndex);

    double spring_constant = 0.0;

    // Determine rest lengths and spring stiffnesses
    if (membraneA)
    {
        if (membraneB)
        {
            if (distance_between_nodes >= mMembraneCutOffRadius)
            {
                return zero_vector<double>(SPACE_DIM);
            }
            rest_length_final = mMembraneRestLength;
            spring_constant = mMembraneSpringStiffness;
        }
        if (stromalB)
        {   
            if (distance_between_nodes >= mMembraneStromalCutOffRadius)
            {
                return zero_vector<double>(SPACE_DIM);
            }
            rest_length_final = mMembraneStromalRestLength;
            spring_constant = mMembraneStromalSpringStiffness;
        }
        if (epiB)
        {
            if (distance_between_nodes >= mEpithelialMembraneCutOffRadius)
            {
                return zero_vector<double>(SPACE_DIM);
            }
            rest_length_final = mEpithelialMembraneRestLength;
            spring_constant = mEpithelialMembraneSpringStiffness;
        }
    }

    if (stromalA)
    {
        if (membraneB)
        {
            if (distance_between_nodes >= mMembraneStromalCutOffRadius)
            {
                return zero_vector<double>(SPACE_DIM);
            }
            rest_length_final = mMembraneStromalRestLength;
            spring_constant = mMembraneStromalSpringStiffness;
        }
        if (stromalB)
        {
            if (distance_between_nodes >= mStromalCutOffRadius)
            {
                return zero_vector<double>(SPACE_DIM);
            }
            rest_length_final = mStromalRestLength;
            spring_constant = mStromalSpringStiffness;
        }
        if (epiB)
        {
            if (distance_between_nodes >= mStromalEpithelialCutOffRadius)
            {
                return zero_vector<double>(SPACE_DIM);
            }
            rest_length_final = mStromalEpithelialRestLength;
            spring_constant = mStromalEpithelialSpringStiffness;
        }
    }

    if (epiA)
    {
        if (membraneB)
        {
            if (distance_between_nodes >= mEpithelialMembraneCutOffRadius)
            {
                return zero_vector<double>(SPACE_DIM);
            }
            rest_length_final = mEpithelialMembraneRestLength;
            spring_constant = mEpithelialMembraneSpringStiffness;
        }
        if (stromalB)
        {
            rest_length_final = mStromalEpithelialRestLength;
            spring_constant = mStromalEpithelialSpringStiffness;
        }
        if (distance_between_nodes >= mStromalEpithelialCutOffRadius)
            {
                return zero_vector<double>(SPACE_DIM);
            }
        if (epiB)
        {
            if (distance_between_nodes >= mEpithelialCutOffRadius)
            {
                return zero_vector<double>(SPACE_DIM);
            }
            rest_length_final = mEpithelialRestLength;
            spring_constant = mEpithelialSpringStiffness;
        }
    }

    assert(spring_constant > 0);
    double rest_length = rest_length_final;

    double ageA = p_cell_A->GetAge();
    double ageB = p_cell_B->GetAge();

    assert(!std::isnan(ageA));
    assert(!std::isnan(ageB));

    /*
     * If the cells are both newly divided, then the rest length of the spring
     * connecting them grows linearly with time, until 1 hour after division.
     */
    if (ageA < mMeinekeSpringGrowthDuration && ageB < mMeinekeSpringGrowthDuration)
    {
        AbstractCentreBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>* p_static_cast_cell_population = static_cast<AbstractCentreBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>*>(&rCellPopulation);

        std::pair<CellPtr,CellPtr> cell_pair = p_static_cast_cell_population->CreateCellPair(p_cell_A, p_cell_B);

        if (p_static_cast_cell_population->IsMarkedSpring(cell_pair))
        {
            // Spring rest length increases from a Force value to the normal rest length over 1 hour
            double lambda = mMeinekeDivisionRestingSpringLength;
            rest_length = lambda + (rest_length_final - lambda) * ageA/mMeinekeSpringGrowthDuration;
        }
        if (ageA + SimulationTime::Instance()->GetTimeStep() >= mMeinekeSpringGrowthDuration)
        {
            // This spring is about to go out of scope
            p_static_cast_cell_population->UnmarkSpring(cell_pair);
        }
    }

    /*
     * For apoptosis, progressively reduce the radius of the cell
     */
    double a_rest_length = rest_length*0.5;
    double b_rest_length = a_rest_length;

    /*
     * If either of the cells has begun apoptosis, then the length of the spring
     * connecting them decreases linearly with time.
     */
    if (p_cell_A->HasApoptosisBegun())
    {
        double time_until_death_a = p_cell_A->GetTimeUntilDeath();
        a_rest_length = a_rest_length * time_until_death_a / p_cell_A->GetApoptosisTime();
    }
    if (p_cell_B->HasApoptosisBegun())
    {
        double time_until_death_b = p_cell_B->GetTimeUntilDeath();
        b_rest_length = b_rest_length * time_until_death_b / p_cell_B->GetApoptosisTime();
    }

    rest_length = a_rest_length + b_rest_length;


    double overlap = distance_between_nodes - rest_length;
    bool is_closer_than_rest_length = (overlap <= 0);

    if (is_closer_than_rest_length) //overlap is negative
    {
        //log(x+1) is undefined for x<=-1
        assert(overlap > -rest_length_final);
        c_vector<double, 2> temp = spring_constant * unitForceDirection * rest_length_final* log(1.0 + overlap/rest_length_final);
        return temp;
    }
    else
    {
        double alpha = 5.0;
        c_vector<double, 2> temp = spring_constant * unitForceDirection * overlap * exp(-alpha * overlap/rest_length_final);
        return temp;
    }

}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void LinearSpringForceMembraneCellNodeBased<ELEMENT_DIM,SPACE_DIM>::SetEpithelialSpringStiffness(double epithelialSpringStiffness)
{
    assert(epithelialSpringStiffness> 0.0);
    mEpithelialSpringStiffness = epithelialSpringStiffness;
}
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void LinearSpringForceMembraneCellNodeBased<ELEMENT_DIM,SPACE_DIM>::SetMembraneSpringStiffness(double membraneSpringStiffness)
{
    assert(membraneSpringStiffness > 0.0);
    mMembraneSpringStiffness = membraneSpringStiffness;
}
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void LinearSpringForceMembraneCellNodeBased<ELEMENT_DIM,SPACE_DIM>::SetStromalSpringStiffness(double stromalSpringStiffness)
{
    assert(stromalSpringStiffness > 0.0);
    mStromalSpringStiffness = stromalSpringStiffness;
}
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void LinearSpringForceMembraneCellNodeBased<ELEMENT_DIM,SPACE_DIM>::SetEpithelialMembraneSpringStiffness(double epithelialMembraneSpringStiffness)
{
    assert(epithelialMembraneSpringStiffness > 0.0);
    mEpithelialMembraneSpringStiffness = epithelialMembraneSpringStiffness;
}
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void LinearSpringForceMembraneCellNodeBased<ELEMENT_DIM,SPACE_DIM>::SetMembraneStromalSpringStiffness(double membraneStromalSpringStiffness)
{
    assert(membraneStromalSpringStiffness > 0.0);
    mMembraneStromalSpringStiffness = membraneStromalSpringStiffness;
}
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void LinearSpringForceMembraneCellNodeBased<ELEMENT_DIM,SPACE_DIM>::SetStromalEpithelialSpringStiffness(double stromalEpithelialSpringStiffness)
{
    assert(stromalEpithelialSpringStiffness > 0.0);
    mStromalEpithelialSpringStiffness = stromalEpithelialSpringStiffness;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void LinearSpringForceMembraneCellNodeBased<ELEMENT_DIM,SPACE_DIM>::SetEpithelialRestLength(double epithelialRestLength)
{
    assert(epithelialRestLength> 0.0);
    mEpithelialRestLength = epithelialRestLength;
}
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void LinearSpringForceMembraneCellNodeBased<ELEMENT_DIM,SPACE_DIM>::SetMembraneRestLength(double membraneRestLength)
{
    assert(membraneRestLength > 0.0);
    mMembraneRestLength = membraneRestLength;
}
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void LinearSpringForceMembraneCellNodeBased<ELEMENT_DIM,SPACE_DIM>::SetStromalRestLength(double stromalRestLength)
{
    assert(stromalRestLength > 0.0);
    mStromalRestLength = stromalRestLength;
}
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void LinearSpringForceMembraneCellNodeBased<ELEMENT_DIM,SPACE_DIM>::SetEpithelialMembraneRestLength(double epithelialMembraneRestLength)
{
    assert(epithelialMembraneRestLength > 0.0);
    mEpithelialMembraneRestLength = epithelialMembraneRestLength;
}
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void LinearSpringForceMembraneCellNodeBased<ELEMENT_DIM,SPACE_DIM>::SetMembraneStromalRestLength(double membraneStromalRestLength)
{
    assert(membraneStromalRestLength > 0.0);
    mMembraneStromalRestLength = membraneStromalRestLength;
}
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void LinearSpringForceMembraneCellNodeBased<ELEMENT_DIM,SPACE_DIM>::SetStromalEpithelialRestLength(double stromalEpithelialRestLength)
{
    assert(stromalEpithelialRestLength > 0.0);
    mStromalEpithelialRestLength = stromalEpithelialRestLength;
}



template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void LinearSpringForceMembraneCellNodeBased<ELEMENT_DIM,SPACE_DIM>::SetEpithelialCutOffRadius(double epithelialCutOffRadius)
{
    assert(epithelialCutOffRadius> 0.0);
    mEpithelialCutOffRadius = epithelialCutOffRadius;
}
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void LinearSpringForceMembraneCellNodeBased<ELEMENT_DIM,SPACE_DIM>::SetMembraneCutOffRadius(double membraneCutOffRadius)
{
    assert(membraneCutOffRadius > 0.0);
    mMembraneCutOffRadius = membraneCutOffRadius;
}
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void LinearSpringForceMembraneCellNodeBased<ELEMENT_DIM,SPACE_DIM>::SetStromalCutOffRadius(double stromalCutOffRadius)
{
    assert(stromalCutOffRadius > 0.0);
    mStromalCutOffRadius = stromalCutOffRadius;
}
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void LinearSpringForceMembraneCellNodeBased<ELEMENT_DIM,SPACE_DIM>::SetEpithelialMembraneCutOffRadius(double epithelialMembraneCutOffRadius)
{
    assert(epithelialMembraneCutOffRadius > 0.0);
    mEpithelialMembraneCutOffRadius = epithelialMembraneCutOffRadius;
}
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void LinearSpringForceMembraneCellNodeBased<ELEMENT_DIM,SPACE_DIM>::SetMembraneStromalCutOffRadius(double membraneStromalCutOffRadius)
{
    assert(membraneStromalCutOffRadius > 0.0);
    mMembraneStromalCutOffRadius = membraneStromalCutOffRadius;
}
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void LinearSpringForceMembraneCellNodeBased<ELEMENT_DIM,SPACE_DIM>::SetStromalEpithelialCutOffRadius(double stromalEpithelialCutOffRadius)
{
    assert(stromalEpithelialCutOffRadius > 0.0);
    mStromalEpithelialCutOffRadius = stromalEpithelialCutOffRadius;
}




template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void LinearSpringForceMembraneCellNodeBased<ELEMENT_DIM,SPACE_DIM>::SetMeinekeDivisionRestingSpringLength(double divisionRestingSpringLength)
{
    assert(divisionRestingSpringLength <= 1.0);
    assert(divisionRestingSpringLength >= 0.0);

    mMeinekeDivisionRestingSpringLength = divisionRestingSpringLength;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void LinearSpringForceMembraneCellNodeBased<ELEMENT_DIM,SPACE_DIM>::SetMeinekeSpringGrowthDuration(double springGrowthDuration)
{
    assert(springGrowthDuration >= 0.0);

    mMeinekeSpringGrowthDuration = springGrowthDuration;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void LinearSpringForceMembraneCellNodeBased<ELEMENT_DIM,SPACE_DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<EpithelialSpringStiffness>" << mEpithelialSpringStiffness << "</EpithelialSpringStiffness>\n";
    *rParamsFile << "\t\t\t<MembraneSpringStiffness>" << mMembraneSpringStiffness << "</MembraneSpringStiffness>\n";
    *rParamsFile << "\t\t\t<StromalSpringStiffness>" << mStromalSpringStiffness << "</StromalSpringStiffness>\n";
    *rParamsFile << "\t\t\t<EpithelialMembraneSpringStiffness>" << mEpithelialMembraneSpringStiffness << "</EpithelialMembraneSpringStiffness>\n";
    *rParamsFile << "\t\t\t<MembranetromalSpringStiffness>" << mMembraneStromalSpringStiffness << "</MembranetromalSpringStiffness>\n";
    *rParamsFile << "\t\t\t<StromalEpithelialSpringStiffness>" << mStromalEpithelialSpringStiffness << "</StromalEpithelialSpringStiffness>\n";

    *rParamsFile << "\t\t\t<EpithelialRestLength>" << mEpithelialRestLength << "</EpithelialRestLength>\n";
    *rParamsFile << "\t\t\t<MembraneRestLength>" << mMembraneRestLength << "</MembraneRestLength>\n";
    *rParamsFile << "\t\t\t<StromalRestLength>" << mStromalRestLength << "</StromalRestLength>\n";
    *rParamsFile << "\t\t\t<EpithelialMembraneRestLength>" << mEpithelialMembraneRestLength << "</EpithelialMembraneRestLength>\n";
    *rParamsFile << "\t\t\t<MembranetromalRestLength>" << mMembraneStromalRestLength << "</MembranetromalRestLength>\n";
    *rParamsFile << "\t\t\t<StromalEpithelialRestLength>" << mStromalEpithelialRestLength << "</StromalEpithelialRestLength>\n";

    *rParamsFile << "\t\t\t<EpithelialCutOffRadius>" << mEpithelialCutOffRadius << "</EpithelialCutOffRadius>\n";
    *rParamsFile << "\t\t\t<MembraneCutOffRadius>" << mMembraneCutOffRadius << "</MembraneCutOffRadius>\n";
    *rParamsFile << "\t\t\t<StromalCutOffRadius>" << mStromalCutOffRadius << "</StromalCutOffRadius>\n";
    *rParamsFile << "\t\t\t<EpithelialMembraneCutOffRadius>" << mEpithelialMembraneCutOffRadius << "</EpithelialMembraneCutOffRadius>\n";
    *rParamsFile << "\t\t\t<MembranetromalCutOffRadius>" << mMembraneStromalCutOffRadius << "</MembranetromalCutOffRadius>\n";
    *rParamsFile << "\t\t\t<StromalEpithelialCutOffRadius>" << mStromalEpithelialCutOffRadius << "</StromalEpithelialCutOffRadius>\n";

    *rParamsFile << "\t\t\t<MeinekeDivisionRestingSpringLength>" << mMeinekeDivisionRestingSpringLength << "</MeinekeDivisionRestingSpringLength>\n";
    *rParamsFile << "\t\t\t<MeinekeSpringGrowthDuration>" << mMeinekeSpringGrowthDuration << "</MeinekeSpringGrowthDuration>\n";

    // Call method on direct parent class
    AbstractTwoBodyInteractionForce<ELEMENT_DIM,SPACE_DIM>::OutputForceParameters(rParamsFile);
    // Call method on direct parent class
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class LinearSpringForceMembraneCellNodeBased<1,1>;
template class LinearSpringForceMembraneCellNodeBased<1,2>;
template class LinearSpringForceMembraneCellNodeBased<2,2>;
template class LinearSpringForceMembraneCellNodeBased<1,3>;
template class LinearSpringForceMembraneCellNodeBased<2,3>;
template class LinearSpringForceMembraneCellNodeBased<3,3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(LinearSpringForceMembraneCellNodeBased)