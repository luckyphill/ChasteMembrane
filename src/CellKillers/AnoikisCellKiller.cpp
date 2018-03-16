/*
 * LAST MODIFIED: 02/10/2015
 * Anoikis cell killer created for epithelial layer model. Removes any epithelial cells that have detached from
 * the non-epithelial region and entered the lumen.
 *
 * Created on: Dec 21 2014
 * Last modified:
 * 			Author: Axel Almet
 */

/*
 * MODIFIED BY PHILLIP BROWN 1/11/2017
 * Added in a check for an anoikis resistant mutation, which stops a cell from dying when it
 * detaches from the non epithelial region i.e. the basement layer of cells
 * Also modified to check if cell is in contact with MembraneCellProliferativeType
 */

#include "AnoikisCellKiller.hpp"
#include "AnoikisCellTagged.hpp"
#include "AbstractCellKiller.hpp"
#include "AbstractCellProperty.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "PanethCellMutationState.hpp"
#include "TransitCellAnoikisResistantMutationState.hpp"
#include "MembraneCellProliferativeType.hpp"

AnoikisCellKiller::AnoikisCellKiller(AbstractCellPopulation<2>* pCellPopulation)
    : AbstractCellKiller<2>(pCellPopulation),
    mCellsRemovedByAnoikis(0),
    mCutOffRadius(1.5),
    mSlowDeath(false)
{
    // Sets up output file
//	OutputFileHandler output_file_handler(mOutputDirectory + "AnoikisData/", false);
//	mAnoikisOutputFile = output_file_handler.OpenOutputFile("results.anoikis");
}

AnoikisCellKiller::~AnoikisCellKiller()
{
//    mAnoikisOutputFile->close();
}

//Method to get mCutOffRadius
double AnoikisCellKiller::GetCutOffRadius()
{
	return mCutOffRadius;
}

//Method to set mCutOffRadius
void AnoikisCellKiller::SetCutOffRadius(double cutOffRadius)
{
	mCutOffRadius = cutOffRadius;
}

std::set<unsigned> AnoikisCellKiller::GetNeighbouringNodeIndices(unsigned nodeIndex)
{
	// Create a set of neighbouring node indices
	std::set<unsigned> neighbouring_node_indices;

	if (dynamic_cast<NodeBasedCellPopulation<2>*>(this->mpCellPopulation))
	{
		// Need access to the mesh but can't get to it because the cell killer only owns a
		// pointer to an AbstractCellPopulation
		NodeBasedCellPopulation<2>* p_tissue = static_cast<NodeBasedCellPopulation<2>*> (this->mpCellPopulation);

		//Update cell population
		p_tissue->Update();

		double radius = GetCutOffRadius();

		neighbouring_node_indices = p_tissue->GetNodesWithinNeighbourhoodRadius(nodeIndex, radius);
	}

    return neighbouring_node_indices;
}

/** Method to determine if an epithelial cell has lost all contacts with the gel cells below
 * TRUE if cell has popped up
 * FALSE if cell remains in the monolayer
 */
bool AnoikisCellKiller::HasCellPoppedUp(unsigned nodeIndex)
{
	bool has_cell_popped_up = false;	// Initialising
	std::set<unsigned> neighbours = GetNeighbouringNodeIndices(nodeIndex);

	if (dynamic_cast<NodeBasedCellPopulation<2>*>(this->mpCellPopulation))
	{
		NodeBasedCellPopulation<2>* p_tissue = static_cast<NodeBasedCellPopulation<2>*> (this->mpCellPopulation);

		unsigned membrane_neighbours = 0;

		// Iterate over the neighbouring cells to check the number of differentiated cell neighbours

		for(std::set<unsigned>::iterator neighbour_iter=neighbours.begin();
				neighbour_iter != neighbours.end();
				++neighbour_iter)
		{
			if (p_tissue->GetCellUsingLocationIndex(*neighbour_iter)->GetCellProliferativeType()->IsType<MembraneCellProliferativeType>() )
			{
				membrane_neighbours += 1;
			}
		}

		if(membrane_neighbours < 1)
		{
			has_cell_popped_up = true;
		}
	}

	return has_cell_popped_up;
}

void AnoikisCellKiller::PopulateAnoikisList()
{
	// Loop through, check if popped up and if so, store the cell pointer and the time

	if (dynamic_cast<NodeBasedCellPopulation<2>*>(this->mpCellPopulation))
    {
    	NodeBasedCellPopulation<2>* p_tissue = static_cast<NodeBasedCellPopulation<2>*> (this->mpCellPopulation);

    	for (AbstractCellPopulation<2>::Iterator cell_iter = p_tissue->Begin();
    			cell_iter != p_tissue->End();
    			++cell_iter)
    	{
    		unsigned node_index = p_tissue->GetNodeCorrespondingToCell(*cell_iter)->GetIndex();
    		if ()
    		if (HasCellPoppedUp(node_index))
    		{
    			//Add to list

    			std::pair<CellPtr, double> cell_data;
    			cell_data[0] = cell_iter;
    			cell_data[1] = SimulationTime::Instance()->GetTime();
    			mCellsForDelayedAnoikis.push_back(cell_data);
    		}
    	}
    }

}

std::std::vector<CellPtr> AnoikisCellKiller::GetCellsReadyToDie()
{
	// Go through the anoikis list, if the lenght of time since it popped up is past a certain
	// threshold, then that cell is ready to be killed
}
/** A method to return a vector that indicates which cells should be killed by anoikis
 * and which by compression-driven apoptosis
 */
std::vector<c_vector<unsigned,2> > AnoikisCellKiller::RemoveByAnoikis()
{

if (dynamic_cast<NodeBasedCellPopulation<2>*>(this->mpCellPopulation))
    {
    	NodeBasedCellPopulation<2>* p_tissue = static_cast<NodeBasedCellPopulation<2>*> (this->mpCellPopulation);

    	c_vector<unsigned,2> individual_node_information;	// Will store the node index and whether to remove or not (1 or 0)

    	for (AbstractCellPopulation<2>::Iterator cell_iter = p_tissue->Begin();
    			cell_iter != p_tissue->End();
    			++cell_iter)
    	{
    		unsigned node_index = p_tissue->GetNodeCorrespondingToCell(*cell_iter)->GetIndex();

    		// Initialise
    		individual_node_information[0] = node_index;
    		individual_node_information[1] = 0;

    		// Examine each epithelial node to see if it should be removed by anoikis and then if it
    		// should be removed by compression-driven apoptosis
    		// Edit by Phillip Brown: Added a check for anoikis resistant mutation to prevent this kind of cell death
    		if (!cell_iter->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>()
    			&& !cell_iter->GetMutationState()->IsType<TransitCellAnoikisResistantMutationState>())
    		{
    			// Determining whether to remove this cell by anoikis

    			if(this->HasCellPoppedUp(node_index))
    			{
    				individual_node_information[1] = 1;
    			}
    		}

    		cells_to_remove.push_back(individual_node_information);
    	}
    }

	return cells_to_remove;
}


/*
 * Cell Killer that kills healthy cells that pop outwards and become detached from
 * the labelled tissue cells, i.e. removal by anoikis
 *
 * Also will remove differentiated cells at the orifice if mSloughOrifice is true
 */
void AnoikisCellKiller::CheckAndLabelCellsForApoptosisOrDeath()
{
	if (dynamic_cast<MeshBasedCellPopulation<2>*>(this->mpCellPopulation))
	{
		MeshBasedCellPopulation<2>* p_tissue = static_cast<MeshBasedCellPopulation<2>*> (this->mpCellPopulation);
		//    assert(p_tissue->GetVoronoiTessellation()!=NULL);	// This fails during archiving of a simulation as Voronoi stuff not archived yet

		// Get the information at this timestep for each node index that says whether to remove by anoikis or random apoptosis
		std::vector<c_vector<unsigned,2> > cells_to_remove = this->RemoveByAnoikis();

		// Keep a record of how many cells have been removed at this timestep
		this->SetNumberCellsRemoved(cells_to_remove);
		this->SetLocationsOfCellsRemovedByAnoikis(cells_to_remove);

		// Need to avoid trying to kill any cells twice (i.e. both by anoikis or sloughing)
		// Loop over these vectors individually and kill any cells that they tell you to

		for (unsigned i=0; i<cells_to_remove.size(); i++)
		{
			if (cells_to_remove[i][1] == 1)
			{
				
				// Get cell associated to this node
				CellPtr p_cell = p_tissue->GetCellUsingLocationIndex(cells_to_remove[i][0]);
				if (mSlowDeath)
				{
					if (!p_cell->HasApoptosisBegun())
					{
						p_cell->StartApoptosis();
					}
				} else {
					p_cell->Kill();
				}
			}
		}
	}
	else if (dynamic_cast<NodeBasedCellPopulation<2>*>(this->mpCellPopulation))
	{
		NodeBasedCellPopulation<2>* p_tissue = static_cast<NodeBasedCellPopulation<2>*> (this->mpCellPopulation);

		// Get the information at this timestep for each node index that says whether to remove by anoikis or random apoptosis
		std::vector<c_vector<unsigned,2> > cells_to_remove = this->RemoveByAnoikis();

		// Keep a record of how many cells have been removed at this timestep
		this->SetNumberCellsRemoved(cells_to_remove);
		this->SetLocationsOfCellsRemovedByAnoikis(cells_to_remove);

		// Need to avoid trying to kill any cells twice (i.e. both by anoikis or sloughing)
		// Loop over these vectors individually and kill any cells that they tell you to

		for (unsigned i=0; i<cells_to_remove.size(); i++)
		{
			if (cells_to_remove[i][1] == 1)
			{
				// Get cell associated to this node
				CellPtr p_cell = p_tissue->GetCellUsingLocationIndex(cells_to_remove[i][0]);
				if (mSlowDeath)
				{
					if (!p_cell->HasApoptosisBegun())
					{
						p_cell->StartApoptosis();
					}
				} else {
					p_cell->Kill();
				}
			}
		}
	}
}

void AnoikisCellKiller::SetNumberCellsRemoved(std::vector<c_vector<unsigned,2> > cellsRemoved)
{
	unsigned num_removed_by_anoikis = 0;

    for (unsigned i=0; i<cellsRemoved.size(); i++)
    {
    	if(cellsRemoved[i][1]==1)
    	{
    		num_removed_by_anoikis+=1;
    	}
    }

    mCellsRemovedByAnoikis += num_removed_by_anoikis;
}

unsigned AnoikisCellKiller::GetNumberCellsRemoved()
{
	return mCellsRemovedByAnoikis;
}

void AnoikisCellKiller::SetLocationsOfCellsRemovedByAnoikis(std::vector<c_vector<unsigned,2> > cellsRemoved)
{
	if (dynamic_cast<MeshBasedCellPopulation<2>*>(this->mpCellPopulation))
	{
		MeshBasedCellPopulation<2>* p_tissue = static_cast<MeshBasedCellPopulation<2>*> (this->mpCellPopulation);
		double x_location, y_location;
		c_vector<double, 3> time_and_location;

		// Need to use the node indices to store the locations of where cells are removed
		for (unsigned i=0; i<cellsRemoved.size(); i++)
		{
			if (cellsRemoved[i][1] == 1)		// This cell has been removed by anoikis
			{
				time_and_location[0] = SimulationTime::Instance()->GetTime();

				unsigned node_index = cellsRemoved[i][0];

				CellPtr p_cell = p_tissue->GetCellUsingLocationIndex(node_index);
				x_location = this->mpCellPopulation->GetLocationOfCellCentre(p_cell)[0];
				y_location = this->mpCellPopulation->GetLocationOfCellCentre(p_cell)[1];

				time_and_location[1] = x_location;
				time_and_location[2] = y_location;

				mLocationsOfAnoikisCells.push_back(time_and_location);
			}
		}
	}
	else if (dynamic_cast<NodeBasedCellPopulation<2>*>(this->mpCellPopulation))
	{
		NodeBasedCellPopulation<2>* p_tissue = static_cast<NodeBasedCellPopulation<2>*> (this->mpCellPopulation);
		double x_location, y_location;
		c_vector<double, 3> time_and_location;

		// Need to use the node indices to store the locations of where cells are removed
		for (unsigned i=0; i<cellsRemoved.size(); i++)
		{
			if (cellsRemoved[i][1] == 1)		// This cell has been removed by anoikis
			{
				time_and_location[0] = SimulationTime::Instance()->GetTime();

				unsigned node_index = cellsRemoved[i][0];

				CellPtr p_cell = p_tissue->GetCellUsingLocationIndex(node_index);
				x_location = this->mpCellPopulation->GetLocationOfCellCentre(p_cell)[0];
				y_location = this->mpCellPopulation->GetLocationOfCellCentre(p_cell)[1];

				time_and_location[1] = x_location;
				time_and_location[2] = y_location;

				mLocationsOfAnoikisCells.push_back(time_and_location);
			}
		}
	}
}

std::vector<c_vector<double,3> > AnoikisCellKiller::GetLocationsOfCellsRemovedByAnoikis()
{
	return mLocationsOfAnoikisCells;
}

void AnoikisCellKiller::SetSlowDeath(bool slowDeath)
{
	mSlowDeath = slowDeath;
}

void AnoikisCellKiller::OutputCellKillerParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<CellsRemovedByAnoikis>" << mCellsRemovedByAnoikis << "</CellsRemovedByAnoikis> \n";
    *rParamsFile << "\t\t\t<CutOffRadius>" << mCutOffRadius << "</CutOffRadius> \n";

    // Call direct parent class
    AbstractCellKiller<2>::OutputCellKillerParameters(rParamsFile);
}




#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(AnoikisCellKiller)
