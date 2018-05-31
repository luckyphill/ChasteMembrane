// Standard includes for tests
#include <cxxtest/TestSuite.h>
#include "AbstractCellBasedTestSuite.hpp"
#include "CheckpointArchiveTypes.hpp"
#include "SmartPointers.hpp"
#include "CommandLineArguments.hpp"

// Simulators and output
#include "OffLatticeSimulation.hpp" //Simulates the evolution of the population

#include "LinearSpringForcePhaseBased.hpp"
#include "LinearSpringForceMembraneCellNodeBased.hpp"

#include "HoneycombMeshGenerator.hpp" //Generates mesh
#include "NodesOnlyMesh.hpp"
#include "CellsGenerator.hpp"

#include "NodeBasedCellPopulation.hpp"

// Proliferative types
#include "MembraneCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "StemCellProliferativeType.hpp"

//Cell cycle models
#include "NoCellCycleModel.hpp"
#include "UniformCellCycleModel.hpp"
#include "GrowingContactInhibitionPhaseBasedCCM.hpp"

#include "WildTypeCellMutationState.hpp"

// Boundary conditions
#include "BoundaryCellProperty.hpp"
#include "CryptBoundaryCondition.hpp"

#include "SimpleSloughingCellKiller.hpp"
#include "AnoikisCellKiller.hpp"

//Division Rules
#include "StickToMembraneDivisionRule.hpp"

// Modifiers
#include "VolumeTrackingModifier.hpp"

// Misc
#include "FakePetscSetup.hpp"
#include "Debug.hpp"

class TestGrowingCellDivision : public AbstractCellBasedTestSuite
{
	public:
	void xTestGrowingDivision() throw(Exception)
	{
		double dt = 0.01;
		double end_time = 20;
		double sampling_multiple = 10;

		double epithelialStiffness = 2.0; 			// 1.5
		double epithelialPreferredRadius = 0.7;			// 1.0
		double epithelialInteractionRadius = 1.5 * epithelialPreferredRadius; // Epithelial covers stem and transit
		double epithelialNewlyDividedRadius = 0.3;

		double stromalStiffness = 2.0; 			// 1.5
		double stromalPreferredRadius = 0.7;			// 1.0
		double stromalInteractionRadius = 1.5 * epithelialPreferredRadius;

		double stromalEpithelialStiffness = 1.0;


		std::vector<Node<2>*> nodes;

		unsigned cells_up = 5;
		unsigned cells_across = 5;
		unsigned node_counter = 0;

		for (unsigned i = 0; i< cells_across; i++)
		{
			for (unsigned j = 0; j< cells_up; j++)
			{
				double x = 0;
				double y = 0;
				if (j == 2* unsigned(j/2))
				{
					x= i;
				} else 
				{
					// stagger for hex mesh
					x = i +0.5;
				}
				y = j * (sqrt(3.0)/2);
				nodes.push_back(new Node<2>(node_counter,  false,  x, y));
				node_counter++;
			}
		}

		NodesOnlyMesh<2> mesh;
		mesh.ConstructNodesWithoutMesh(nodes, 3.0);

		std::vector<CellPtr> cells;

		MAKE_PTR(MembraneCellProliferativeType, p_membrane_type);
		MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
		MAKE_PTR(StemCellProliferativeType, p_stem_type);
		MAKE_PTR(TransitCellProliferativeType, p_trans_type);

		MAKE_PTR(WildTypeCellMutationState, p_state);

		for (unsigned i = 0; i < nodes.size(); i++)
		{
			if (i==12)
			{
				// Set the middle cell to be proliferating
				GrowingContactInhibitionPhaseBasedCCM* p_cycle_model = new GrowingContactInhibitionPhaseBasedCCM();
				p_cycle_model->SetNewlyDividedRadius(epithelialNewlyDividedRadius);
				p_cycle_model->SetPreferredRadius(epithelialPreferredRadius);
				p_cycle_model->SetInteractionRadius(epithelialInteractionRadius);

				CellPtr p_cell(new Cell(p_state, p_cycle_model));
				p_cell->SetCellProliferativeType(p_trans_type);

				p_cell->InitialiseCellCycleModel();

				cells.push_back(p_cell);
			}
			else
			{
				NoCellCycleModel* p_cycle_model = new NoCellCycleModel();
			
				CellPtr p_cell(new Cell(p_state, p_cycle_model));
				p_cell->SetCellProliferativeType(p_diff_type);
	
				p_cell->InitialiseCellCycleModel();
	
				cells.push_back(p_cell);
			}
		}

        NodeBasedCellPopulation<2> cell_population(mesh, cells);

        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestGrowingDivision");
        simulator.SetSamplingTimestepMultiple(sampling_multiple);
        simulator.SetEndTime(end_time);
        simulator.SetDt(dt);

        MAKE_PTR(LinearSpringForcePhaseBased<2>, p_force);
        //MAKE_PTR(LinearSpringForceMembraneCellNodeBased<2>, p_force);
		p_force->SetEpithelialSpringStiffness(epithelialStiffness);
		p_force->SetEpithelialPreferredRadius(epithelialPreferredRadius);
		p_force->SetEpithelialInteractionRadius(epithelialInteractionRadius);

		p_force->SetStromalSpringStiffness(stromalStiffness);
		p_force->SetStromalPreferredRadius(stromalPreferredRadius);
		p_force->SetStromalInteractionRadius(stromalInteractionRadius);

		p_force->SetStromalEpithelialSpringStiffness(stromalEpithelialStiffness);

		p_force->SetMeinekeSpringGrowthDuration(1);
		p_force->SetMeinekeDivisionRestingSpringLength(0.1);

		p_force->SetDebugMode(false);

        simulator.AddForce(p_force);

        MAKE_PTR(VolumeTrackingModifier<2>, p_mod);
		simulator.AddSimulationModifier(p_mod);

        simulator.Solve();
	};

	void TestGrowingDivisionWithMembrane() throw(Exception)
	{

		//TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-wh"));
        double wh = 11;// CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-wh");


		bool debugging = false;

        double epithelialStiffness = 2.0;
		
        double epithelialMembraneStiffness = 5.0;
		
		std::vector<Node<2>*> nodes;
		std::vector<unsigned> transit_nodes;
		std::vector<unsigned> membrane_nodes;
		std::vector<unsigned> location_indices;
		std::vector<std::vector<CellPtr>> membraneSections;

		unsigned node_counter = 0;

		double dt = 0.05;
		double end_time = 50;
		double sampling_multiple = 1;

		// Values that produce a working simulation in the comments
		double membraneStiffness = 5; 			// 5.0

		double epithelialPreferredRadius = 0.75;			// 1.0
		double membranePreferredRadius = 0.2;			// 0.2

		double epithelialInteractionRadius = 1.5 * epithelialPreferredRadius; // Epithelial covers stem and transit
		double membraneInteractionRadius = 7.0 * membranePreferredRadius;
		double epithelialNewlyDividedRadius = 0.5;
		double maxInteractionRadius = 3.0;

		double minCellCycleDuration = 5;

		double mDuration = 1;
		double g1Duration = 8;
		double sDuration = 7.5;
		double g2Duration = 1.5;

		double equilibriumVolume = 0.7;
		double volumeFraction = 0.88;


		double springGrowthDuration = 1.0;


		double membrane_spacing = 0.2;
		double epithelial_spacing = 1.5 * epithelialPreferredRadius;
		double wall_height = wh;
		double left_side = 0;
		double wall_top = wall_height;
		double wall_bottom = 0;

		// Drawing the membrane
		for (double y = wall_bottom; y <= wall_top; y+=membrane_spacing)
		{
			double x = left_side;
			nodes.push_back(new Node<2>(node_counter,  false,  x, y));
			membrane_nodes.push_back(node_counter);
			location_indices.push_back(node_counter);
			node_counter++;
		}

		// Drawing the epithelium
		// The transit amplifying cells
		for (double y = wall_bottom; y <= wall_top; y+= epithelial_spacing)
		{
			double x = 0.88; // Note this value is determined from observing simulations//left_side + epithelialPreferredRadius;
			nodes.push_back(new Node<2>(node_counter,  false,  x, y));
			transit_nodes.push_back(node_counter);
			location_indices.push_back(node_counter);
			node_counter++;
		}

		NodesOnlyMesh<2> mesh;
		mesh.ConstructNodesWithoutMesh(nodes, maxInteractionRadius);

		std::vector<CellPtr> cells;
		std::vector<CellPtr> membrane_cells;

		MAKE_PTR(MembraneCellProliferativeType, p_membrane_type);
		MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
		MAKE_PTR(TransitCellProliferativeType, p_trans_type);
		MAKE_PTR(WildTypeCellMutationState, p_state);
		MAKE_PTR(BoundaryCellProperty, p_boundary);

		//Initialise membrane nodes
		for (unsigned i = 0; i < membrane_nodes.size(); i++)
		{
			NoCellCycleModel* p_cycle_model = new NoCellCycleModel();

			CellPtr p_cell(new Cell(p_state, p_cycle_model));
			p_cell->SetCellProliferativeType(p_membrane_type);
			
			p_cell->AddCellProperty(p_boundary);

			p_cell->InitialiseCellCycleModel();

			cells.push_back(p_cell);
			membrane_cells.push_back(p_cell);
		}

		// First node is fixed
		{
			GrowingContactInhibitionPhaseBasedCCM* p_cycle_model = new GrowingContactInhibitionPhaseBasedCCM();

			CellPtr p_cell(new Cell(p_state, p_cycle_model));
			p_cell->SetCellProliferativeType(p_diff_type);
			p_cell->AddCellProperty(p_boundary);

			p_cell->InitialiseCellCycleModel();

			cells.push_back(p_cell);
		}
		//Initialise trans nodes
		for (unsigned i = 1; i < transit_nodes.size(); i++)
		{
			GrowingContactInhibitionPhaseBasedCCM* p_cycle_model = new GrowingContactInhibitionPhaseBasedCCM();
			p_cycle_model->SetNewlyDividedRadius(epithelialNewlyDividedRadius);
			p_cycle_model->SetPreferredRadius(epithelialPreferredRadius);
			p_cycle_model->SetInteractionRadius(epithelialInteractionRadius);
			p_cycle_model->SetMDuration(mDuration);
			p_cycle_model->SetSDuration(sDuration);
			p_cycle_model->SetTransitCellG1Duration(g1Duration);
			p_cycle_model->SetG2Duration(g2Duration);
			p_cycle_model->SetEquilibriumVolume(equilibriumVolume);
			p_cycle_model->SetQuiescentVolumeFraction(volumeFraction);

			double birth_time = minCellCycleDuration * RandomNumberGenerator::Instance()->ranf(); //Randomly set birth time to stop pulsing behaviour
			p_cycle_model->SetBirthTime(-birth_time);
			//p_cycle_model->SetMinCellCycleDuration(minCellCycleDuration);

			CellPtr p_cell(new Cell(p_state, p_cycle_model));
			p_cell->SetCellProliferativeType(p_trans_type);

			p_cell->InitialiseCellCycleModel();

			cells.push_back(p_cell);
		}

		membraneSections.push_back(membrane_cells);

		NodeBasedCellPopulation<2> cell_population(mesh, cells, location_indices);

		{ //Division vector rules
			c_vector<double, 2> membraneAxis;
			membraneAxis(0) = 0;
			membraneAxis(1) = 1;

			MAKE_PTR(StickToMembraneDivisionRule<2>, pCentreBasedDivisionRule);
			pCentreBasedDivisionRule->SetMembraneAxis(membraneAxis);
			cell_population.SetCentreBasedDivisionRule(pCentreBasedDivisionRule);
		}

		OffLatticeSimulation<2> simulator(cell_population);

		simulator.SetOutputDirectory("TestGrowingDivisionWithMembrane");

		simulator.SetEndTime(end_time);
		simulator.SetDt(dt);
		simulator.SetSamplingTimestepMultiple(sampling_multiple);

		MAKE_PTR(LinearSpringForcePhaseBased<2>, p_force);
		p_force->SetEpithelialSpringStiffness(epithelialStiffness);
		p_force->SetMembraneSpringStiffness(membraneStiffness);
		p_force->SetEpithelialMembraneSpringStiffness(epithelialMembraneStiffness);

		p_force->SetEpithelialPreferredRadius(epithelialPreferredRadius);
		p_force->SetMembranePreferredRadius(membranePreferredRadius);

		p_force->SetEpithelialInteractionRadius(epithelialInteractionRadius);
		p_force->SetMembraneInteractionRadius(membraneInteractionRadius);
		
		p_force->SetMeinekeSpringGrowthDuration(springGrowthDuration);
		p_force->SetMeinekeDivisionRestingSpringLength(0.1);

		p_force->SetDebugMode(debugging);

		simulator.AddForce(p_force);

		MAKE_PTR_ARGS(AnoikisCellKiller, p_anoikis_killer, (&cell_population));
		simulator.AddCellKiller(p_anoikis_killer);

		MAKE_PTR_ARGS(CryptBoundaryCondition, p_bc, (&cell_population));
		simulator.AddCellPopulationBoundaryCondition(p_bc);

		MAKE_PTR_ARGS(SimpleSloughingCellKiller, p_sloughing_killer, (&cell_population));
		p_sloughing_killer->SetCryptTop(wall_top);
		simulator.AddCellKiller(p_sloughing_killer);

		MAKE_PTR(VolumeTrackingModifier<2>, p_mod);
		simulator.AddSimulationModifier(p_mod);

		simulator.Solve();

	};
};