// Standard includes for tests
#include <cxxtest/TestSuite.h>
#include "AbstractCellBasedTestSuite.hpp"
#include "CheckpointArchiveTypes.hpp"
#include "SmartPointers.hpp"
#include "CommandLineArguments.hpp"

// Simulators and output
#include "OffLatticeSimulation.hpp" //Simulates the evolution of the population
#include "OffLatticeSimulationTearOffStoppingEvent.hpp"

#include "LinearSpringForcePhaseBased.hpp"
#include "LinearSpringForceMembraneCellNodeBased.hpp"
#include "PushForce.hpp"

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

#include "PushForceModifier.hpp"
#include "BasicNonLinearSpringForce.hpp"


// Misc
#include "FakePetscSetup.hpp"
#include "Debug.hpp"

class TestMembraneCellWall : public AbstractCellBasedTestSuite
{

public:
	void xTestSingleCellOnMembrane() throw(Exception)
	{
		// This test can be used to observe how a cell interacts with the membrane layer
		// You can add a force to the cell to see how it moves along the wall
		// You can change how far it starts from the wall to see how it is pulled in


		TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-x"));
        double x_distance = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-x");

        TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-y"));
        double y_distance = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-y");

        TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-yf"));
        double y_force = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-yf");

        TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-xf"));
        double x_force = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-xf");

        bool multiple_cells = false;
        unsigned n = 0;
        if(CommandLineArguments::Instance()->OptionExists("-n"))
        {	
        	multiple_cells = true;
        	n = CommandLineArguments::Instance()->GetUnsignedCorrespondingToOption("-n");

        }

		bool debugging = true;

        double epithelialStiffness = 2.0;
		
        double epithelialMembraneStiffness = 5.0;
		
		std::vector<Node<2>*> nodes;
		std::vector<unsigned> transit_nodes;
		std::vector<unsigned> membrane_nodes;
		std::vector<unsigned> location_indices;
		std::vector<std::vector<CellPtr>> membraneSections;

		unsigned node_counter = 0;

		double dt = 0.01;
		double end_time = 5;
		double sampling_multiple = 1;

		// Values that produce a working simulation in the comments
		double membraneStiffness = 5.0; 			// 5.0

		double epithelialPreferredRadius = 0.75;			// 1.0
		double membranePreferredRadius = 0.2;			// 0.2

		double epithelialInteractionRadius = 1.5 * epithelialPreferredRadius; // Epithelial covers stem and transit
		double membraneInteractionRadius = 10.0 * membranePreferredRadius;
		double epithelialNewlyDividedRadius = 0.5;
		double maxInteractionRadius = 4.0;


		double springGrowthDuration = 1.0;


		double membrane_spacing = 0.2;
		double epithelial_spacing = 1.5 * epithelialPreferredRadius;
		double wall_height = 10;
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

		// Placing a single cell on the wall
		
		double x = x_distance;
		double y = y_distance;
		Node<2>* single_node =  new Node<2>(node_counter,  false,  x, y);
		nodes.push_back(single_node);
		transit_nodes.push_back(node_counter);
		location_indices.push_back(node_counter);
		node_counter++;

		if(multiple_cells)
		{
			for(unsigned i=1; i<=n; i++)
			{
				x = x_distance;
				y = y_distance + 2 * i * epithelialPreferredRadius;
				Node<2>* single_node_2 =  new Node<2>(node_counter,  false,  x, y);
				nodes.push_back(single_node_2);
				transit_nodes.push_back(node_counter);
				location_indices.push_back(node_counter);
				node_counter++;
			}
			
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

		// Make the single cell

		NoCellCycleModel* p_cycle_model = new NoCellCycleModel();

		CellPtr p_cell(new Cell(p_state, p_cycle_model));
		p_cell->SetCellProliferativeType(p_diff_type);

		cells.push_back(p_cell);

		// Add any additional static cells
		if(multiple_cells)
		{
			for(unsigned i=1; i<=n; i++)
			{
				CellPtr p_cell_2(new Cell(p_state, p_cycle_model));
				p_cell_2->SetCellProliferativeType(p_diff_type);

				cells.push_back(p_cell_2);
			}
		}
		
		membraneSections.push_back(membrane_cells);

		NodeBasedCellPopulation<2> cell_population(mesh, cells, location_indices);

		OffLatticeSimulation<2> simulator(cell_population);

		simulator.SetOutputDirectory("TestSingleCellOnMembrane");

		simulator.SetEndTime(end_time);
		simulator.SetDt(dt);
		simulator.SetSamplingTimestepMultiple(sampling_multiple);

		MAKE_PTR(LinearSpringForceMembraneCellNodeBased<2>, p_force);
		p_force->SetEpithelialSpringStiffness(epithelialStiffness);
		p_force->SetStromalSpringStiffness(epithelialStiffness);
		p_force->SetMembraneSpringStiffness(membraneStiffness);
		
		p_force->SetEpithelialMembraneSpringStiffness(epithelialMembraneStiffness);
		p_force->SetMembraneStromalSpringStiffness(epithelialMembraneStiffness);
		p_force->SetStromalEpithelialSpringStiffness(epithelialStiffness);

		p_force->SetEpithelialPreferredRadius(epithelialPreferredRadius);
		p_force->SetStromalPreferredRadius(epithelialPreferredRadius);
		p_force->SetMembranePreferredRadius(membranePreferredRadius);

		p_force->SetEpithelialInteractionRadius(epithelialInteractionRadius);
		p_force->SetStromalInteractionRadius(epithelialInteractionRadius);
		p_force->SetMembraneInteractionRadius(membraneInteractionRadius);
		
		p_force->SetMeinekeSpringGrowthDuration(springGrowthDuration);
		p_force->SetMeinekeDivisionRestingSpringLength(0.1);

		p_force->SetDebugMode(debugging);

		simulator.AddForce(p_force);

		MAKE_PTR_ARGS(CryptBoundaryCondition, p_bc, (&cell_population));
		simulator.AddCellPopulationBoundaryCondition(p_bc);

		MAKE_PTR(PushForce, p_push_force);
		p_push_force->SetCell(p_cell);
		c_vector<double, 2> force;
		force[0] = x_force;
		force[1] = y_force;
		p_push_force->SetForce(force);
		simulator.AddForce(p_push_force);

		simulator.Solve();

	};


	void xTestSingleCellRestPosition() throw(Exception)
	{
		// This test is to make sure that the force calculations match
		// the hand calculations as provided by matlab


		TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-x"));
        double x_distance = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-x");

        TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-y"));
        double y_distance = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-y");

        TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-ms"));
        double membrane_spacing = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-ms"); // Distance between membrane nodes

        TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-ems"));
        double epithelialMembraneStiffness = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-ems"); // The only spring stiffness to worry about

        TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-mir"));
        double membraneInteractionRadius = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-mir"); // The furthest that a membrane node can detect the epithelial cell

        TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-mpr"));
        double membranePreferredRadius = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-mpr"); // The way that the natural spring length of the membrane-epithelial connection is controlled

		
		std::vector<Node<2>*> nodes;
		std::vector<unsigned> transit_nodes;
		std::vector<unsigned> membrane_nodes;
		std::vector<unsigned> location_indices;
		std::vector<std::vector<CellPtr>> membraneSections;

		unsigned node_counter = 0;

		double dt = 0.01;
		double end_time = 2;
		double sampling_multiple = 1;
		
		double maxInteractionRadius = 4.0;

		double wall_height = 5;
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

		// Placing a single cell on the wall
		
		double x = x_distance;
		double y = y_distance;
		Node<2>* single_node =  new Node<2>(node_counter,  false,  x, y);
		nodes.push_back(single_node);
		transit_nodes.push_back(node_counter);
		location_indices.push_back(node_counter);
		node_counter++;
	

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

		// Make the single cell

		NoCellCycleModel* p_cycle_model = new NoCellCycleModel();

		CellPtr p_cell(new Cell(p_state, p_cycle_model));
		p_cell->SetCellProliferativeType(p_diff_type);

		cells.push_back(p_cell);

		membraneSections.push_back(membrane_cells);

		NodeBasedCellPopulation<2> cell_population(mesh, cells, location_indices);

		OffLatticeSimulation<2> simulator(cell_population);

		simulator.SetOutputDirectory("TestSingleCellOnMembrane");

		simulator.SetEndTime(end_time);
		simulator.SetDt(dt);
		simulator.SetSamplingTimestepMultiple(sampling_multiple);

		MAKE_PTR(BasicNonLinearSpringForce<2>, p_force);


		p_force->SetEpithelialMembraneSpringStiffness(epithelialMembraneStiffness);

		p_force->SetEpithelialMembraneRestLength(membranePreferredRadius);

		p_force->SetEpithelialMembraneCutOffLength(membraneInteractionRadius);

		simulator.AddForce(p_force);

		MAKE_PTR_ARGS(CryptBoundaryCondition, p_bc, (&cell_population));
		simulator.AddCellPopulationBoundaryCondition(p_bc);

		simulator.Solve();

		Node<2>* p_node =  cell_population.GetNodeCorrespondingToCell(p_cell);
        c_vector<double, 2> location = p_node->rGetLocation();

		ofstream myfile;
		myfile.open ("rest_position.txt", ios::app);
		myfile << x_distance << ","<< y_distance << ","<< membrane_spacing << ","<< epithelialMembraneStiffness << ","<< membraneInteractionRadius << ","<< membranePreferredRadius << "|";
		myfile << location[0] << "," << location[1] << "\n";
		myfile.close();

	};

	void xTestSingleCellMigrationForce() throw(Exception)
	{
		// Determine the minimum force for any movement, and the minimum force for total movement
	};

	void TestSingleCellTearOffForce() throw(Exception)
	{
		// Determine the minimum force to tear a cell off the wall
		TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-x"));
        double x_distance = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-x");

        TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-y"));
        double y_distance = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-y");

        TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-xf"));
        double x_force = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-xf");

        TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-ms"));
        double membrane_spacing = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-ms"); // Distance between membrane nodes

        TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-ems"));
        double epithelialMembraneStiffness = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-ems"); // The only spring stiffness to worry about

        TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-mir"));
        double membraneInteractionRadius = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-mir"); // The furthest that a membrane node can detect the epithelial cell

        TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-mpr"));
        double membranePreferredRadius = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-mpr"); // The way that the natural spring length of the membrane-epithelial connection is controlled

		
		std::vector<Node<2>*> nodes;
		std::vector<unsigned> transit_nodes;
		std::vector<unsigned> membrane_nodes;
		std::vector<unsigned> location_indices;
		std::vector<std::vector<CellPtr>> membraneSections;

		unsigned node_counter = 0;

		double dt = 0.01;
		double end_time = 20;
		double sampling_multiple = 1;
		
		double maxInteractionRadius = 4.0;

		double wall_height = 5;
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

		// Placing a single cell on the wall
		
		double x = x_distance;
		double y = y_distance;
		Node<2>* single_node =  new Node<2>(node_counter,  false,  x, y);
		nodes.push_back(single_node);
		transit_nodes.push_back(node_counter);
		location_indices.push_back(node_counter);
	

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

		// Make the single cell

		NoCellCycleModel* p_cycle_model = new NoCellCycleModel();

		CellPtr p_cell(new Cell(p_state, p_cycle_model));
		p_cell->SetCellProliferativeType(p_diff_type);

		cells.push_back(p_cell);

		membraneSections.push_back(membrane_cells);

		NodeBasedCellPopulation<2> cell_population(mesh, cells, location_indices);

		OffLatticeSimulationTearOffStoppingEvent simulator(cell_population);

		simulator.SetOutputDirectory("TestSingleCellOnMembrane");

		simulator.SetEndTime(end_time);
		simulator.SetDt(dt);
		simulator.SetSamplingTimestepMultiple(sampling_multiple);

		MAKE_PTR(BasicNonLinearSpringForce<2>, p_force);


		p_force->SetEpithelialMembraneSpringStiffness(epithelialMembraneStiffness);

		p_force->SetEpithelialMembraneRestLength(membranePreferredRadius);

		p_force->SetEpithelialMembraneCutOffLength(membraneInteractionRadius);

		simulator.AddForce(p_force);

		MAKE_PTR_ARGS(CryptBoundaryCondition, p_bc, (&cell_population));
		simulator.AddCellPopulationBoundaryCondition(p_bc);

		MAKE_PTR(PushForce, p_push_force);
		p_push_force->SetCell(p_cell);
		c_vector<double, 2> force;
		force[0] = x_force;
		force[1] = 0;
		p_push_force->SetForce(force);
		simulator.AddForce(p_push_force);

		Node<2>* p_node =  cell_population.GetNodeCorrespondingToCell(p_cell);

		simulator.SetSingleNode(p_node);
		simulator.SetPushForce(force);

		double difference = 1;
		double tol = 1e-6;
		double iterations = 0;
		double it_limit = 100;

		double x_force_upper = 50;
		double x_force_lower = 0;

		bool failed = false;

		while (difference > tol && iterations < it_limit)
		{
			// Reset the simulation parameters
			force[0] = (x_force_upper + x_force_lower)/2;
			
			//simulator.RemoveAllForces();
			p_push_force->SetForce(force);
			//simulator.AddForce(p_push_force);
			simulator.SetPushForce(force);
			simulator.SetSimulationStartTime(SimulationTime::Instance()->GetTime());
			simulator.SetEndTime( SimulationTime::Instance()->GetTime() + end_time ); // Can't reset starting time, so just continue on
			c_vector<double, 2>& modifable_location = p_node->rGetModifiableLocation();
			modifable_location[0] = x_distance;
			modifable_location[1] = y_distance;
			// Solve again
			simulator.Solve();
			PRINT_VARIABLE(SimulationTime::Instance()->GetTime())
			
			c_vector<double, 2> end_force = p_node->rGetAppliedForce();
			PRINT_2_VARIABLES(end_force[0], end_force[1]);
			if (end_force[0] > .99 * force[0]) // If the force is significantly non-zero then the cell has broken off. It's usually exactly the push force
			{
				// This means the push force is too large, so make it the new upper bound
				TRACE("breaks")
				x_force_upper = force[0];
				
				
			}

			if (end_force[0] < .01) // If the force is very small, then the cell hasn't broken off. It's generally 
			{
				TRACE("too weak")
				x_force_lower = force[0];

			}

			if (x_force_upper != force[0] && x_force_lower != force[0])
			{
				failed = true;
				break;
			}

			iterations += 1;
			difference = x_force_upper - x_force_lower;
			TRACE(" ")
			PRINT_VARIABLE(x_force_upper)
			PRINT_VARIABLE(x_force_lower)
			PRINT_2_VARIABLES(end_force[0], end_force[1]);

		}

		if (failed)
		{
			TRACE("The forces at the end of the run were neither close to zero, or increasing")
			c_vector<double, 2> end_force = p_node->rGetAppliedForce();
			PRINT_2_VARIABLES(end_force[0], end_force[1]);
			PRINT_VARIABLE(difference)
		} else 
		{
			TRACE(" ")
			PRINT_VARIABLE(x_force_lower)
			PRINT_VARIABLE(iterations)
			PRINT_VARIABLE(difference)
		}

		ofstream myfile;
		myfile.open ("tear_off_force.txt", ios::app);
		myfile << x_distance << ","<< y_distance << ","<< membrane_spacing << ","<< epithelialMembraneStiffness << ","<< membraneInteractionRadius << ","<< membranePreferredRadius << "|" << x_force_lower << "\n";
		myfile.close();

	};

};
