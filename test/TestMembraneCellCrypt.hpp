// Standard includes for tests
#include <cxxtest/TestSuite.h>
#include "AbstractCellBasedTestSuite.hpp"
#include "CheckpointArchiveTypes.hpp"
#include "SmartPointers.hpp"

// Mesh stuff
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "HoneycombMeshGenerator.hpp" //Generates mesh
#include "CylindricalHoneycombMeshGenerator.hpp"

// Simulators and output
#include "OffLatticeSimulation.hpp" //Simulates the evolution of the population
#include "VoronoiDataWriter.hpp" //Allows us to visualise output in Paraview

// Forces and BCs
#include "GeneralisedLinearSpringForce.hpp"
#include "MembraneCellForce.hpp"
#include "MembraneCellForceNodeBased.hpp"
#include "CryptBoundaryCondition.hpp"
#include "LinearSpringForceMembraneCell.hpp"
#include "LinearSpringForceMembraneCellNodeBased.hpp"

// Cell details
#include "MembraneCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "StemCellProliferativeType.hpp"
#include "BoundaryCellProperty.hpp"

#include "WildTypeCellMutationState.hpp"

//Cell cycle models
#include "NoCellCycleModel.hpp"
#include "UniformCellCycleModel.hpp"

// Misc
#include "AnoikisCellKillerMembraneCell.hpp"
#include "FakePetscSetup.hpp"
#include "Debug.hpp"


class TestMembraneCellCrypt : public AbstractCellBasedTestSuite
{
	public:
	void TestInsertCloseMembraneNodeBased() throw(Exception)
	{
		// In this we introduce a row of membrane point cells with a small rest length
		std::vector<Node<2>*> nodes;
		std::vector<unsigned> stem_nodes;
		std::vector<unsigned> transit_nodes;
		std::vector<unsigned> membrane_nodes;
		std::vector<unsigned> location_indices;
		std::vector<unsigned> ghost_nodes;
		std::vector<std::vector<CellPtr>> membraneSections;

		double dt = 0.01;
		double end_time = 50;
		double sampling_multiple = 10;

		unsigned cells_up = 40;
		unsigned cells_across = 40;
		unsigned ghosts = 3;
		unsigned node_counter = 0;
		unsigned num_membrane_nodes = 60;			// 60

		// Values that produce a working simulation in the comments
		double epithelialStiffness = 1.50; 			// 1.5
		double membraneStiffness = 0; 			// 5.0
		double stromalStiffness = 5.0; 				// 2.0

		double epithelialMembraneStiffness = 5.0; 	// 1.0
		double membraneStromalStiffness = 5.0; 		// 5.0
		double stromalEpithelialStiffness = 1.0;	// 1.0

		double epithelialPreferredRadius = 1.0;			// 1.0
		double membranePreferredRadius = 0.2;			// 0.2
		double stromalPreferredRadius = 0.5;			// 1.0

		double epithelialInteractionRadius = 3.0 * epithelialPreferredRadius; // Epithelial covers stem and transit
		double membraneInteractionRadius = 1.5 * membranePreferredRadius;
		double stromalInteractionRadius = 2.0 * stromalPreferredRadius; // Stromal is the differentiated "filler" cells

		double maxInteractionRadius = 2.5;

		double torsional_stiffness = 10;			// 10.0

		double targetCurvatureStemTrans = 0;
		double targetCurvatureTransTrans = 0;

		TRACE("The assertion for non-zero membrane spring force in LinearSpringForceMembraneCellNodeBased has been silenced at line 203 and 294")

		double centre_x = 5.0;
		double centre_y = 15.0;
		double base_radius = 3.0;
		double membrane_spacing = 0.2;
		double wall_height = 10;
		double left_side = centre_x - base_radius;
		double right_side = centre_x + base_radius;
		double wall_top = centre_y + wall_height;
		double wall_bottom = centre_y;

		double targetCurvatureStemStem = 1/base_radius;

		// Drawing the membrane
		// Need to follow this order to ensure membrane nodes are registered in order
		for (double y = wall_top; y >= wall_bottom; y-=membrane_spacing)
		{
			double x = left_side;
			nodes.push_back(new Node<2>(node_counter,  false,  x, y));
			membrane_nodes.push_back(node_counter);
			location_indices.push_back(node_counter);
			node_counter++;
		}

		double acos_arg = (2*pow(base_radius,2) - pow(membrane_spacing,2))/(2*pow(base_radius,2));
		double  d_theta = acos(acos_arg);
		PRINT_VARIABLE(acos_arg)

		for (double theta = d_theta; theta <= M_PI; theta += d_theta)
		{
			double x = centre_x - base_radius * cos(theta);
			double y = centre_y - base_radius * sin(theta);
			nodes.push_back(new Node<2>(node_counter,  false,  x, y));
			membrane_nodes.push_back(node_counter);
			location_indices.push_back(node_counter);
			node_counter++;
		}

		for (double y = wall_bottom; y <= wall_top; y+=membrane_spacing)
		{
			double x = right_side;
			nodes.push_back(new Node<2>(node_counter,  false,  x, y));
			membrane_nodes.push_back(node_counter);
			location_indices.push_back(node_counter);
			node_counter++;
		}

		//Drawing the epithelium
		//Start with the stem cells
		base_radius -= epithelialPreferredRadius;
		acos_arg = (2*pow(base_radius,2) - pow(epithelialPreferredRadius,2))/(2*pow(base_radius,2));
		d_theta = acos(acos_arg);
		PRINT_VARIABLE(acos_arg)

		for (double theta = d_theta; theta <= M_PI; theta += d_theta)
		{
			double x = centre_x - base_radius * cos(theta);
			double y = centre_y - base_radius * sin(theta);
			nodes.push_back(new Node<2>(node_counter,  false,  x, y));
			stem_nodes.push_back(node_counter);
			location_indices.push_back(node_counter);
			node_counter++;
		}
		// Next the transit amplifying cells
		for (double y = wall_bottom; y <= wall_top; y+= epithelialPreferredRadius)
		{
			double x = left_side + epithelialPreferredRadius;
			nodes.push_back(new Node<2>(node_counter,  false,  x, y));
			transit_nodes.push_back(node_counter);
			location_indices.push_back(node_counter);
			node_counter++;

			x = right_side - epithelialPreferredRadius;
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
		MAKE_PTR(StemCellProliferativeType, p_stem_type);
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

		//Initialise stem nodes
		for (unsigned i = 0; i < stem_nodes.size(); i++)
		{
			UniformCellCycleModel* p_cycle_model = new UniformCellCycleModel();
			double birth_time = 12.0*RandomNumberGenerator::Instance()->ranf(); //Randomly set birth time to stop pulsing behaviour
			p_cycle_model->SetBirthTime(-birth_time);

			CellPtr p_cell(new Cell(p_state, p_cycle_model));
			p_cell->SetCellProliferativeType(p_stem_type);

			p_cell->InitialiseCellCycleModel();

			cells.push_back(p_cell);
		}

		//Initialise trans nodes
		for (unsigned i = 0; i < transit_nodes.size(); i++)
		{
			UniformCellCycleModel* p_cycle_model = new UniformCellCycleModel();
			double birth_time = 12.0*RandomNumberGenerator::Instance()->ranf(); //Randomly set birth time to stop pulsing behaviour
			p_cycle_model->SetBirthTime(-birth_time);

			CellPtr p_cell(new Cell(p_state, p_cycle_model));
			p_cell->SetCellProliferativeType(p_trans_type);

			p_cell->InitialiseCellCycleModel();

			cells.push_back(p_cell);
		}

		membraneSections.push_back(membrane_cells);

		NodeBasedCellPopulation<2> cell_population(mesh, cells, location_indices);

		OffLatticeSimulation<2> simulator(cell_population);
		simulator.SetOutputDirectory("MembraneCellCryptNodeBased");
		simulator.SetEndTime(end_time);
		simulator.SetDt(dt);
		simulator.SetSamplingTimestepMultiple(sampling_multiple);

		MAKE_PTR(LinearSpringForceMembraneCellNodeBased<2>, p_force);
		p_force->SetEpithelialSpringStiffness(epithelialStiffness);
		p_force->SetMembraneSpringStiffness(membraneStiffness);
		p_force->SetStromalSpringStiffness(stromalStiffness);
		p_force->SetEpithelialMembraneSpringStiffness(epithelialMembraneStiffness);
		p_force->SetMembraneStromalSpringStiffness(membraneStromalStiffness);
		p_force->SetStromalEpithelialSpringStiffness(stromalEpithelialStiffness);

		p_force->SetEpithelialPreferredRadius(epithelialPreferredRadius);
		p_force->SetMembranePreferredRadius(membranePreferredRadius);
		p_force->SetStromalPreferredRadius(stromalPreferredRadius);

		p_force->SetEpithelialInteractionRadius(epithelialInteractionRadius);
		p_force->SetMembraneInteractionRadius(membraneInteractionRadius);
		p_force->SetStromalInteractionRadius(stromalInteractionRadius);

		simulator.AddForce(p_force);

		MAKE_PTR(MembraneCellForceNodeBased, p_membrane_force);
		p_membrane_force->SetBasementMembraneTorsionalStiffness(torsional_stiffness);
		p_membrane_force->SetTargetCurvatures(targetCurvatureStemStem, targetCurvatureStemTrans, targetCurvatureTransTrans);
		p_membrane_force->SetMembraneSections(membraneSections);
		//p_membrane_force->SetCalculationToTorsion(true);
		//simulator.AddForce(p_membrane_force);

		MAKE_PTR_ARGS(CryptBoundaryCondition, p_bc, (&cell_population));
		simulator.AddCellPopulationBoundaryCondition(p_bc);

		MAKE_PTR_ARGS(AnoikisCellKillerMembraneCell, p_anoikis_killer, (&cell_population));
		simulator.AddCellKiller(p_anoikis_killer);

		simulator.Solve();
	};

};