// Standard includes for tests
#include <cxxtest/TestSuite.h>
#include "AbstractCellBasedTestSuite.hpp"
#include "CheckpointArchiveTypes.hpp"
#include "SmartPointers.hpp"
#include "CommandLineArguments.hpp"

// Simulators and output
#include "OffLatticeSimulation.hpp" //Simulates the evolution of the population

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
#include "WntUniformCellCycleModel.hpp"
#include "WntCellCycleModelMembraneCell.hpp"

// Misc
#include "FakePetscSetup.hpp"
#include "Debug.hpp"

class TestContactNeighbour : public AbstractCellBasedTestSuite
{
	public:
	void TestContactNeighbours2DWithDivision() throw(Exception)
	{
		std::vector<Node<2>*> nodes;

		double dt = 0.05;
		double end_time = 20;
		double sampling_multiple = 1;

		double epithelialStiffness = 2.0; 			// 1.5

		double epithelialPreferredRadius = 0.7;			// 1.0

		double epithelialInteractionRadius = 1.5 * epithelialPreferredRadius; // Epithelial covers stem and transit

		HoneycombMeshGenerator generator(4,4);
        MutableMesh<2,2>* p_generating_mesh = generator.GetMesh();

        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(*p_generating_mesh, 1.5);

        std::vector<CellPtr> cells;
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        CellsGenerator<UniformCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_transit_type);

        NodeBasedCellPopulation<2> cell_population(mesh, cells);

        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestContactNeighbours2D");
        simulator.SetSamplingTimestepMultiple(sampling_multiple);
        simulator.SetEndTime(end_time);
        simulator.SetDt(dt);

        MAKE_PTR(LinearSpringForceMembraneCellNodeBased<2>, p_force);
		p_force->SetEpithelialSpringStiffness(epithelialStiffness);

		p_force->SetEpithelialPreferredRadius(epithelialPreferredRadius);

		p_force->SetEpithelialInteractionRadius(epithelialInteractionRadius);

		p_force->SetMeinekeSpringGrowthDuration(1);
		p_force->SetMeinekeDivisionRestingSpringLength(0.1);

		p_force->SetDebugMode(true);

        simulator.AddForce(p_force);

        simulator.Solve();
	};

};