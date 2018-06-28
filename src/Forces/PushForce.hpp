#ifndef PUSHFORCE_HPP_
#define PUSHFORCE_HPP_

#include "ChasteSerialization.hpp"
#include "ClassIsAbstract.hpp"

#include "AbstractCellPopulation.hpp"
#include "AbstractCellBasedSimulationModifier.hpp"
#include "AbstractForce.hpp"

class PushForce: public AbstractForce<2,2>
{
	private:

		CellPtr mpCell;
		c_vector<double, 2> mForce;

	public:

		PushForce();
		~PushForce();


		void SetCell(CellPtr cell);

		void SetForce(c_vector<double, 2> force);

		void AddForceContribution(AbstractCellPopulation<2>& rCellPopulation);

		void OutputForceParameters(out_stream& rParamsFile);

};

#endif /*PUSHFORCE_HPP_*/