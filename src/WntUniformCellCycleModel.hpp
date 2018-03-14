

#ifndef WNTUNIFORMCELLCYCLEMODEL_HPP_
#define WNTUNIFORMCELLCYCLEMODEL_HPP_

#include "UniformCellCycleModel.hpp"
#include "WntConcentrationXSection.hpp"


class WntUniformCellCycleModel : public UniformCellCycleModel
{
	
public:
	WntUniformCellCycleModel();
	WntUniformCellCycleModel(const WntUniformCellCycleModel& rModel);

    /** Empty virtual destructor so archiving works with static libraries. */
    ~WntUniformCellCycleModel();

    bool IsAbovetWntThreshold();
    AbstractCellCycleModel* CreateCellCycleModel();
    void OutputCellCycleModelParameters(out_stream& rParamsFile);
    bool ReadyToDivide();

};

#include "SerializationExportWrapper.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(WntUniformCellCycleModel)

#endif /*WNTUNIFORMCELLCYCLEMODEL_HPP_*/