/*

Copyright (c) 2005-2018, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#include "UniformContactInhibition.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "CellLabel.hpp"
#include "Debug.hpp"

UniformContactInhibition::UniformContactInhibition()
    : AbstractSimpleCellCycleModel(),
      mMinCellCycleDuration(12.0), // Hours
      mMaxCellCycleDuration(14.0)  // Hours
{
}

UniformContactInhibition::UniformContactInhibition(const UniformContactInhibition& rModel)
   : AbstractSimpleCellCycleModel(rModel),
     mMinCellCycleDuration(rModel.mMinCellCycleDuration),
     mMaxCellCycleDuration(rModel.mMaxCellCycleDuration)
{
    /*
     * Initialize only those member variables defined in this class.
     *
     * The member variable mCellCycleDuration is initialized in the
     * AbstractSimpleCellCycleModel constructor.
     *
     * The member variables mBirthTime, mReadyToDivide and mDimension
     * are initialized in the AbstractCellCycleModel constructor.
     *
     * Note that mCellCycleDuration is (re)set as soon as
     * InitialiseDaughterCell() is called on the new cell-cycle model.
     */
}

AbstractCellCycleModel* UniformContactInhibition::CreateCellCycleModel()
{
    return new UniformContactInhibition(*this);
}


bool UniformContactInhibition::ReadyToDivide()
{
    assert(mpCell != nullptr);


    if (!mReadyToDivide)
    {
        double cell_volume = mpCell->GetCellData()->GetItem("volume");

        //PRINT_VARIABLE(cell_volume)

        // For some reason I am unable to set the cell radius that the method GetItem("volume") returns
        // It always sets the radius as 0.5
        // Currently I have no work around, so we have to use 0.5

        double quiescent_volume = mEquilibriumVolume * mQuiescentVolumeFraction;
        
        if (GetAge() >= mCellCycleDuration && cell_volume > quiescent_volume)
        {
            mReadyToDivide = true;
        }
    }
    return mReadyToDivide;
};

void UniformContactInhibition::SetQuiescentVolumeFraction(double quiescentVolumeFraction)
{
    mQuiescentVolumeFraction = quiescentVolumeFraction;
};

void UniformContactInhibition::SetEquilibriumVolume(double equilibriumVolume)
{
    mEquilibriumVolume = equilibriumVolume;
};

void UniformContactInhibition::SetProliferativeRegion(double proliferativeRegion)
{
    mProliferativeRegion = proliferativeRegion;
};


void UniformContactInhibition::SetCellCycleDuration()
{
    RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();

    if (mpCell->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>())
    {
        mCellCycleDuration = DBL_MAX;
    }
    else
    {
        mCellCycleDuration = mMinCellCycleDuration + (mMaxCellCycleDuration - mMinCellCycleDuration) * p_gen->ranf(); // U[MinCCD,MaxCCD]
    }
}

double UniformContactInhibition::GetMinCellCycleDuration()
{
    return mMinCellCycleDuration;
}

void UniformContactInhibition::SetMinCellCycleDuration(double minCellCycleDuration)
{
    mMinCellCycleDuration = minCellCycleDuration;
}

double UniformContactInhibition::GetMaxCellCycleDuration()
{
    return mMaxCellCycleDuration;
}

void UniformContactInhibition::SetMaxCellCycleDuration(double maxCellCycleDuration)
{
    mMaxCellCycleDuration = maxCellCycleDuration;
}

double UniformContactInhibition::GetAverageTransitCellCycleTime()
{
    return 0.5*(mMinCellCycleDuration + mMaxCellCycleDuration);
}

double UniformContactInhibition::GetAverageStemCellCycleTime()
{
    return 0.5*(mMinCellCycleDuration + mMaxCellCycleDuration);
}

void UniformContactInhibition::OutputCellCycleModelParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<MinCellCycleDuration>" << mMinCellCycleDuration << "</MinCellCycleDuration>\n";
    *rParamsFile << "\t\t\t<MaxCellCycleDuration>" << mMaxCellCycleDuration << "</MaxCellCycleDuration>\n";

    // Call method on direct parent class
    AbstractSimpleCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(UniformContactInhibition)
