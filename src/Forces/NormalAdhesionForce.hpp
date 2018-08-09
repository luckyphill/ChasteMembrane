// Provides a normal force to restrain the epithelial cells to the membrane
// Requires that the membrane cells are along the y axis

#ifndef NormalAdhesionForce_HPP_
#define NormalAdhesionForce_HPP_

#include "AbstractForce.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "DifferentiatedCellProliferativeType.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>


template<unsigned  ELEMENT_DIM, unsigned SPACE_DIM=ELEMENT_DIM>
class NormalAdhesionForce : public AbstractForce<ELEMENT_DIM, SPACE_DIM>
{
    friend class TestForces;

private:

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the object and its member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractForce<ELEMENT_DIM, SPACE_DIM> >(*this);
        archive & mStromalSpringStiffness; // Stromal is the differentiated "filler" cells
        archive & mMembraneStromalSpringStiffness;
    }

protected:


    double mStromalSpringStiffness; // Stromal is the differentiated "filler" cells

    double mMembraneStromalSpringStiffness;


    double mMembranePreferredRadius;
    double mStromalPreferredRadius; // Stromal is the differentiated "filler" cells

    double mMembraneInteractionRadius;
    double mStromalInteractionRadius; // Stromal is the differentiated "filler" cells


    bool mDebugMode = false;

public:

    /**
     * Constructor.
     */
    NormalAdhesionForce();

    /**
     * Destructor.
     */
    virtual ~NormalAdhesionForce();

    /**
     * Overridden CalculateForceBetweenNodes() method.
     *
     * Calculates the force between two nodes.
     *
     * Note that this assumes they are connected and is called by AddForceContribution()
     *
     * @param nodeAGlobalIndex index of one neighbouring node
     * @param nodeBGlobalIndex index of the other neighbouring node
     * @param rCellPopulation the cell population
     * @return The force exerted on Node A by Node B.
     */
    c_vector<double, SPACE_DIM> CalculateForceBetweenNodes(unsigned nodeAGlobalIndex,
                                                     unsigned nodeBGlobalIndex,
                                                     AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation);

    void AddForceContribution(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation);


    void SetMembraneSpringStiffness(double membraneSpringStiffness);
    void SetStromalSpringStiffness(double stromalSpringStiffness); // Stromal is the differentiated "filler" cells
    void SetMembraneStromalSpringStiffness(double membraneStromalSpringStiffness);

    void SetMembranePreferredRadius(double membranePreferredRadius);
    void SetStromalPreferredRadius(double stromalPreferredRadius); // Stromal is the differentiated "filler" cells

    void SetMembraneInteractionRadius(double membraneInteractionRadius);
    void SetStromalInteractionRadius(double stromalInteractionRadius); // Stromal is the differentiated "filler" cells

    
    /**
     * Overridden OutputForceParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */

    void SetDebugMode(bool debugStatus);
    
    virtual void OutputForceParameters(out_stream& rParamsFile);
};


#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(NormalAdhesionForce)

#endif /*NormalAdhesionForce_HPP_*/

#ifndef ND_SORT_FUNCTION
#define ND_SORT_FUNCTION
// Need to declare this sort function outide the class, otherwise it won't work
template<unsigned  ELEMENT_DIM, unsigned SPACE_DIM=ELEMENT_DIM>
bool nd_sort(std::tuple< Node<SPACE_DIM>*, c_vector<double, SPACE_DIM>, double > i,
                 std::tuple< Node<SPACE_DIM>*, c_vector<double, SPACE_DIM>, double > j)
{ 
    return (std::get<2>(i)<std::get<2>(j)); 
};
#endif