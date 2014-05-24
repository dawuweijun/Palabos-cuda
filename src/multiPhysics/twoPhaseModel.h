#ifndef TWOPHASEMODEL_H
#define TWOPHASEMODEL_H
#include "core/globalDefs.h"
#include "core/plbDebug.h"
namespace plb
{

/** kinetic: The pressure on the interface is entirely
 *               determined by the mass exchange between either phase.
 *  dynamic:     The pressure is considered to be a fluctuation around
 *               "outsideDensity", which can be either equal to rhoEmpty, or
 *               to the result of the pressure-correction model.
 *  constRho:    The pressure is constant throughout the bubble interface
 *               and takes its value from "outsideDensity", which can be either
 *               equal to rhoEmpty, or to the result of the pressure-correction model.
 **/
typedef enum {kinetic=1, dynamic=2, bubblePressure=3, constRho=4, freeSurface=5} TwoPhaseModel;

TwoPhaseModel stringToTwoPhaseModel ( std::string modelName )
{
    if ( modelName=="kinetic" )
    {
        return kinetic;
    }
    else if ( modelName=="dynamic" )
    {
        return dynamic;
    }
    else if ( modelName=="bubblePressure" )
    {
        return bubblePressure;
    }
    else if ( modelName=="constRho" )
    {
        return constRho;
    }
    else if ( modelName=="freeSurface" )
    {
        return freeSurface;
    }
    else
    {
        PLB_ASSERT ( false );
    }
};
/// Data structure for holding lists of cells along the free surface in an AtomicContainerBlock.
template< typename T,template<typename U> class Descriptor>
struct TwoPhaseInterfaceLists : public ContainerBlockData {
	typedef Array<plint,Descriptor<T>::d> Node;
	struct ExtrapolInfo {
		ExtrapolInfo() : density(T())
		{
			j.resetToZero(); PiNeq.resetToZero();
		}
		T density;
		Array<T,Descriptor<T>::d> j;
		Array<T,SymmetricTensor<T,Descriptor>::n> PiNeq;
	};
	/// Holds all nodes which have excess mass.
	std::map<Node,T> massExcess;
	/// Holds all nodes which have excess mass for fluid 2.
	std::map<Node,T> massExcess2;
	/// Holds all nodes that need to change status from interface to fluid.
	std::set<Node>   interfaceToFluid;
	/// Holds all nodes that need to change status from interface to empty.
	std::set<Node>   interfaceToEmpty;
	/// Holds all nodes that need to change status from empty to interface.
	std::map<Node,ExtrapolInfo> emptyToInterface;
	/// Holds all nodes that need to change status from fluid to interface.
	std::map<Node,ExtrapolInfo> fluidToInterface;

	virtual TwoPhaseInterfaceLists<T,Descriptor>* clone() const {
		return new TwoPhaseInterfaceLists<T,Descriptor>(*this);
	}
};


}//namespace plb
#endif
