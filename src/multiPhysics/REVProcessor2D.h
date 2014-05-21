#ifndef REV_PROCESSORS_2D_H
#define REV_PROCESSORS_2D_H

#include "core/globalDefs.h"
#include "core/block2D.h"
#include "atomicBlock/dataProcessor2D.h"
#include "atomicBlock/blockLattice2D.h"
#include "multiBlock/multiDataField2D.h"
namespace plb
{
template<typename T1, template<typename U> class Descriptor, typename T2>
class BrinkmanProcessor2D: public BoxProcessingFunctional2D_LT<T1,Descriptor,T2,4>
{
public:
    virtual void process ( Box2D domain, BlockLattice2D<T1,Descriptor>& lattice,
                           TensorField2D<T2,4>& negNiuInvsK );
    virtual BrinkmanProcessor2D<T1,Descriptor,T2> *clone() const;
    virtual void getTypeOfModification ( std::vector<modif::ModifT>& modified ) const
    {
        modified[0] = modif::staticVariables;
        modified[1] = modif::staticVariables;
    }

};


} //namespace plb

#endif
