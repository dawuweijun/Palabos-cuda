#ifndef REV_PROCESSORS_3D_H
#define REV_PROCESSORS_3D_H

#include "core/globalDefs.h"
#include "core/block2D.h"
#include "atomicBlock/dataProcessor3D.h"
#include "atomicBlock/blockLattice3D.h"
#include "multiBlock/multiDataField3D.h"
namespace plb{
template<typename T1, template<typename U> class Descriptor, typename T2>
class BrinkmanProcessor3D: public BoxProcessingFunctional3D_LT<T1,Descriptor,T2,9>
{
public:
    virtual void process ( Box3D domain, BlockLattice3D<T1,Descriptor>& lattice,
                           TensorField3D<T2,9>& negNiuInvsK );
    virtual BrinkmanProcessor3D<T1,Descriptor,T2> *clone() const;
    virtual void getTypeOfModification ( std::vector<modif::ModifT>& modified ) const
    {
        modified[0] = modif::staticVariables;
        modified[1] = modif::staticVariables;
    }

};


}// namespace plb
#endif
