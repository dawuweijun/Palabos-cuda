/* This file is part of the Palabos library.
 *
 * Copyright (C) 2011-2013 FlowKit Sarl
 * Route d'Oron 2
 * 1010 Lausanne, Switzerland
 * E-mail contact: contact@flowkit.com
 *
 * The most recent release of Palabos can be downloaded at 
 * <http://www.palabos.org/>
 *
 * The library Palabos is free software: you can redistribute it and/or
 * modify it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * The library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef OFF_LATTICE_BOUNDARY_PROCESSOR_2D_H
#define OFF_LATTICE_BOUNDARY_PROCESSOR_2D_H

#include "core/globalDefs.h"
#include "core/array.h"
#include "atomicBlock/dataProcessingFunctional2D.h"
#include "atomicBlock/atomicContainerBlock2D.h"
#include "dataProcessors/metaStuffHelper.h"
#include "offLattice/boundaryShapes2D.h"
#include <vector>
#include <utility>

namespace plb {

template<typename T>
class CheckVoxelizationFunctional2D : public PlainReductiveBoxProcessingFunctional2D
{
public:
    CheckVoxelizationFunctional2D(BoundaryShape2D<T,Array<T,3> >* shape_, int flowType_);
    virtual ~CheckVoxelizationFunctional2D();
    CheckVoxelizationFunctional2D(CheckVoxelizationFunctional2D const& rhs);
    CheckVoxelizationFunctional2D& operator= (
            CheckVoxelizationFunctional2D const& rhs );

    /// First AtomicBlock: Flag matrix.
    ///   If there are more atomic-blocks then they are forwarded to the
    ///   shape function, to provide additional read-only parameters.
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D*> fields);
    virtual CheckVoxelizationFunctional2D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
    plint getNumErrors() const;
private:
    bool isFluid(Dot2D const& location) const;
    void computeCell(Dot2D const& cellLocation,
                     ScalarField2D<int>& flagMatrix,
                     Dot2D const& offset);
private:
    plint numErrorsId;
    BoundaryShape2D<T,Array<T,3> >* shape;
    int flowType;
};

}  // namespace plb

#endif  // OFF_LATTICE_BOUNDARY_PROCESSOR_2D_H
