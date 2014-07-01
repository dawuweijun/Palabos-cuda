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

#ifndef OFF_LATTICE_BOUNDARY_CONDITION_2D_H
#define OFF_LATTICE_BOUNDARY_CONDITION_2D_H

#include "core/globalDefs.h"
#include "offLattice/triangularSurfaceMesh.h"
#include "offLattice/offLatticeBoundaryProfiles2D.h"
#include "offLattice/segmentBoundary2D.h"
#include "offLattice/triangleToDef.h"
#include "offLattice/guoOffLatticeModel2D.h"

namespace plb {

template< typename T,
          template<typename U> class Descriptor,
          class BoundaryType >
class OffLatticeBoundaryCondition2D {
public:
    OffLatticeBoundaryCondition2D (
            OffLatticeModel2D<T,BoundaryType>* offLatticeModel_,
            VoxelizedDomain2D<T>& voxelizedDomain_,
            MultiBlockLattice2D<T,Descriptor>& lattice_ );
    OffLatticeBoundaryCondition2D (
            OffLatticeModel2D<T,BoundaryType>* offLatticeModel_,
            VoxelizedDomain2D<T>& voxelizedDomain_,
            MultiBlockLattice2D<T,Descriptor>& lattice_,
            MultiParticleField2D<DenseParticleField2D<T,Descriptor> >& particleField_ );
    OffLatticeBoundaryCondition2D (
            OffLatticeBoundaryCondition2D<T,Descriptor,BoundaryType> const& rhs );
    ~OffLatticeBoundaryCondition2D();
    MultiBlockLattice2D<T,Descriptor> const& getLattice() const { return lattice; }
    VoxelizedDomain2D<T> const& getVoxelizedDomain() const { return voxelizedDomain; }
    VoxelizedDomain2D<T>& getVoxelizedDomain() { return voxelizedDomain; }
    void apply();
    void insert();
    void apply(std::vector<MultiBlock2D*> const& completionArg);
    void insert(std::vector<MultiBlock2D*> const& completionArg);
    Array<T,2> getForceOnObject();
    std::auto_ptr<MultiTensorField2D<T,2> > computeVelocity(Box2D domain);
    std::auto_ptr<MultiTensorField2D<T,2> > computeVelocity();
    std::auto_ptr<MultiTensorField2D<T,2> > computeVorticity(Box2D domain);
    std::auto_ptr<MultiTensorField2D<T,2> > computeVorticity();
    std::auto_ptr<MultiScalarField2D<T> > computeVelocityNorm(Box2D domain);
    std::auto_ptr<MultiScalarField2D<T> > computeVelocityNorm();
    std::auto_ptr<MultiScalarField2D<T> > computeVelocityComponent(Box2D domain, plint iComp);
    std::auto_ptr<MultiScalarField2D<T> > computeVelocityComponent(plint iComp);
    std::auto_ptr<MultiScalarField2D<T> > computePressure(Box2D domain);
    std::auto_ptr<MultiScalarField2D<T> > computePressure();
    std::auto_ptr<MultiScalarField2D<T> > computeDensity(Box2D domain, T solidDensity=T());
    std::auto_ptr<MultiScalarField2D<T> > computeDensity(T solidDensity=T());
    std::auto_ptr<MultiScalarField2D<T> > computeStrainRateNorm();
    std::auto_ptr<MultiScalarField2D<T> > computeStrainRateNorm(Box2D domain);
    std::auto_ptr<MultiTensorField2D<T,SymmetricTensor<T,Descriptor>::n> > computeStrainRate();
    std::auto_ptr<MultiTensorField2D<T,SymmetricTensor<T,Descriptor>::n> > computeStrainRate(Box2D domain);
    std::auto_ptr<MultiScalarField2D<T> > computeShearStressNorm();
    std::auto_ptr<MultiScalarField2D<T> > computeShearStressNorm(Box2D domain);
    std::auto_ptr<MultiTensorField2D<T,SymmetricTensor<T,Descriptor>::n> > computeShearStress();
    std::auto_ptr<MultiTensorField2D<T,SymmetricTensor<T,Descriptor>::n> > computeShearStress(Box2D domain);
    T computeAverageVelocityComponent(Box2D domain, plint iComponent);
    T computeAverageDensity(Box2D domain);
    T computeAverageDensity();
    T computeAverageEnergy(Box2D domain);
    T computeAverageEnergy();
    T computeRMSvorticity(Box2D domain);
    T computeRMSvorticity();
private:
    VoxelizedDomain2D<T>& voxelizedDomain;
    MultiBlockLattice2D<T,Descriptor>& lattice;
    MultiBlock2D& boundaryShapeArg;
    OffLatticeModel2D<T,BoundaryType>* offLatticeModel;
    MultiContainerBlock2D offLatticePattern;
};

}  // namespace plb

#endif  // OFF_LATTICE_BOUNDARY_CONDITION_2D_H
