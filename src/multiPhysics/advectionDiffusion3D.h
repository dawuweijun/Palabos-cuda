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

#ifndef ADVECTION_DIFFUSION_3D_H
#define ADVECTION_DIFFUSION_3D_H

#include "core/globalDefs.h"
#include "core/block3D.h"
#include "atomicBlock/dataProcessor3D.h"
#include "multiBlock/multiBlockLattice3D.h"
#include "multiBlock/multiDataField3D.h"

namespace plb {

/**
* Multiphysics class for one-way coupling between Navier-Stokes and
* advection-diffusion equations: the fluid velocity is copied
* to the advection-diffusion field, which is advected passively.
*/
template< typename T, template<typename U> class TemperatureDescriptor >
class VelocityToPassiveAdvDiff3D :
    public BoxProcessingFunctional3D_LT<T,TemperatureDescriptor,T,3>
{
public:
    
    virtual void process( Box3D domain,
                          BlockLattice3D<T,TemperatureDescriptor>& temperature,
                          TensorField3D<T,3>& velocity );
    virtual VelocityToPassiveAdvDiff3D<T,TemperatureDescriptor>* clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
};

template< typename T, template<typename U> class TemperatureDescriptor >
class N_VelocityToPassiveAdvDiff3D :
    public BoxProcessingFunctional3D_LN<T,TemperatureDescriptor,T>
{
public:
    
    virtual void process( Box3D domain,
                          BlockLattice3D<T,TemperatureDescriptor>& temperature,
                          NTensorField3D<T>& velocity );
    virtual N_VelocityToPassiveAdvDiff3D<T,TemperatureDescriptor>* clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
};


/**
* Multiphysics class for one-way coupling between Navier-Stokes and
* advection-diffusion equations: the fluid velocity is copied
* to the advection-diffusion field, which is advected passively.
*/
template< typename T,
          template<typename U1> class FluidDescriptor, 
          template<typename U2> class ScalarDescriptor
        >
class LatticeToPassiveAdvDiff3D :
    public BoxProcessingFunctional3D_LL<T,FluidDescriptor,T,ScalarDescriptor>
{
public:
    LatticeToPassiveAdvDiff3D(T scaling_=1.);
    virtual void process( Box3D domain,
                          BlockLattice3D<T,FluidDescriptor>& fluid,
                          BlockLattice3D<T,ScalarDescriptor>& scalar );
    virtual LatticeToPassiveAdvDiff3D<T,FluidDescriptor,ScalarDescriptor>* clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
private:
    T scaling;
};

template< typename T, template<typename U> class Descriptor >
class CrystallizeAndAggregate :
    public BoxProcessingFunctional3D_L<T,Descriptor>
{
public:
    CrystallizeAndAggregate(T Ncr_, T Nag_);
    virtual void process( Box3D domain,
                          BlockLattice3D<T,Descriptor>& lattice );
    virtual CrystallizeAndAggregate<T,Descriptor>* clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
private:
    T Ncr, Nag;
};

template< typename T, template<typename U> class Descriptor >
void crystallizeAndAggregate(MultiBlockLattice3D<T,Descriptor>& lattice, T Ncr, T Nag, Box3D domain);

template< typename T,
          template<typename U1> class FluidDescriptor, 
          template<typename U2> class ScalarDescriptor
        >
void latticeToPassiveAdvDiff(MultiBlockLattice3D<T,FluidDescriptor>& fluid, MultiBlockLattice3D<T,ScalarDescriptor>& scalar, Box3D domain);

template< typename T, template<typename U> class TemperatureDescriptor >
void NVelocityToPassiveAdvDiff(MultiBlockLattice3D<T,TemperatureDescriptor>& temperature, MultiNTensorField3D<T>& velocity, Box3D domain);

}  // namespace plb

#endif  // ADVECTION_DIFFUSION_3D_H

