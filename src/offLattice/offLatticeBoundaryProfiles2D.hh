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

#ifndef OFF_LATTICE_BOUNDARY_PROFILES_2D_HH
#define OFF_LATTICE_BOUNDARY_PROFILES_2D_HH

#include "core/globalDefs.h"
#include "offLattice/offLatticeBoundaryProfiles2D.h"

namespace plb {

template<typename T>
struct DefaultWallProfile2D<T, Array<T,2> > {
    BoundaryProfile2D<T,Array<T,2> >* generate() {
        return new NoSlipProfile2D<T>();
    }
};
#if 0
template<typename T>
struct DefaultWallProfile2D<T, Array<T,2> > {
    BoundaryProfile2D<T, Array<T,2> >* generate() {
        return new ScalarNeumannProfile2D<T>();
    }
};
#endif

/********** NoSlipProfile2D ********************************************/

template<typename T>
void NoSlipProfile2D<T>::setNormal(Array<T,2> const& normal_)
{ }

template<typename T>
void NoSlipProfile2D<T>::defineCircularShape(Array<T,2> const& center_, T radius_)
{ }

template<typename T>
void NoSlipProfile2D<T>::getData (
        Array<T,2> const& pos, plint id, AtomicBlock2D const* argument,
        Array<T,2>& data, OffBoundary::Type& bdType ) const
{
    data.resetToZero();
    bdType = OffBoundary::dirichlet;
}

template<typename T>
NoSlipProfile2D<T>*
    NoSlipProfile2D<T>::clone() const
{
    return new NoSlipProfile2D<T>(*this);
}


/********** FreeSlipProfile2D ********************************************/

template<typename T>
void FreeSlipProfile2D<T>::setNormal(Array<T,2> const& normal_)
{ }

template<typename T>
void FreeSlipProfile2D<T>::defineCircularShape(Array<T,2> const& center_, T radius_)
{ }

template<typename T>
void FreeSlipProfile2D<T>::getData (
        Array<T,2> const& pos, plint id, AtomicBlock2D const* argument,
        Array<T,2>& data, OffBoundary::Type& bdType ) const
{
    data.resetToZero();
    bdType = OffBoundary::freeSlip;
}

template<typename T>
FreeSlipProfile2D<T>*
    FreeSlipProfile2D<T>::clone() const
{
    return new FreeSlipProfile2D<T>(*this);
}


/********** PoiseuilleProfile2D ********************************************/

template<typename T>
PoiseuilleProfile2D<T>::PoiseuilleProfile2D(T uAverage_)
    : uAverage(uAverage_),
      normal(T(),T(),T()),
      center(T(),T(),T()),
      radius((T)1)
{ }

template<typename T>
void PoiseuilleProfile2D<T>::setNormal(Array<T,2> const& normal_) {
    normal = normal_;
}

template<typename T>
void PoiseuilleProfile2D<T>::defineCircularShape(Array<T,2> const& center_, T radius_)
{
    center = center_;
    radius = radius_;
}

template<typename T>
void PoiseuilleProfile2D<T>::getData (
        Array<T,2> const& pos, plint id, AtomicBlock2D const* argument,
        Array<T,2>& data, OffBoundary::Type& bdType ) const
{
    bdType = OffBoundary::dirichlet;
    Array<T,2> radial = pos-center;
    T r = norm(radial) / radius;
    if (r<=(T)1.) {
        data = 2*uAverage*(1-r*r)*(-normal);
    }
    else {
        data.resetToZero();
    }
}

template<typename T>
PoiseuilleProfile2D<T>*
    PoiseuilleProfile2D<T>::clone() const
{
    return new PoiseuilleProfile2D<T>(*this);
}


/********** IncreasingPoiseuilleProfile2D ********************************************/

template< typename T, template<typename U> class Descriptor>
IncreasingPoiseuilleProfile2D<T,Descriptor>::IncreasingPoiseuilleProfile2D (
        T uAverage_, plint maxT_ )
    : uAverage(uAverage_),
      maxT(maxT_),
      normal(T(),T(),T()),
      center(T(),T(),T()),
      radius((T)1)
{ }

template< typename T, template<typename U> class Descriptor>
void IncreasingPoiseuilleProfile2D<T,Descriptor>::setNormal(Array<T,2> const& normal_) {
    normal = normal_;
}

template< typename T, template<typename U> class Descriptor>
void IncreasingPoiseuilleProfile2D<T,Descriptor>::defineCircularShape (
        Array<T,2> const& center_, T radius_ )
{
    center = center_;
    radius = radius_;
}

template< typename T, template<typename U> class Descriptor>
void IncreasingPoiseuilleProfile2D<T,Descriptor>::getData (
        Array<T,2> const& pos, plint id, AtomicBlock2D const* argument,
        Array<T,2>& data, OffBoundary::Type& bdType ) const
{
    bdType = OffBoundary::dirichlet;
    BlockLattice2D<T,Descriptor> const* lattice =
        dynamic_cast<BlockLattice2D<T,Descriptor> const* >(argument);
    plint t = lattice->getTimeCounter().getTime();
    T signal=T();
    if (t<maxT) {
        signal = (T)0.5*((T)1+tanh((T)6/(T)maxT*((T)t-(T)maxT/(T)2)));
    }
    else {
        signal = (T)1;
    }
    Array<T,2> radial = pos-center;
    T r = norm(radial) / radius;
    if (r<=(T)1.) {
        data = signal*2*uAverage*(1-r*r)*(-normal);
    }
    else {
        data.resetToZero();
    }
}

template< typename T, template<typename U> class Descriptor>
IncreasingPoiseuilleProfile2D<T,Descriptor>*
    IncreasingPoiseuilleProfile2D<T,Descriptor>::clone() const
{
    return new IncreasingPoiseuilleProfile2D<T,Descriptor>(*this);
}


/********** IncreasingVelocityProfile2D ********************************************/

template< typename T, template<typename U> class Descriptor>
IncreasingVelocityProfile2D<T,Descriptor>::IncreasingVelocityProfile2D (
        Array<T,2> const& u_, plint maxT_ )
    : u(u_),
      maxT(maxT_)
{ }

template< typename T, template<typename U> class Descriptor>
void IncreasingVelocityProfile2D<T,Descriptor>::setNormal(Array<T,2> const& normal_)
{ }

template< typename T, template<typename U> class Descriptor>
void IncreasingVelocityProfile2D<T,Descriptor>::defineCircularShape (
        Array<T,2> const& center_, T radius_ )
{ }

template< typename T, template<typename U> class Descriptor>
void IncreasingVelocityProfile2D<T,Descriptor>::getData (
        Array<T,2> const& pos, plint id, AtomicBlock2D const* argument,
        Array<T,2>& data, OffBoundary::Type& bdType ) const
{
    bdType = OffBoundary::dirichlet;
    BlockLattice2D<T,Descriptor> const* lattice =
        dynamic_cast<BlockLattice2D<T,Descriptor> const* >(argument);
    plint t = lattice->getTimeCounter().getTime();
    T signal=T();
    if (t<maxT) {
        signal = (T)0.5*((T)1+tanh((T)6/(T)maxT*((T)t-(T)maxT/(T)2)));
    }
    else {
        signal = (T)1;
    }
    data = u*signal;
}

template< typename T, template<typename U> class Descriptor>
IncreasingVelocityProfile2D<T,Descriptor>*
    IncreasingVelocityProfile2D<T,Descriptor>::clone() const
{
    return new IncreasingVelocityProfile2D<T,Descriptor>(*this);
}


/********** TimeDependentVelocityProfile2D ********************************************/

template< typename T, template<typename U> class Descriptor>
TimeDependentVelocityProfile2D<T,Descriptor>::TimeDependentVelocityProfile2D (
        util::TimeDependentFunction<T,2>* velocity_)
    : velocity(velocity_)
{ }

template< typename T, template<typename U> class Descriptor>
TimeDependentVelocityProfile2D<T,Descriptor>::TimeDependentVelocityProfile2D (
        TimeDependentVelocityProfile2D<T,Descriptor> const& rhs)
    : velocity(rhs.velocity->clone())
{ }

template< typename T, template<typename U> class Descriptor>
TimeDependentVelocityProfile2D<T,Descriptor>& TimeDependentVelocityProfile2D<T,Descriptor>::operator=
        (TimeDependentVelocityProfile2D<T,Descriptor> const& rhs)
{
    TimeDependentVelocityProfile2D<T,Descriptor>(rhs).swap(*this);
}

template< typename T, template<typename U> class Descriptor>
void TimeDependentVelocityProfile2D<T,Descriptor>::swap(TimeDependentVelocityProfile2D<T,Descriptor>& rhs)
{
    std::swap(velocity, rhs.velocity);
}

template< typename T, template<typename U> class Descriptor>
void TimeDependentVelocityProfile2D<T,Descriptor>::setNormal(Array<T,2> const& normal_)
{ }

template< typename T, template<typename U> class Descriptor>
void TimeDependentVelocityProfile2D<T,Descriptor>::defineCircularShape (
        Array<T,2> const& center_, T radius_ )
{ }

template< typename T, template<typename U> class Descriptor>
void TimeDependentVelocityProfile2D<T,Descriptor>::getData (
        Array<T,2> const& pos, plint id, AtomicBlock2D const* argument,
        Array<T,2>& data, OffBoundary::Type& bdType ) const
{
    bdType = OffBoundary::dirichlet;
    BlockLattice2D<T,Descriptor> const* lattice =
        dynamic_cast<BlockLattice2D<T,Descriptor> const* >(argument);
    plint t = lattice->getTimeCounter().getTime();
    data = (*velocity)(t);
}

template< typename T, template<typename U> class Descriptor>
TimeDependentVelocityProfile2D<T,Descriptor>*
    TimeDependentVelocityProfile2D<T,Descriptor>::clone() const
{
    return new TimeDependentVelocityProfile2D<T,Descriptor>(*this);
}


/********** OscillatingPoiseuilleProfile2D ********************************************/

template< typename T, template<typename U> class Descriptor>
OscillatingPoiseuilleProfile2D<T,Descriptor>::OscillatingPoiseuilleProfile2D (
        T minUave_, T maxUave_, T period_ )
    : minUave(minUave_),
      maxUave(maxUave_),
      period(period_),
      normal(T(),T(),T()),
      center(T(),T(),T()),
      radius((T)1)
{ }

template< typename T, template<typename U> class Descriptor>
void OscillatingPoiseuilleProfile2D<T,Descriptor>::setNormal(Array<T,2> const& normal_)
{
    normal = normal_;
}

template< typename T, template<typename U> class Descriptor>
void OscillatingPoiseuilleProfile2D<T,Descriptor>::defineCircularShape (
        Array<T,2> const& center_, T radius_ )
{
    center = center_;
    radius = radius_;
}

template< typename T, template<typename U> class Descriptor>
void OscillatingPoiseuilleProfile2D<T,Descriptor>::getData (
        Array<T,2> const& pos, plint id, AtomicBlock2D const* argument,
        Array<T,2>& data, OffBoundary::Type& bdType ) const
{
    bdType = OffBoundary::dirichlet;
    static const T pi = (T)4.*atan(1.);
    BlockLattice2D<T,Descriptor> const* lattice =
        dynamic_cast<BlockLattice2D<T,Descriptor> const* >(argument);
    PLB_ASSERT(lattice);
    T t = (T) lattice->getTimeCounter().getTime();
    T signal = (sin((T)2.*pi*t/period)+(T)1.)*(T)0.5;

    Array<T,2> radial = pos-center;
    T r = norm(radial) / radius;
    if (r<=(T)1.) {
        data = ( minUave+(maxUave-minUave)*signal* (T)2*((T)1-r*r) ) * (-normal);
    }
    else {
        data.resetToZero();
    }
}

template< typename T, template<typename U> class Descriptor>
OscillatingPoiseuilleProfile2D<T,Descriptor>*
    OscillatingPoiseuilleProfile2D<T,Descriptor>::clone() const
{
    return new OscillatingPoiseuilleProfile2D<T,Descriptor>(*this);
}


/********** ConstantVelocityProfile2D ********************************************/

template<typename T>
ConstantVelocityProfile2D<T>::ConstantVelocityProfile2D(Array<T,2> const& u_)
    : u(u_)
{ }

template<typename T>
void ConstantVelocityProfile2D<T>::setNormal(Array<T,2> const& normal_)
{ }

template<typename T>
void ConstantVelocityProfile2D<T>::defineCircularShape(Array<T,2> const& center_, T radius_)
{ }

template<typename T>
void ConstantVelocityProfile2D<T>::getData (
        Array<T,2> const& pos, plint id, AtomicBlock2D const* argument,
        Array<T,2>& data, OffBoundary::Type& bdType ) const
{
    bdType = OffBoundary::dirichlet;
    data = u;
}

template<typename T>
ConstantVelocityProfile2D<T>*
    ConstantVelocityProfile2D<T>::clone() const
{
    return new ConstantVelocityProfile2D<T>(*this);
}


/********** VelocityPlugProfile2D ********************************************/

template<typename T>
VelocityPlugProfile2D<T>::VelocityPlugProfile2D(T uMax_)
    : uMax(uMax_),
      normal(T(),T(),T())
{ }

template<typename T>
void VelocityPlugProfile2D<T>::setNormal(Array<T,2> const& normal_) {
    normal = normal_;
}

template<typename T>
void VelocityPlugProfile2D<T>::defineCircularShape(Array<T,2> const& center_, T radius_)
{ }

template<typename T>
void VelocityPlugProfile2D<T>::getData (
        Array<T,2> const& pos, plint id, AtomicBlock2D const* argument,
        Array<T,2>& data, OffBoundary::Type& bdType ) const
{
    bdType = OffBoundary::dirichlet;
    data = -uMax*normal;
}

template<typename T>
VelocityPlugProfile2D<T>*
    VelocityPlugProfile2D<T>::clone() const
{
    return new VelocityPlugProfile2D<T>(*this);
}


/********** NeumannBoundaryProfile2D ******************************************/

template<typename T>
void NeumannBoundaryProfile2D<T>::setNormal(Array<T,2> const& normal_)
{
    normal = normal_;
}

template<typename T>
void NeumannBoundaryProfile2D<T>::defineCircularShape (
        Array<T,2> const& center_, T radius_ )
{ }

template<typename T>
void NeumannBoundaryProfile2D<T>::getData (
        Array<T,2> const& pos, plint id, AtomicBlock2D const* argument,
        Array<T,2>& data, OffBoundary::Type& bdType ) const
{
    bdType = OffBoundary::neumann;
    // In a Neumann condition, the velocity will need to be
    //   computed in the actual off-lattice boundary algorithm.
    //   Nothing to be done here.
    data.resetToZero();
}

template<typename T>
NeumannBoundaryProfile2D<T>*
    NeumannBoundaryProfile2D<T>::clone() const
{
    return new NeumannBoundaryProfile2D<T>(*this);
}


/********** DensityNeumannBoundaryProfile2D ******************************************/

template<typename T>
DensityNeumannBoundaryProfile2D<T>::DensityNeumannBoundaryProfile2D(T rho_)
    : rho(rho_)
{ }

template<typename T>
void DensityNeumannBoundaryProfile2D<T>::setNormal(Array<T,2> const& normal_)
{
    normal = normal_;
}

template<typename T>
void DensityNeumannBoundaryProfile2D<T>::defineCircularShape (
        Array<T,2> const& center_, T radius_ )
{ }

template<typename T>
void DensityNeumannBoundaryProfile2D<T>::getData (
        Array<T,2> const& pos, plint id, AtomicBlock2D const* argument,
        Array<T,2>& data, OffBoundary::Type& bdType ) const
{
    bdType = OffBoundary::densityNeumann;
    // In a Neumann condition, the velocity will need to be
    //   computed in the actual off-lattice boundary algorithm.
    //   Nothing to be done here.
    data.resetToZero();
    data[0] = rho;
}

template<typename T>
DensityNeumannBoundaryProfile2D<T>*
    DensityNeumannBoundaryProfile2D<T>::clone() const
{
    return new DensityNeumannBoundaryProfile2D<T>(*this);
}

/********** ScalarDirichletProfile2D ******************************************/

template<typename T>
ScalarDirichletProfile2D<T>::ScalarDirichletProfile2D(T value_)
    : value(value_)
{ }

template<typename T>
void ScalarDirichletProfile2D<T>::setNormal(Array<T,2> const& normal_)
{ }

template<typename T>
void ScalarDirichletProfile2D<T>::defineCircularShape (
        Array<T,2> const& center_, T radius_ )
{ }

template<typename T>
void ScalarDirichletProfile2D<T>::getData (
        Array<T,2> const& pos, plint id, AtomicBlock2D const* argument,
        Array<T,2>& data, OffBoundary::Type& bdType ) const
{
    bdType = OffBoundary::dirichlet;
    data[0] = value;
    data[1] = T();  // Second argument is not used here.
}

template<typename T>
ScalarDirichletProfile2D<T>*
    ScalarDirichletProfile2D<T>::clone() const
{
    return new ScalarDirichletProfile2D<T>(*this);
}


/********** ScalarNeumannProfile2D ******************************************/

template<typename T>
void ScalarNeumannProfile2D<T>::setNormal(Array<T,2> const& normal_)
{
    normal = normal_;
}

template<typename T>
void ScalarNeumannProfile2D<T>::defineCircularShape (
        Array<T,2> const& center_, T radius_ )
{ }

template<typename T>
void ScalarNeumannProfile2D<T>::getData (
        Array<T,2> const& pos, plint id, AtomicBlock2D const* argument,
        Array<T,2>& data, OffBoundary::Type& bdType ) const
{
    bdType = OffBoundary::neumann;
    // Neumann means zero-gradient, the two arguments are unused.
    data[0] = T();
    data[1] = T();
}

template<typename T>
ScalarNeumannProfile2D<T>*
    ScalarNeumannProfile2D<T>::clone() const
{
    return new ScalarNeumannProfile2D<T>(*this);
}



/********** ScalarFluxProfile2D ******************************************/

template<typename T>
ScalarFluxProfile2D<T>::ScalarFluxProfile2D(T gradVal_)
    : gradVal(gradVal_)
{ }

template<typename T>
void ScalarFluxProfile2D<T>::setNormal(Array<T,2> const& normal_)
{
    normal = normal_;
}

template<typename T>
void ScalarFluxProfile2D<T>::defineCircularShape (
        Array<T,2> const& center_, T radius_ )
{ }

template<typename T>
void ScalarFluxProfile2D<T>::getData (
        Array<T,2> const& pos, plint id, AtomicBlock2D const* argument,
        Array<T,2>& data, OffBoundary::Type& bdType ) const
{
    bdType = OffBoundary::flux;
    data[0] = gradVal;
    data[1] = T();
}

template<typename T>
ScalarFluxProfile2D<T>*
    ScalarFluxProfile2D<T>::clone() const
{
    return new ScalarFluxProfile2D<T>(*this);
}



/********** ScalarIsolationProfile2D ******************************************/

template<typename T>
ScalarIsolationProfile2D<T>::ScalarIsolationProfile2D(T asymptoticRho_, T kappa_)
    : asymptoticRho(asymptoticRho_),
      kappa(kappa_)
{ }

template<typename T>
void ScalarIsolationProfile2D<T>::setNormal(Array<T,2> const& normal_)
{
    normal = normal_;
}

template<typename T>
void ScalarIsolationProfile2D<T>::defineCircularShape (
        Array<T,2> const& center_, T radius_ )
{ }

template<typename T>
void ScalarIsolationProfile2D<T>::getData (
        Array<T,2> const& pos, plint id, AtomicBlock2D const* argument,
        Array<T,2>& data, OffBoundary::Type& bdType ) const
{
    bdType = OffBoundary::isolation;
    data[0] = asymptoticRho;
    data[1] = kappa;
}

template<typename T>
ScalarIsolationProfile2D<T>*
    ScalarIsolationProfile2D<T>::clone() const
{
    return new ScalarIsolationProfile2D<T>(*this);
}

}  // namespace plb

#endif  // OFF_LATTICE_BOUNDARY_PROFILES_2D_HH

