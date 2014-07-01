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

#ifndef OFF_LATTICE_BOUNDARY_PROFILES_2D_H
#define OFF_LATTICE_BOUNDARY_PROFILES_2D_H

#include "core/globalDefs.h"
#include "offLattice/triangularSurfaceMesh.h"
#include "algorithm/functions.h"

namespace plb {


template<typename T, class SurfaceData>
class BoundaryProfile2D {
public:
    virtual ~BoundaryProfile2D() { }
    virtual void setNormal(Array<T,2> const& normal_) =0;
    virtual void defineCircularShape(Array<T,2> const& radius_, T center_) =0;
    virtual void getData( Array<T,2> const& pos, plint id, AtomicBlock2D const* argument,
                          SurfaceData& data, OffBoundary::Type& bdType ) const =0;
    virtual BoundaryProfile2D<T,SurfaceData>* clone() const =0;
};

template<typename T, class SurfaceData>
struct DefaultWallProfile2D {
    BoundaryProfile2D<T,SurfaceData>* generate() {
        // A default wall profile needs yet to be implemented for this case.
        PLB_ASSERT( false );
        return 0;
    }
}; 

template<typename T, class SurfaceData>
BoundaryProfile2D<T, SurfaceData>* generateDefaultWallProfile2D() {
    return DefaultWallProfile2D<T,SurfaceData>().generate();
}

template<typename T>
class NoSlipProfile2D : public BoundaryProfile2D<T, Array<T,2> >
{
public:
    virtual void setNormal(Array<T,2> const& normal_);
    virtual void defineCircularShape(Array<T,2> const& center_, T radius_);
    virtual void getData( Array<T,2> const& pos, plint id, AtomicBlock2D const* argument,
                          Array<T,2>& data, OffBoundary::Type& bdType ) const;
    virtual NoSlipProfile2D<T>* clone() const;
};

template<typename T>
class FreeSlipProfile2D : public BoundaryProfile2D<T, Array<T,2> >
{
public:
    virtual void setNormal(Array<T,2> const& normal_);
    virtual void defineCircularShape(Array<T,2> const& center_, T radius_);
    virtual void getData( Array<T,2> const& pos, plint id, AtomicBlock2D const* argument,
                          Array<T,2>& data, OffBoundary::Type& bdType ) const;
    virtual FreeSlipProfile2D<T>* clone() const;
};

template<typename T>
class ConstantVelocityProfile2D : public BoundaryProfile2D<T, Array<T,2> >
{
public:
    ConstantVelocityProfile2D(Array<T,2> const& u_);
    virtual void setNormal(Array<T,2> const& normal_);
    virtual void defineCircularShape(Array<T,2> const& center_, T radius_);
    virtual void getData( Array<T,2> const& pos, plint id, AtomicBlock2D const* argument,
                          Array<T,2>& data, OffBoundary::Type& bdType ) const;
    virtual ConstantVelocityProfile2D<T>* clone() const;
private:
    Array<T,2> u;
};

template<typename T>
class VelocityPlugProfile2D : public BoundaryProfile2D<T, Array<T,2> >
{
public:
    VelocityPlugProfile2D(T uMax_);
    virtual void setNormal(Array<T,2> const& normal_);
    virtual void defineCircularShape(Array<T,2> const& center_, T radius_);
    virtual void getData( Array<T,2> const& pos, plint id, AtomicBlock2D const* argument,
                          Array<T,2>& data, OffBoundary::Type& bdType ) const;
    virtual VelocityPlugProfile2D<T>* clone() const;
private:
    T uMax;
    Array<T,2> normal;
};

template< typename T, template<typename U> class Descriptor>
class OscillatingPoiseuilleProfile2D : public BoundaryProfile2D<T, Array<T,2> >
{
public:
    OscillatingPoiseuilleProfile2D(T minUave_, T maxUave_, T period_);
    virtual void setNormal(Array<T,2> const& normal_);
    virtual void defineCircularShape(Array<T,2> const& center_, T radius_);
    virtual void getData( Array<T,2> const& pos, plint id, AtomicBlock2D const* argument,
                          Array<T,2>& data, OffBoundary::Type& bdType ) const;
    virtual OscillatingPoiseuilleProfile2D<T,Descriptor>* clone() const;
private:
    T minUave, maxUave, period;
    Array<T,2> normal;
    Array<T,2> center;
    T radius;
};

template< typename T, template<typename U> class Descriptor>
class IncreasingPoiseuilleProfile2D : public BoundaryProfile2D<T, Array<T,2> >
{
public:
    IncreasingPoiseuilleProfile2D(T uAverage_, plint maxT_);
    virtual void setNormal(Array<T,2> const& normal_);
    virtual void defineCircularShape(Array<T,2> const& center_, T radius_);
    virtual void getData( Array<T,2> const& pos, plint id, AtomicBlock2D const* argument,
                          Array<T,2>& data, OffBoundary::Type& bdType ) const;
    virtual IncreasingPoiseuilleProfile2D<T,Descriptor>* clone() const;
private:
    T uAverage;
    plint maxT;
    Array<T,2> normal;
    Array<T,2> center;
    T radius;
};

template< typename T, template<typename U> class Descriptor>
class IncreasingVelocityProfile2D : public BoundaryProfile2D<T, Array<T,2> >
{
public:
    IncreasingVelocityProfile2D(Array<T,2> const& u_, plint maxT_);
    virtual void setNormal(Array<T,2> const& normal_);
    virtual void defineCircularShape(Array<T,2> const& center_, T radius_);
    virtual void getData( Array<T,2> const& pos, plint id, AtomicBlock2D const* argument,
                          Array<T,2>& data, OffBoundary::Type& bdType ) const;
    virtual IncreasingVelocityProfile2D<T,Descriptor>* clone() const;
private:
    Array<T,2> u;
    plint maxT;
};

template< typename T, template<typename U> class Descriptor>
class TimeDependentVelocityProfile2D : public BoundaryProfile2D<T, Array<T,2> >
{
public:
    TimeDependentVelocityProfile2D(util::TimeDependentFunction<T,2>* velocity_);
    TimeDependentVelocityProfile2D(TimeDependentVelocityProfile2D<T,Descriptor> const& rhs);
    TimeDependentVelocityProfile2D<T,Descriptor>& operator= (
            TimeDependentVelocityProfile2D<T,Descriptor> const& rhs );
    void swap(TimeDependentVelocityProfile2D<T,Descriptor>& rhs);
    virtual void setNormal(Array<T,2> const& normal_);
    virtual void defineCircularShape(Array<T,2> const& center_, T radius_);
    virtual void getData( Array<T,2> const& pos, plint id, AtomicBlock2D const* argument,
                          Array<T,2>& data, OffBoundary::Type& bdType ) const;
    virtual TimeDependentVelocityProfile2D<T,Descriptor>* clone() const;
private:
    util::TimeDependentFunction<T,2>* velocity;
};

template<typename T>
class PoiseuilleProfile2D : public BoundaryProfile2D<T, Array<T,2> >
{
public:
    PoiseuilleProfile2D(T uAverage_);
    virtual void setNormal(Array<T,2> const& normal_);
    virtual void defineCircularShape(Array<T,2> const& center_, T radius_);
    virtual void getData( Array<T,2> const& pos, plint id, AtomicBlock2D const* argument,
                          Array<T,2>& data, OffBoundary::Type& bdType ) const;
    virtual PoiseuilleProfile2D<T>* clone() const;
private:
    T uAverage;
    Array<T,2> normal;
    Array<T,2> center;
    T radius;
};

template<typename T>
class NeumannBoundaryProfile2D : public BoundaryProfile2D<T, Array<T,2> >
{
public:
    virtual void setNormal(Array<T,2> const& normal_);
    virtual void defineCircularShape(Array<T,2> const& center_, T radius_);
    virtual void getData( Array<T,2> const& pos, plint id, AtomicBlock2D const* argument,
                          Array<T,2>& data, OffBoundary::Type& bdType ) const;
    virtual NeumannBoundaryProfile2D<T>* clone() const;
private:
    Array<T,2> normal;
};

template<typename T>
class DensityNeumannBoundaryProfile2D : public BoundaryProfile2D<T, Array<T,2> >
{
public:
    DensityNeumannBoundaryProfile2D(T rho_ = (T)1);
    virtual void setNormal(Array<T,2> const& normal_);
    virtual void defineCircularShape(Array<T,2> const& center_, T radius_);
    virtual void getData( Array<T,2> const& pos, plint id, AtomicBlock2D const* argument,
                          Array<T,2>& data, OffBoundary::Type& bdType ) const;
    virtual DensityNeumannBoundaryProfile2D<T>* clone() const;
private:
    Array<T,2> normal;
    T rho;
};

template<typename T>
class ScalarNeumannProfile2D : public BoundaryProfile2D<T,Array<T,2> >
{
public:
    virtual void setNormal(Array<T,2> const& normal_);
    virtual void defineCircularShape(Array<T,2> const& center_, T radius_);
    virtual void getData( Array<T,2> const& pos, plint id, AtomicBlock2D const* argument,
                          Array<T,2>& data, OffBoundary::Type& bdType ) const;
    virtual ScalarNeumannProfile2D<T>* clone() const;
private:
    Array<T,2> normal;
};

template<typename T>
class ScalarDirichletProfile2D : public BoundaryProfile2D<T,Array<T,2> >
{
public:
    ScalarDirichletProfile2D(T value_);
    virtual void setNormal(Array<T,2> const& normal_);
    virtual void defineCircularShape(Array<T,2> const& center_, T radius_);
    virtual void getData( Array<T,2> const& pos, plint id, AtomicBlock2D const* argument,
                          Array<T,2>& data, OffBoundary::Type& bdType ) const;
    virtual ScalarDirichletProfile2D<T>* clone() const;
private:
    T value;
};

template<typename T>
class ScalarFluxProfile2D : public BoundaryProfile2D<T,Array<T,2> >
{
public:
    ScalarFluxProfile2D(T gradVal_);
    virtual void setNormal(Array<T,2> const& normal_);
    virtual void defineCircularShape(Array<T,2> const& center_, T radius_);
    virtual void getData( Array<T,2> const& pos, plint id, AtomicBlock2D const* argument,
                          Array<T,2>& data, OffBoundary::Type& bdType ) const;
    virtual ScalarFluxProfile2D<T>* clone() const;
private:
    Array<T,2> normal;
    T gradVal;
};

/// Implements the condition grad(rho) = kappa(asymptoticRho-rho).
template<typename T>
class ScalarIsolationProfile2D : public BoundaryProfile2D<T,Array<T,2> >
{
public:
    ScalarIsolationProfile2D(T asymptoticRho_, T kappa_);
    virtual void setNormal(Array<T,2> const& normal_);
    virtual void defineCircularShape(Array<T,2> const& center_, T radius_);
    virtual void getData( Array<T,2> const& pos, plint id, AtomicBlock2D const* argument,
                          Array<T,2>& data, OffBoundary::Type& bdType ) const;
    virtual ScalarIsolationProfile2D<T>* clone() const;
private:
    Array<T,2> normal;
    T asymptoticRho, kappa;
};

}  // namespace plb

#endif  // OFF_LATTICE_BOUNDARY_PROFILES_2D_H

