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

#ifndef FREE_SURFACE_INITIALIZER_2D_H
#define FREE_SURFACE_INITIALIZER_2D_H

#include "core/globalDefs.h"
#include "atomicBlock/dataProcessingFunctional2D.h"
#include "multiPhysics/freeSurfaceUtil2D.h"

namespace plb {

template< typename T,template<typename U> class Descriptor>
class DefaultInitializeFreeSurface2D : public BoxProcessingFunctional2D {
public:
    DefaultInitializeFreeSurface2D(Dynamics<T,Descriptor>* dynamicsTemplate_, Array<T,2> g_, T rhoIni_=(T)1.)
        : dynamicsTemplate(dynamicsTemplate_), g(g_), rhoIni(rhoIni_)
    { }
    DefaultInitializeFreeSurface2D(DefaultInitializeFreeSurface2D<T,Descriptor> const& rhs)
        : dynamicsTemplate(rhs.dynamicsTemplate->clone()),
          g(rhs.g), rhoIni(rhs.rhoIni)
    { }
    DefaultInitializeFreeSurface2D<T,Descriptor>* operator=(DefaultInitializeFreeSurface2D<T,Descriptor> const& rhs)
    { 
        DefaultInitializeFreeSurface2D<T,Descriptor>(rhs).swap(*this);
        return *this;
    }
    void swap(DefaultInitializeFreeSurface2D<T,Descriptor>& rhs) {
        std::swap(dynamicsTemplate, rhs.dynamicsTemplate);
        std::swap(g, rhs.g);
        std::swap(rhoIni, rhs.rhoIni);
    }
    virtual ~DefaultInitializeFreeSurface2D() {
        delete dynamicsTemplate;
    }
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D*> atomicBlocks);
    virtual DefaultInitializeFreeSurface2D<T,Descriptor>* clone() const {
        return new DefaultInitializeFreeSurface2D(*this);
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        std::fill(modified.begin(), modified.end(), modif::nothing);
        modified[0] = modif::dataStructure;    // Fluid
        modified[1] = modif::staticVariables;  // rhoBar.
        modified[2] = modif::staticVariables;  // j.
        modified[3] = modif::staticVariables;  // Mass
        modified[4] = modif::staticVariables;  // Volume-fraction
        modified[5] = modif::staticVariables;  // Flag-status
        modified[6] = modif::nothing;          // Normal.
        modified[7] = modif::nothing;          // Interface lists
        modified[8] = modif::nothing;          // Curvature.
        modified[9] = modif::nothing;          // Outside density.
    }
private:
    Dynamics<T,Descriptor>* dynamicsTemplate;
    Array<T,2> g;
    T rhoIni;
};

// Same as DefaultInitializeFreeSurface, but without initializing the Volume-fraction.
template< typename T,template<typename U> class Descriptor>
class PartiallyDefaultInitializeFreeSurface2D : public BoxProcessingFunctional2D {
public:
    PartiallyDefaultInitializeFreeSurface2D(Dynamics<T,Descriptor>* dynamicsTemplate_, Array<T,2> g_, T rhoIni_=(T)1.)
        : dynamicsTemplate(dynamicsTemplate_), g(g_), rhoIni(rhoIni_)
    { }
    PartiallyDefaultInitializeFreeSurface2D(PartiallyDefaultInitializeFreeSurface2D<T,Descriptor> const& rhs)
        : dynamicsTemplate(rhs.dynamicsTemplate->clone()),
          g(rhs.g),
          rhoIni(rhs.rhoIni)
    { }
    PartiallyDefaultInitializeFreeSurface2D<T,Descriptor>* operator= (
            PartiallyDefaultInitializeFreeSurface2D<T,Descriptor> const& rhs )
    { 
        PartiallyDefaultInitializeFreeSurface2D<T,Descriptor>(rhs).swap(*this);
        return *this;
    }
    void swap(PartiallyDefaultInitializeFreeSurface2D<T,Descriptor>& rhs) {
        std::swap(dynamicsTemplate, rhs.dynamicsTemplate);
        std::swap(g, rhs.g);
        std::swap(rhoIni, rhs.rhoIni);
    }
    virtual ~PartiallyDefaultInitializeFreeSurface2D() {
        delete dynamicsTemplate;
    }
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D*> atomicBlocks);
    virtual PartiallyDefaultInitializeFreeSurface2D<T,Descriptor>* clone() const {
        return new PartiallyDefaultInitializeFreeSurface2D(*this);
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        std::fill(modified.begin(), modified.end(), modif::nothing);
        modified[0] = modif::dataStructure;    // Fluid
        modified[1] = modif::staticVariables;  // rhoBar.
        modified[2] = modif::staticVariables;  // j.
        modified[3] = modif::staticVariables;  // Mass
        modified[4] = modif::staticVariables;  // Volume-fraction
        modified[5] = modif::staticVariables;  // Flag-status
        modified[6] = modif::nothing;          // Normal.
        modified[7] = modif::nothing;          // Interface lists
        modified[8] = modif::nothing;          // Curvature.
        modified[9] = modif::nothing;          // Outside density.
    }
private:
    Dynamics<T,Descriptor>* dynamicsTemplate;
    Array<T,2> g;
    T rhoIni;
};

template<typename T,class InsideFunction>
class AnalyticalIniVolumeFraction2D : public BoxProcessingFunctional2D {
public:
    AnalyticalIniVolumeFraction2D(InsideFunction const& insideFunction_, plint subDivision_=5)
        : insideFunction(insideFunction_),
          subDivision(subDivision_)
    { PLB_ASSERT( subDivision > 1 ); }
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D*> atomicBlocks);
    virtual AnalyticalIniVolumeFraction2D<T,InsideFunction>* clone() const {
        return new AnalyticalIniVolumeFraction2D<T,InsideFunction>(*this);
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::staticVariables;  // Volume-fraction
        modified[1] = modif::staticVariables;  // Flag-status
    }
private:
    void subDomainVolumeFraction (
            plint iX, plint iY, plint iZ, int& flag, T& volumeFraction );
private:
    InsideFunction const& insideFunction;
    plint subDivision;
};

template<typename T, class InsideFunction>
void analyticalIniVolumeFraction (
        MultiScalarField2D<T>& volumeFraction, MultiScalarField2D<int>& flagStatus,
        InsideFunction const& insideFunction, plint subDivision = 5 );


template<typename T, template<typename U> class Descriptor>
class ConstantIniVelocityFreeSurface2D : public BoxProcessingFunctional2D {
public:
    ConstantIniVelocityFreeSurface2D(Array<T,2> velocity_, T rhoIni_)
        : velocity(velocity_),
          rhoIni(rhoIni_)
    { }
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D*> atomicBlocks);
    virtual ConstantIniVelocityFreeSurface2D<T,Descriptor>* clone() const {
        return new ConstantIniVelocityFreeSurface2D<T,Descriptor>(*this);
    }
    virtual void getTypeOfModification (std::vector<modif::ModifT>& modified) const {
        std::fill(modified.begin(), modified.end(), modif::nothing);
        modified[0] = modif::staticVariables; // Fluid.
        modified[1] = modif::nothing;         // rhoBar.
        modified[2] = modif::staticVariables; // j.
        modified[3] = modif::nothing;         // Mass.
        modified[4] = modif::nothing;         // Volume-fraction.
        modified[5] = modif::nothing;         // Flag-status.
        modified[6] = modif::nothing;         // Normal.
        modified[7] = modif::nothing;         // Interface-lists.
        modified[8] = modif::nothing;         // Curvature.
        modified[9] = modif::nothing;         // Outside density.
    }
private:
    Array<T,2> velocity;
    T rhoIni;
};

template<typename T, template<typename U> class Descriptor>
class InletConstVolumeFraction2D : public BoxProcessingFunctional2D {
public:
    InletConstVolumeFraction2D(T volumeFraction_)
        : volumeFraction(volumeFraction_)
    { }
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D*> atomicBlocks);
    virtual InletConstVolumeFraction2D<T,Descriptor>* clone() const {
        return new InletConstVolumeFraction2D<T,Descriptor>(*this);
    }
    virtual void getTypeOfModification (std::vector<modif::ModifT>& modified) const {
        std::fill(modified.begin(), modified.end(), modif::nothing);
        modified[0] = modif::nothing;         // Fluid.
        modified[1] = modif::nothing;         // rhoBar.
        modified[2] = modif::nothing;         // j.
        modified[3] = modif::staticVariables; // Mass.
        modified[4] = modif::nothing;         // Volume-fraction.
        modified[5] = modif::nothing;         // Flag-status.
        modified[6] = modif::nothing;         // Normal.
        modified[7] = modif::nothing;         // Interface-lists.
        modified[8] = modif::nothing;         // Curvature.
        modified[9] = modif::nothing;         // Outside density.
    }
private:
    T volumeFraction;
};

template<typename T, template<typename U> class Descriptor>
class OutletMaximumVolumeFraction2D : public BoxProcessingFunctional2D {
public:
    OutletMaximumVolumeFraction2D(T volumeFraction_)
        : volumeFraction(volumeFraction_)
    { }
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D*> atomicBlocks);
    virtual OutletMaximumVolumeFraction2D<T,Descriptor>* clone() const {
        return new OutletMaximumVolumeFraction2D<T,Descriptor>(*this);
    }
    virtual void getTypeOfModification (std::vector<modif::ModifT>& modified) const {
        std::fill(modified.begin(), modified.end(), modif::nothing);
        modified[0] = modif::nothing;         // Fluid.
        modified[1] = modif::nothing;         // rhoBar.
        modified[2] = modif::nothing;         // j.
        modified[3] = modif::staticVariables; // Mass.
        modified[4] = modif::nothing;         // Volume-fraction.
        modified[5] = modif::nothing;         // Flag-status.
        modified[6] = modif::nothing;         // Normal.
        modified[7] = modif::nothing;         // Interface-lists.
        modified[8] = modif::nothing;         // Curvature.
        modified[9] = modif::nothing;         // Outside density.
    }
private:
    T volumeFraction;
};

template<typename T, template<typename U> class Descriptor>
class OutletMaximumVolumeFraction2_2D : public BoxProcessingFunctional2D {
public:
    OutletMaximumVolumeFraction2_2D(T volumeFraction_)
        : volumeFraction(volumeFraction_)
    { }
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D*> atomicBlocks);
    virtual OutletMaximumVolumeFraction2_2D<T,Descriptor>* clone() const {
        return new OutletMaximumVolumeFraction2_2D<T,Descriptor>(*this);
    }
    virtual void getTypeOfModification (std::vector<modif::ModifT>& modified) const {
        std::fill(modified.begin(), modified.end(), modif::nothing);
        modified[0] = modif::staticVariables; // Fluid.
        modified[1] = modif::nothing;         // rhoBar.
        modified[2] = modif::nothing;         // j.
        modified[3] = modif::staticVariables; // Mass.
        modified[4] = modif::nothing;         // Volume-fraction.
        modified[5] = modif::nothing;         // Flag-status.
        modified[6] = modif::nothing;         // Normal.
        modified[7] = modif::nothing;         // Interface-lists.
        modified[8] = modif::nothing;         // Curvature.
        modified[9] = modif::nothing;         // Outside density.
    }
private:
    T volumeFraction;
};

template<typename T, template<typename U> class Descriptor>
class NoSlipMaximumVolumeFraction2D : public BoxProcessingFunctional2D {
public:
    NoSlipMaximumVolumeFraction2D(T volumeFraction_)
        : volumeFraction(volumeFraction_)
    { }
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D*> atomicBlocks);
    virtual NoSlipMaximumVolumeFraction2D<T,Descriptor>* clone() const {
        return new NoSlipMaximumVolumeFraction2D<T,Descriptor>(*this);
    }
    virtual void getTypeOfModification (std::vector<modif::ModifT>& modified) const {
        std::fill(modified.begin(), modified.end(), modif::nothing);
        modified[0] = modif::nothing;         // Fluid.
        modified[1] = modif::staticVariables; // rhoBar.
        modified[2] = modif::nothing;         // j.
        modified[3] = modif::nothing;         // Mass.
        modified[4] = modif::nothing;         // Volume-fraction.
        modified[5] = modif::nothing;         // Flag-status.
        modified[6] = modif::nothing;         // Normal.
        modified[7] = modif::nothing;         // Interface-lists.
        modified[8] = modif::nothing;         // Curvature.
        modified[9] = modif::nothing;         // Outside density.
    }
private:
    T volumeFraction;
};

template<typename T, template<typename U> class Descriptor>
class PunchSphere2D : public BoxProcessingFunctional2D {
public:
    PunchSphere2D(Array<T,2> const& center_, T radius_, T rho0_)
        : center(center_),
          radius(radius_),
          rho0(rho0_)
    { }
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D*> atomicBlocks);
    virtual PunchSphere2D<T,Descriptor>* clone() const {
        return new PunchSphere2D<T,Descriptor>(*this);
    }
    virtual void getTypeOfModification (std::vector<modif::ModifT>& modified) const {
        std::fill(modified.begin(), modified.end(), modif::nothing);
        modified[0] = modif::staticVariables; // Fluid.
        modified[1] = modif::staticVariables; // rhoBar.
        modified[2] = modif::staticVariables; // j.
        modified[3] = modif::staticVariables; // Mass.
        modified[4] = modif::staticVariables; // Volume-fraction.
        modified[5] = modif::staticVariables; // Flag-status.
        modified[6] = modif::nothing;         // Normal.
        modified[7] = modif::nothing;         // Interface-lists.
        modified[8] = modif::staticVariables; // Curvature.
        modified[9] = modif::staticVariables; // Outside density.
    }
private:
    Array<T,2> center;
    T radius;
    T rho0;
};

template<typename T, template<typename U> class Descriptor>
class CalculateAverageSphereDensity2D : public PlainReductiveBoxProcessingFunctional2D {
public:
    CalculateAverageSphereDensity2D(Array<T,2> const& center_, T radius_)
        : center(center_),
          radius(radius_),
          averageDensityId(this->getStatistics().subscribeAverage())
    { }
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D*> atomicBlocks);
    virtual CalculateAverageSphereDensity2D<T,Descriptor>* clone() const {
        return new CalculateAverageSphereDensity2D<T,Descriptor>(*this);
    }
    T getAverageDensity() const {
        return this->getStatistics().getAverage(averageDensityId);
    }
    virtual void getTypeOfModification (std::vector<modif::ModifT>& modified) const {
        std::fill(modified.begin(), modified.end(), modif::nothing);
        modified[0] = modif::staticVariables; // Fluid.
        modified[1] = modif::staticVariables; // rhoBar.
        modified[2] = modif::staticVariables; // j.
        modified[3] = modif::staticVariables; // Mass.
        modified[4] = modif::staticVariables; // Volume-fraction.
        modified[5] = modif::staticVariables; // Flag-status.
        modified[6] = modif::nothing;         // Normal.
        modified[7] = modif::nothing;         // Interface-lists.
        modified[8] = modif::staticVariables; // Curvature.
        modified[9] = modif::staticVariables; // Outside density.
    }
private:
    Array<T,2> center;
    T radius;
    plint averageDensityId;
};


}  // namespace plb

#endif  // FREE_SURFACE_INITIALIZER_2D_H

