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

#ifndef MULTI_FREE_SURFACE_MODEL_2D_H
#define MULTI_FREE_SURFACE_MODEL_2D_H

#include <algorithm>
#include "core/globalDefs.h"
#include "atomicBlock/dataProcessingFunctional2D.h"
#include "multiBlock/defaultMultiBlockPolicy2D.h"
#include "multiPhysics/freeSurfaceModel2D.h"
#include "multiPhysics/freeSurfaceUtil2D.h"
#include "multiPhysics/freeSurfaceInitializer2D.h"
#include "dataProcessors/dataInitializerWrapper2D.h"

namespace plb {

// Couples the velocities of the two fluids. Fluid 1 sees fluid 2, but fluid 2
// does not see fluid 1.
template<typename T, template<typename U> class Descriptor>
class MultiFreeSurfaceOneWayCoupling2D : public BoxProcessingFunctional2D {
public:
    MultiFreeSurfaceOneWayCoupling2D(T interactionStrength_, T rhoDefault1_, T rhoDefault2_)
        : interactionStrength(interactionStrength_),
          rhoDefault1(rhoDefault1_),
          rhoDefault2(rhoDefault2_)
    { }
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D*> atomicBlocks);
    virtual MultiFreeSurfaceOneWayCoupling2D<T,Descriptor>* clone() const {
        return new MultiFreeSurfaceOneWayCoupling2D<T,Descriptor>(*this);
    }
    virtual void getTypeOfModification (std::vector<modif::ModifT>& modified) const {
        std::fill(modified.begin(), modified.end(), modif::nothing);
        modified[0 ] = modif::nothing;         // 1 Fluid.
        modified[1 ] = modif::nothing;         // 1 rhoBar.
        modified[2 ] = modif::staticVariables; // 1 j.
        modified[3 ] = modif::nothing;         // 1 Mass.
        modified[4 ] = modif::nothing;         // 1 Volume fraction.
        modified[5 ] = modif::nothing;         // 1 Flag-status.
        modified[6 ] = modif::nothing;         // 1 Normal.
        modified[7 ] = modif::nothing;         // 1 Interface-lists.
        modified[8 ] = modif::nothing;         // 1 Curvature.
        modified[9 ] = modif::nothing;         // 1 Outside density.

        modified[10] = modif::nothing;         // 2 Fluid.
        modified[11] = modif::nothing;         // 2 rhoBar.
        modified[12] = modif::nothing;         // 2 j.
        modified[13] = modif::nothing;         // 2 Mass.
        modified[14] = modif::nothing;         // 2 Volume fraction.
        modified[15] = modif::nothing;         // 2 Flag-status.
        modified[16] = modif::nothing;         // 2 Normal.
        modified[17] = modif::nothing;         // 2 Interface-lists.
        modified[18] = modif::nothing;         // 2 Curvature.
        modified[19] = modif::nothing;         // 2 Outside density.
    }
private:
    T interactionStrength;
    T rhoDefault1, rhoDefault2;
};

template<typename T, template<typename U> class Descriptor>
class MultiFreeSurfaceVelocityContinuityCoupling2D : public BoxProcessingFunctional2D {
public:
    MultiFreeSurfaceVelocityContinuityCoupling2D(T rhoDefault1_, T rhoDefault2_)
        : rhoDefault1(rhoDefault1_),
          rhoDefault2(rhoDefault2_)
    { }
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D*> atomicBlocks);
    virtual MultiFreeSurfaceVelocityContinuityCoupling2D<T,Descriptor>* clone() const {
        return new MultiFreeSurfaceVelocityContinuityCoupling2D<T,Descriptor>(*this);
    }
    virtual void getTypeOfModification (std::vector<modif::ModifT>& modified) const {
        std::fill(modified.begin(), modified.end(), modif::nothing);
        modified[0 ] = modif::nothing;         // 1 Fluid.
        modified[1 ] = modif::staticVariables; // 1 rhoBar.
        modified[2 ] = modif::staticVariables; // 1 j.
        modified[3 ] = modif::nothing;         // 1 Mass.
        modified[4 ] = modif::staticVariables; // 1 Volume fraction.
        modified[5 ] = modif::nothing;         // 1 Flag-status.
        modified[6 ] = modif::nothing;         // 1 Normal.
        modified[7 ] = modif::nothing;         // 1 Interface-lists.
        modified[8 ] = modif::nothing;         // 1 Curvature.
        modified[9 ] = modif::nothing;         // 1 Outside density.

        modified[10] = modif::nothing;         // 2 Fluid.
        modified[11] = modif::staticVariables; // 2 rhoBar.
        modified[12] = modif::staticVariables; // 2 j.
        modified[13] = modif::nothing;         // 2 Mass.
        modified[14] = modif::staticVariables; // 2 Volume fraction.
        modified[15] = modif::nothing;         // 2 Flag-status.
        modified[16] = modif::nothing;         // 2 Normal.
        modified[17] = modif::nothing;         // 2 Interface-lists.
        modified[18] = modif::nothing;         // 2 Curvature.
        modified[19] = modif::nothing;         // 2 Outside density.
    }
private:
    T rhoDefault1, rhoDefault2;
};

template<typename T, template<typename U> class Descriptor>
class MultiFreeSurfaceRepellingForceCoupling2D : public BoxProcessingFunctional2D {
public:
    MultiFreeSurfaceRepellingForceCoupling2D(T interactionStrength_, T rhoDefault1_, T rhoDefault2_)
        : interactionStrength(interactionStrength_),
          rhoDefault1(rhoDefault1_),
          rhoDefault2(rhoDefault2_)
    { }
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D*> atomicBlocks);
    virtual MultiFreeSurfaceRepellingForceCoupling2D<T,Descriptor>* clone() const {
        return new MultiFreeSurfaceRepellingForceCoupling2D<T,Descriptor>(*this);
    }
    virtual void getTypeOfModification (std::vector<modif::ModifT>& modified) const {
        std::fill(modified.begin(), modified.end(), modif::nothing);
        modified[0 ] = modif::nothing;         // 1 Fluid.
        modified[1 ] = modif::staticVariables; // 1 rhoBar.
        modified[2 ] = modif::staticVariables; // 1 j.
        modified[3 ] = modif::nothing;         // 1 Mass.
        modified[4 ] = modif::staticVariables; // 1 Volume fraction.
        modified[5 ] = modif::nothing;         // 1 Flag-status.
        modified[6 ] = modif::nothing;         // 1 Normal.
        modified[7 ] = modif::nothing;         // 1 Interface-lists.
        modified[8 ] = modif::nothing;         // 1 Curvature.
        modified[9 ] = modif::nothing;         // 1 Outside density.

        modified[10] = modif::nothing;         // 2 Fluid.
        modified[11] = modif::staticVariables; // 2 rhoBar.
        modified[12] = modif::staticVariables; // 2 j.
        modified[13] = modif::nothing;         // 2 Mass.
        modified[14] = modif::staticVariables; // 2 Volume fraction.
        modified[15] = modif::nothing;         // 2 Flag-status.
        modified[16] = modif::nothing;         // 2 Normal.
        modified[17] = modif::nothing;         // 2 Interface-lists.
        modified[18] = modif::nothing;         // 2 Curvature.
        modified[19] = modif::nothing;         // 2 Outside density.
    }
private:
    T deltaFunction(T r, T h);
private:
    T interactionStrength;
    T rhoDefault1, rhoDefault2;
};

template<typename T, template<typename U> class Descriptor>
class MultiFreeSurfaceComplexCoupling2D : public BoxProcessingFunctional2D {
public:
    MultiFreeSurfaceComplexCoupling2D(T interactionStrength_, T rhoDefault1_, T rhoDefault2_)
        : interactionStrength(interactionStrength_),
          rhoDefault1(rhoDefault1_),
          rhoDefault2(rhoDefault2_)
    { }
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D*> atomicBlocks);
    virtual MultiFreeSurfaceComplexCoupling2D<T,Descriptor>* clone() const {
        return new MultiFreeSurfaceComplexCoupling2D<T,Descriptor>(*this);
    }
    virtual void getTypeOfModification (std::vector<modif::ModifT>& modified) const {
        std::fill(modified.begin(), modified.end(), modif::nothing);
        modified[0 ] = modif::nothing;         // 1 Fluid.
        modified[1 ] = modif::staticVariables; // 1 rhoBar.
        modified[2 ] = modif::staticVariables; // 1 j.
        modified[3 ] = modif::nothing;         // 1 Mass.
        modified[4 ] = modif::staticVariables; // 1 Volume fraction.
        modified[5 ] = modif::nothing;         // 1 Flag-status.
        modified[6 ] = modif::nothing;         // 1 Normal.
        modified[7 ] = modif::nothing;         // 1 Interface-lists.
        modified[8 ] = modif::nothing;         // 1 Curvature.
        modified[9 ] = modif::nothing;         // 1 Outside density.

        modified[10] = modif::nothing;         // 2 Fluid.
        modified[11] = modif::staticVariables; // 2 rhoBar.
        modified[12] = modif::staticVariables; // 2 j.
        modified[13] = modif::nothing;         // 2 Mass.
        modified[14] = modif::staticVariables; // 2 Volume fraction.
        modified[15] = modif::nothing;         // 2 Flag-status.
        modified[16] = modif::nothing;         // 2 Normal.
        modified[17] = modif::nothing;         // 2 Interface-lists.
        modified[18] = modif::nothing;         // 2 Curvature.
        modified[19] = modif::nothing;         // 2 Outside density.
    }
private:
    T deltaFunction(T r, T h);
private:
    T interactionStrength;
    T rhoDefault1, rhoDefault2;
};

template<typename T, template<typename U> class Descriptor>
struct MultiFreeSurfaceFields2D {
    //static const int envelopeWidth = 2;
    static const int envelopeWidth = 3;
    static const int smallEnvelopeWidth = 1;

    MultiFreeSurfaceFields2D(SparseBlockStructure2D const& blockStructure,
                             Dynamics<T,Descriptor> *dynamics1_, Dynamics<T,Descriptor> *dynamics2_,
                             T rhoDefault1_, T rhoDefault2_, T surfaceTension1_, T surfaceTension2_,
                             T contactAngle1_, T contactAngle2_, Array<T,2> force1_, Array<T,2> force2_,
                             T interactionStrength_,
                             ThreadAttribution* threadAttribution = defaultMultiBlockPolicy2D().getThreadAttribution() )
        : dynamics1(dynamics1_),
          dynamics2(dynamics2_),
          rhoDefault1(rhoDefault1_),
          rhoDefault2(rhoDefault2_),
          surfaceTension1(surfaceTension1_),
          surfaceTension2(surfaceTension2_),
          contactAngle1(contactAngle1_),
          contactAngle2(contactAngle2_),
          force1(force1_),
          force2(force2_),
          lattice1 (
                MultiBlockManagement2D (
                        blockStructure, threadAttribution->clone(), smallEnvelopeWidth ),
                defaultMultiBlockPolicy2D().getBlockCommunicator(),
                defaultMultiBlockPolicy2D().getCombinedStatistics(),
                defaultMultiBlockPolicy2D().getMultiCellAccess<T,Descriptor>(), dynamics1->clone() ),
          lattice2 (
                MultiBlockManagement2D (
                        blockStructure, threadAttribution->clone(), smallEnvelopeWidth ),
                defaultMultiBlockPolicy2D().getBlockCommunicator(),
                defaultMultiBlockPolicy2D().getCombinedStatistics(),
                defaultMultiBlockPolicy2D().getMultiCellAccess<T,Descriptor>(), dynamics2->clone() ),
          helperLists1(lattice1),
          helperLists2(lattice2),
          mass1(lattice1),
          mass2(lattice2),
          flag1 (
                MultiBlockManagement2D (
                        blockStructure, threadAttribution->clone(), envelopeWidth ),
                defaultMultiBlockPolicy2D().getBlockCommunicator(),
                defaultMultiBlockPolicy2D().getCombinedStatistics(),
                defaultMultiBlockPolicy2D().getMultiScalarAccess<int>() ),
          flag2 (
                MultiBlockManagement2D (
                        blockStructure, threadAttribution->clone(), envelopeWidth ),
                defaultMultiBlockPolicy2D().getBlockCommunicator(),
                defaultMultiBlockPolicy2D().getCombinedStatistics(),
                defaultMultiBlockPolicy2D().getMultiScalarAccess<int>() ),
          volumeFraction1((MultiBlock2D&) flag1),
          volumeFraction2((MultiBlock2D&) flag2),
          curvature1((MultiBlock2D&) flag1),
          curvature2((MultiBlock2D&) flag2),
          outsideDensity1((MultiBlock2D&) flag1),
          outsideDensity2((MultiBlock2D&) flag2),
          rhoBar1(lattice1),
          rhoBar2(lattice2),
          j1(lattice1),
          j2(lattice2),
          normal1((MultiBlock2D&) flag1),
          normal2((MultiBlock2D&) flag2),
          interactionStrength(interactionStrength_)
    {
        delete threadAttribution;
        Precision precision;
        if (sizeof(T) == sizeof(float))
            precision = FLT;
        else if (sizeof(T) == sizeof(double))
            precision = DBL;
        else if (sizeof(T) == sizeof(long double))
            precision = LDBL;
        else
            PLB_ASSERT(false);

        T eps = getEpsilon<T>(precision);
        // The contact angles take values between 0 and 180 degrees. If they
        // are negative, this means that contact angle effects will not be
        // modeled.
        PLB_ASSERT(contactAngle1 < (T) 180.0 || fabs(contactAngle1 - (T) 180.0) <= eps);
        PLB_ASSERT(contactAngle2 < (T) 180.0 || fabs(contactAngle2 - (T) 180.0) <= eps);

        if (fabs(surfaceTension1) <= eps) {
            useSurfaceTension1 = 0;
        } else {
            useSurfaceTension1 = 1;
        }

        if (fabs(surfaceTension2) <= eps) {
            useSurfaceTension2 = 0;
        } else {
            useSurfaceTension2 = 1;
        }

        twoPhaseArgs1 = aggregateFreeSurfaceParams2D(lattice1, rhoBar1, j1, mass1, volumeFraction1,
                    flag1, normal1, helperLists1, curvature1, outsideDensity1);

        twoPhaseArgs2 = aggregateFreeSurfaceParams2D(lattice2, rhoBar2, j2, mass2, volumeFraction2,
                    flag2, normal2, helperLists2, curvature2, outsideDensity2);

        multiFreeSurfaceArgs = twoPhaseArgs1;
        multiFreeSurfaceArgs.insert(multiFreeSurfaceArgs.end(), twoPhaseArgs2.begin(), twoPhaseArgs2.end());

        initializeInterfaceLists2D<T,Descriptor>(helperLists1);
        initializeInterfaceLists2D<T,Descriptor>(helperLists2);
        lattice1.periodicity().toggleAll(true);
        lattice2.periodicity().toggleAll(true);
        mass1.periodicity().toggleAll(true);
        mass2.periodicity().toggleAll(true);
        flag1.periodicity().toggleAll(true);
        flag2.periodicity().toggleAll(true);
        volumeFraction1.periodicity().toggleAll(true);
        volumeFraction2.periodicity().toggleAll(true);
        curvature1.periodicity().toggleAll(true);
        curvature2.periodicity().toggleAll(true);
        outsideDensity1.periodicity().toggleAll(true);
        outsideDensity2.periodicity().toggleAll(true);

        rhoBar1.periodicity().toggleAll(true);
        rhoBar2.periodicity().toggleAll(true);
        j1.periodicity().toggleAll(true);
        j2.periodicity().toggleAll(true);
        normal1.periodicity().toggleAll(true);
        normal2.periodicity().toggleAll(true);
        setToConstant(flag1, flag1.getBoundingBox(), (int) twoPhaseFlag::empty);
        setToConstant(flag2, flag2.getBoundingBox(), (int) twoPhaseFlag::empty);
        setToConstant(outsideDensity1, outsideDensity1.getBoundingBox(), rhoDefault1);
        setToConstant(outsideDensity2, outsideDensity2.getBoundingBox(), rhoDefault2);
        rhoBarJparam1.push_back(&lattice1);
        rhoBarJparam1.push_back(&rhoBar1);
        rhoBarJparam1.push_back(&j1);
        rhoBarJparam2.push_back(&lattice2);
        rhoBarJparam2.push_back(&rhoBar2);
        rhoBarJparam2.push_back(&j2);

        lattice1.internalStatSubscription().subscribeSum();     // Total mass.
        lattice1.internalStatSubscription().subscribeSum();     // Lost mass.
        lattice1.internalStatSubscription().subscribeIntSum();  // Num interface cells.

        lattice2.internalStatSubscription().subscribeSum();     // Total mass.
        lattice2.internalStatSubscription().subscribeSum();     // Lost mass.
        lattice2.internalStatSubscription().subscribeIntSum();  // Num interface cells.

        freeSurfaceDataProcessors(rhoDefault1, rhoDefault2, force1, force2, *dynamics1, *dynamics2);
    }

    MultiFreeSurfaceFields2D(MultiFreeSurfaceFields2D<T,Descriptor> const& rhs)
        : dynamics1(rhs.dynamics1->clone()),
          dynamics2(rhs.dynamics2->clone()),
          rhoDefault1(rhs.rhoDefault1),
          rhoDefault2(rhs.rhoDefault2),
          surfaceTension1(rhs.surfaceTension1),
          surfaceTension2(rhs.surfaceTension2),
          contactAngle1(rhs.contactAngle1),
          contactAngle2(rhs.contactAngle2),
          useSurfaceTension1(rhs.useSurfaceTension1),
          useSurfaceTension2(rhs.useSurfaceTension2),
          force1(rhs.force1),
          force2(rhs.force2),
          lattice1(rhs.lattice1),
          lattice2(rhs.lattice2),
          helperLists1(rhs.helperLists1),
          helperLists2(rhs.helperLists2),
          mass1(rhs.mass1),
          mass2(rhs.mass2),
          flag1(rhs.flag1),
          flag2(rhs.flag2),
          volumeFraction1(rhs.volumeFraction1),
          volumeFraction2(rhs.volumeFraction2),
          curvature1(rhs.curvature1),
          curvature2(rhs.curvature2),
          outsideDensity1(rhs.outsideDensity1),
          outsideDensity2(rhs.outsideDensity2),
          rhoBar1(rhs.rhoBar1),
          rhoBar2(rhs.rhoBar2),
          j1(rhs.j1),
          j2(rhs.j2),
          normal1(rhs.normal1),
          normal2(rhs.normal2),
          rhoBarJparam1(rhs.rhoBarJparam1),
          rhoBarJparam2(rhs.rhoBarJparam2),
          twoPhaseArgs1(rhs.twoPhaseArgs1),
          twoPhaseArgs2(rhs.twoPhaseArgs2),
          multiFreeSurfaceArgs(rhs.multiFreeSurfaceArgs),
          interactionStrength(rhs.interactionStrength)
    { }

    void swap(MultiFreeSurfaceFields2D<T,Descriptor>& rhs)
    {
        std::swap(dynamics1, rhs.dynamics1);
        std::swap(dynamics2, rhs.dynamics2);
        std::swap(rhoDefault1, rhs.rhoDefault1);
        std::swap(rhoDefault2, rhs.rhoDefault2);
        std::swap(surfaceTension1, rhs.surfaceTension1);
        std::swap(surfaceTension2, rhs.surfaceTension2);
        std::swap(contactAngle1, rhs.contactAngle1);
        std::swap(contactAngle2, rhs.contactAngle2);
        std::swap(useSurfaceTension1, rhs.useSurfaceTension1);
        std::swap(useSurfaceTension2, rhs.useSurfaceTension2);
        std::swap(force1, rhs.force1);
        std::swap(force2, rhs.force2);
        std::swap(lattice1, rhs.lattice1);
        std::swap(lattice2, rhs.lattice2);
        std::swap(helperLists1, rhs.helperLists1);
        std::swap(helperLists2, rhs.helperLists2);
        std::swap(mass1, rhs.mass1);
        std::swap(mass2, rhs.mass2);
        std::swap(flag1, rhs.flag1);
        std::swap(flag2, rhs.flag2);
        std::swap(volumeFraction1, rhs.volumeFraction1);
        std::swap(volumeFraction2, rhs.volumeFraction2);
        std::swap(curvature1, rhs.curvature1);
        std::swap(curvature2, rhs.curvature2);
        std::swap(outsideDensity1, rhs.outsideDensity1);
        std::swap(outsideDensity2, rhs.outsideDensity2);
        std::swap(rhoBar1, rhs.rhoBar1);
        std::swap(rhoBar2, rhs.rhoBar2);
        std::swap(j1, rhs.j1);
        std::swap(j2, rhs.j2);
        std::swap(normal1, rhs.normal1);
        std::swap(normal2, rhs.normal2);
        std::swap(rhoBarJparam1, rhs.rhoBarJparam1);
        std::swap(rhoBarJparam2, rhs.rhoBarJparam2);
        std::swap(twoPhaseArgs1, rhs.twoPhaseArgs1);
        std::swap(twoPhaseArgs2, rhs.twoPhaseArgs2);
        std::swap(multiFreeSurfaceArgs, rhs.multiFreeSurfaceArgs);
        std::swap(interactionStrength, rhs.interactionStrength);
    }

    MultiFreeSurfaceFields2D<T,Descriptor>& operator=(MultiFreeSurfaceFields2D<T,Descriptor> const& rhs)
    {
        MultiFreeSurfaceFields2D<T,Descriptor>(rhs).swap(*this);
        return *this;
    }

    MultiFreeSurfaceFields2D<T,Descriptor>* clone() const
    {
        return new MultiFreeSurfaceFields2D<T,Descriptor>(*this);
    }

    ~MultiFreeSurfaceFields2D() {
        delete dynamics1;
        delete dynamics2;
    }

    void periodicityToggle(plint direction, bool periodic)
    {
        PLB_ASSERT(direction == 0 || direction == 1 || direction == 2);

        lattice1.periodicity().toggle(direction, periodic);
        lattice2.periodicity().toggle(direction, periodic);
        mass1.periodicity().toggle(direction, periodic);
        mass2.periodicity().toggle(direction, periodic);
        flag1.periodicity().toggle(direction, periodic);
        flag2.periodicity().toggle(direction, periodic);
        volumeFraction1.periodicity().toggle(direction, periodic);
        volumeFraction2.periodicity().toggle(direction, periodic);
        curvature1.periodicity().toggle(direction, periodic);
        curvature2.periodicity().toggle(direction, periodic);
        outsideDensity1.periodicity().toggle(direction, periodic);
        outsideDensity2.periodicity().toggle(direction, periodic);
        rhoBar1.periodicity().toggle(direction, periodic);
        rhoBar2.periodicity().toggle(direction, periodic);
        j1.periodicity().toggle(direction, periodic);
        j2.periodicity().toggle(direction, periodic);
        normal1.periodicity().toggle(direction, periodic);
        normal2.periodicity().toggle(direction, periodic);
    }

    void periodicityToggleAll(bool periodic)
    {
        lattice1.periodicity().toggleAll(periodic);
        lattice2.periodicity().toggleAll(periodic);
        mass1.periodicity().toggleAll(periodic);
        mass2.periodicity().toggleAll(periodic);
        flag1.periodicity().toggleAll(periodic);
        flag2.periodicity().toggleAll(periodic);
        volumeFraction1.periodicity().toggleAll(periodic);
        volumeFraction2.periodicity().toggleAll(periodic);
        curvature1.periodicity().toggleAll(periodic);
        curvature2.periodicity().toggleAll(periodic);
        outsideDensity1.periodicity().toggleAll(periodic);
        outsideDensity2.periodicity().toggleAll(periodic);
        rhoBar1.periodicity().toggleAll(periodic);
        rhoBar2.periodicity().toggleAll(periodic);
        j1.periodicity().toggleAll(periodic);
        j2.periodicity().toggleAll(periodic);
        normal1.periodicity().toggleAll(periodic);
        normal2.periodicity().toggleAll(periodic);
    }

    void defaultInitialize() {
        applyProcessingFunctional (
           new DefaultInitializeFreeSurface2D<T,Descriptor>(dynamics1->clone(), force1, rhoDefault1),
                   lattice1.getBoundingBox(), twoPhaseArgs1 );

        applyProcessingFunctional (
           new DefaultInitializeFreeSurface2D<T,Descriptor>(dynamics2->clone(), force2, rhoDefault2),
                   lattice2.getBoundingBox(), twoPhaseArgs2 );
    }

    void partiallyDefaultInitialize() {
        applyProcessingFunctional (
           new PartiallyDefaultInitializeFreeSurface2D<T,Descriptor>(dynamics1->clone(), force1, rhoDefault1),
                   lattice1.getBoundingBox(), twoPhaseArgs1 );

        applyProcessingFunctional (
           new PartiallyDefaultInitializeFreeSurface2D<T,Descriptor>(dynamics2->clone(), force2, rhoDefault2),
                   lattice2.getBoundingBox(), twoPhaseArgs2 );
    }

    void freeSurfaceDataProcessors(T rhoDefault1, T rhoDefault2, Array<T,2> force1, Array<T,2> force2,
            Dynamics<T,Descriptor>& dynamics1, Dynamics<T,Descriptor>& dynamics2)
    {
        MultiBlock2D& actor1 = *twoPhaseArgs2[0];
        MultiBlock2D& actor2 = *twoPhaseArgs2[0];

        plint pl; // Processor level.

        /***** Initial level ******/
        pl = 0;

        integrateProcessingFunctional (
                new ExternalRhoJcollideAndStream2D<T,Descriptor>,
                lattice1.getBoundingBox(), rhoBarJparam1, pl );
        integrateProcessingFunctional (
                new ExternalRhoJcollideAndStream2D<T,Descriptor>,
                lattice2.getBoundingBox(), rhoBarJparam2, pl );

        integrateProcessingFunctional (
                new TwoPhaseComputeNormals2D<T,Descriptor>,
                lattice1.getBoundingBox(), actor1, twoPhaseArgs1, pl );
        integrateProcessingFunctional (
                new TwoPhaseComputeNormals2D<T,Descriptor>,
                lattice2.getBoundingBox(), actor2, twoPhaseArgs2, pl );

        /***** New level ******/
        pl++;

        if (useSurfaceTension1) {
            integrateProcessingFunctional (
                    new TwoPhaseComputeCurvature2D<T,Descriptor>(contactAngle1, lattice1.getBoundingBox()),
                    lattice1.getBoundingBox(), actor1, twoPhaseArgs1, pl );
        }
        if (useSurfaceTension2) {
            integrateProcessingFunctional (
                    new TwoPhaseComputeCurvature2D<T,Descriptor>(contactAngle2, lattice2.getBoundingBox()),
                    lattice2.getBoundingBox(), actor2, twoPhaseArgs2, pl );
        }

        integrateProcessingFunctional (
            new FreeSurfaceMassChange2D<T,Descriptor>, lattice1.getBoundingBox(),
            actor1, twoPhaseArgs1, pl );
        integrateProcessingFunctional (
            new FreeSurfaceMassChange2D<T,Descriptor>, lattice2.getBoundingBox(),
            actor2, twoPhaseArgs2, pl );
       
        integrateProcessingFunctional (
            new FreeSurfaceCompletion2D<T,Descriptor>,
            lattice1.getBoundingBox(), actor1, twoPhaseArgs1, pl );
        integrateProcessingFunctional (
            new FreeSurfaceCompletion2D<T,Descriptor>,
            lattice2.getBoundingBox(), actor2, twoPhaseArgs2, pl );
                                    
        integrateProcessingFunctional (
            new FreeSurfaceMacroscopic2D<T,Descriptor>(rhoDefault1),
            lattice1.getBoundingBox(), actor1, twoPhaseArgs1, pl );
        integrateProcessingFunctional (
            new FreeSurfaceMacroscopic2D<T,Descriptor>(rhoDefault2),
            lattice2.getBoundingBox(), actor2, twoPhaseArgs2, pl );

        /***** New level ******/
        pl++;

        Precision precision;
        if (sizeof(T) == sizeof(float))
            precision = FLT;
        else if (sizeof(T) == sizeof(double))
            precision = DBL;
        else if (sizeof(T) == sizeof(long double))
            precision = LDBL;
        else
            PLB_ASSERT(false);

        T eps = getEpsilon<T>(precision);

        int useRepellingForceCoupling = 1;
        if (fabs(interactionStrength) <= eps) {
            useRepellingForceCoupling = 0;
        }

        if (useRepellingForceCoupling) {
            //integrateProcessingFunctional (
            //    new MultiFreeSurfaceRepellingForceCoupling2D<T,Descriptor>(interactionStrength, rhoDefault1, rhoDefault2),
            //    lattice2.getBoundingBox(), actor2, multiFreeSurfaceArgs, pl );
            integrateProcessingFunctional (
                new MultiFreeSurfaceComplexCoupling2D<T,Descriptor>(interactionStrength, rhoDefault1, rhoDefault2),
                lattice2.getBoundingBox(), actor2, multiFreeSurfaceArgs, pl );
        } else {
            integrateProcessingFunctional (
                new MultiFreeSurfaceVelocityContinuityCoupling2D<T,Descriptor>(rhoDefault1, rhoDefault2),
                lattice2.getBoundingBox(), actor2, multiFreeSurfaceArgs, pl );
        }

        //integrateProcessingFunctional (
        //    new MultiFreeSurfaceOneWayCoupling2D<T,Descriptor>(interactionStrength,rhoDefault1, rhoDefault2),
        //    lattice2.getBoundingBox(), actor2, multiFreeSurfaceArgs, pl );

        /***** New level ******/
        if (useSurfaceTension1 || useSurfaceTension2)
            pl++;

        if (useSurfaceTension1) {
            integrateProcessingFunctional (
                new TwoPhaseAddSurfaceTension2D<T,Descriptor>(surfaceTension1),
                lattice1.getBoundingBox(), actor1, twoPhaseArgs1, pl );
        }
        if (useSurfaceTension2) {
            integrateProcessingFunctional (
                new TwoPhaseAddSurfaceTension2D<T,Descriptor>(surfaceTension2),
                lattice2.getBoundingBox(), actor2, twoPhaseArgs2, pl );
        }

        /***** New level ******/ // Maybe no new level is necessary here ...
        pl++;

        integrateProcessingFunctional (
            new FreeSurfaceComputeInterfaceLists2D<T,Descriptor>(),
            lattice1.getBoundingBox(), actor1, twoPhaseArgs1, pl );
        integrateProcessingFunctional (
            new FreeSurfaceComputeInterfaceLists2D<T,Descriptor>(),
            lattice2.getBoundingBox(), actor2, twoPhaseArgs2, pl );

        integrateProcessingFunctional (
            new FreeSurfaceIniInterfaceToAnyNodes2D<T,Descriptor>(rhoDefault1),
            lattice1.getBoundingBox(), actor1, twoPhaseArgs1, pl );
        integrateProcessingFunctional (
            new FreeSurfaceIniInterfaceToAnyNodes2D<T,Descriptor>(rhoDefault2),
            lattice2.getBoundingBox(), actor2, twoPhaseArgs2, pl );
            
        integrateProcessingFunctional (
            new FreeSurfaceIniEmptyToInterfaceNodes2D<T,Descriptor>(dynamics1.clone(), force1),
                                    lattice1.getBoundingBox(),
                                    actor1, twoPhaseArgs1, pl ); 
        integrateProcessingFunctional (
            new FreeSurfaceIniEmptyToInterfaceNodes2D<T,Descriptor>(dynamics2.clone(), force2),
                                    lattice2.getBoundingBox(),
                                    actor2, twoPhaseArgs2, pl ); 

        /***** New level ******/
        pl++;

        integrateProcessingFunctional (
            new FreeSurfaceRemoveFalseInterfaceCells2D<T,Descriptor>(rhoDefault1),
            lattice1.getBoundingBox(), actor1, twoPhaseArgs1, pl);
        integrateProcessingFunctional (
            new FreeSurfaceRemoveFalseInterfaceCells2D<T,Descriptor>(rhoDefault2),
            lattice2.getBoundingBox(), actor2, twoPhaseArgs2, pl);

        /***** New level ******/
        pl++;

        integrateProcessingFunctional (
            new FreeSurfaceEqualMassExcessReDistribution2D<T,Descriptor>(),
            lattice1.getBoundingBox(), actor1, twoPhaseArgs1, pl );
        integrateProcessingFunctional (
            new FreeSurfaceEqualMassExcessReDistribution2D<T,Descriptor>(),
            lattice2.getBoundingBox(), actor2, twoPhaseArgs2, pl );

        integrateProcessingFunctional (
            new TwoPhaseComputeStatistics2D<T,Descriptor>,
            lattice1.getBoundingBox(), actor1, twoPhaseArgs1, pl );
        integrateProcessingFunctional (
            new TwoPhaseComputeStatistics2D<T,Descriptor>,
            lattice2.getBoundingBox(), actor2, twoPhaseArgs2, pl );
    }

    Dynamics<T,Descriptor> *dynamics1, *dynamics2;
    T rhoDefault1, rhoDefault2;
    T surfaceTension1, surfaceTension2;
    T contactAngle1, contactAngle2;
    int useSurfaceTension1, useSurfaceTension2;
    Array<T,2> force1, force2;
    MultiBlockLattice2D<T, Descriptor> lattice1, lattice2;
    MultiContainerBlock2D helperLists1, helperLists2;
    MultiScalarField2D<T> mass1, mass2;
    MultiScalarField2D<int> flag1, flag2;
    MultiScalarField2D<T> volumeFraction1, volumeFraction2;
    MultiScalarField2D<T> curvature1, curvature2;
    MultiScalarField2D<T> outsideDensity1, outsideDensity2;
    MultiScalarField2D<T> rhoBar1, rhoBar2;
    MultiTensorField2D<T,2> j1, j2;
    MultiTensorField2D<T,2> normal1, normal2;
    T interactionStrength;
    std::vector<MultiBlock2D*> rhoBarJparam1, rhoBarJparam2;
    std::vector<MultiBlock2D*> twoPhaseArgs1, twoPhaseArgs2, multiFreeSurfaceArgs;
};


}  // namespace plb

#endif  // MULTI_FREE_SURFACE_MODEL_2D_H

