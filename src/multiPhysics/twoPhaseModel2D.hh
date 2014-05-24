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

#ifndef TWO_PHASE_MODEL_2D_HH
#define TWO_PHASE_MODEL_2D_HH

#include "core/globalDefs.h"
#include "multiPhysics/twoPhaseModel2D.h"
#include "core/block2D.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include "atomicBlock/dataProcessor2D.h"
#include "atomicBlock/blockLattice2D.h"
#include "atomicBlock/atomicContainerBlock2D.h"
#include "multiPhysics/freeSurfaceModel2D.h"
#include "multiPhysics/freeSurfaceTemplates.h"

#include <limits>

namespace plb
{

template<typename T, template<typename U> class Descriptor>
struct twoPhaseTemplates
{

    static void massExchangeFluidCell1 (
        TwoPhaseProcessorParam2D<T,Descriptor>& param,
        plint iX, plint iY )
    {
        typedef Descriptor<T> D;
        using namespace twoPhaseFlag;
        // Calculate mass at time t+1 --> eq 6 Thurey's paper.
        for ( plint iPop=0; iPop < D::q; ++iPop )
        {
            plint nextX = iX + D::c[iPop][0];
            plint nextY = iY + D::c[iPop][1];

            int nextFlag = param.flag ( nextX,nextY );
            if ( nextFlag==fluid || nextFlag==interface )
            {
                // In Thuerey's paper, the mass balance is computed locally on one cell, but
                // N. Thuerey uses outgoing populations. Palabos works with incoming populations
                // and uses the relation f_out_i(x,t) = f_in_opp(i)(x+c_i,t+1).
                plint opp = indexTemplates::opposite<D> ( iPop );
                param.mass ( iX,iY ) += param.cell ( iX,iY ) [opp] - param.cell ( nextX,nextY ) [iPop];
            }
        }
    }

    static void massExchangeFluidCell2 (
        TwoPhaseProcessorParam2D<T,Descriptor>& param,
        plint iX, plint iY )
    {
        typedef Descriptor<T> D;
        using namespace twoPhaseFlag;
        // Calculate mass at time t+1 --> eq 6 Thurey's paper.
        for ( plint iPop=0; iPop < D::q; ++iPop )
        {
            plint nextX = iX + D::c[iPop][0];
            plint nextY = iY + D::c[iPop][1];

            int nextFlag = param.flag ( nextX,nextY );
            if ( isEmpty ( nextFlag ) || nextFlag==interface )
            {
                // In Thuerey's paper, the mass balance is computed locally on one cell, but
                // N. Thuerey uses outgoing populations. Palabos works with incoming populations
                // and uses the relation f_out_i(x,t) = f_in_opp(i)(x+c_i,t+1).
                plint opp = indexTemplates::opposite<D> ( iPop );
                param.mass2 ( iX,iY ) += param.cell2 ( iX,iY ) [opp] - param.cell2 ( nextX,nextY ) [iPop];
            }
        }
    }
};


template<typename T, template<typename U> class Descriptor>
void twoPhasePunchSphere ( std::vector<MultiBlock2D*> const& twoPhaseArgs, Array<T,2> const& center, T radius,
                           T rhoEmpty, T referenceDensity, T densityRatio, Dynamics<T,Descriptor>& dynamics,
                           TwoPhaseModel model, Box2D domain )
{
    applyProcessingFunctional (
        new TwoPhasePunchSphere2D<T,Descriptor> (
            center, radius, rhoEmpty, referenceDensity, densityRatio, dynamics.clone(), model ),
        domain, twoPhaseArgs );
}

template<typename T, template<typename U> class Descriptor>
void twoPhasePunchRectangle ( std::vector<MultiBlock2D*> const& twoPhaseArgs, Box2D rectangle,
                              T rhoEmpty, T referenceDensity, T densityRatio, Dynamics<T,Descriptor>& dynamics,
                              TwoPhaseModel model, Box2D domain )
{
    applyProcessingFunctional (
        new TwoPhasePunchRectangle2D<T,Descriptor> (
            rectangle, rhoEmpty, referenceDensity, densityRatio, dynamics.clone(), model ),
        domain, twoPhaseArgs );
}

/* *************** Class TwoPhasePunchSphere2D ******************************************* */

template< typename T,template<typename U> class Descriptor>
void TwoPhasePunchSphere2D<T,Descriptor>
::processGenericBlocks ( Box2D domain, std::vector<AtomicBlock2D*> atomicBlocks )
{
    using namespace twoPhaseFlag;
    TwoPhaseProcessorParam2D<T,Descriptor> param ( atomicBlocks );
    typedef Descriptor<T> D;

    Dot2D offset = param.absOffset();
    Array<T,2> localCenter ( center-Array<T,2> ( offset.x,offset.y ) );

    for ( plint iX=domain.x0; iX<=domain.x1; ++iX )
    {
        for ( plint iY=domain.y0; iY<=domain.y1; ++iY )
        {
            int fl = param.flag ( iX,iY );
            if ( normSqr ( Array<T,2> ( iX,iY )-localCenter ) <radius*radius &&
                    ( fl==interface || fl==fluid ) )
            {
                bool isInterface=false;
                if ( fl==interface || fl==fluid )
                {
                    for ( plint iPop=0; iPop<D::q; ++iPop )
                    {
                        plint nextX = iX+D::c[iPop][0];
                        plint nextY = iY+D::c[iPop][1];
                        if ( normSqr ( ( Array<T,2> ( nextX,nextY )-localCenter ) ) >=radius*radius )
                        {
                            isInterface=true;
                        }
                    }
                }
                T density = param.getDensity ( iX,iY );
                Array<T,2> momentum ( param.getMomentum ( iX,iY ) );
                if ( model==constRho )
                {
                    momentum *= referenceDensity/density;
                    density = referenceDensity;
                }
                T density2 ( density );
                Array<T,2> momentum2 ( momentum );
                if ( model==dynamic )
                {
                    T deltaP = density-referenceDensity;
                    if ( model!=freeSurface )
                    {
                        density2 = referenceDensity + deltaP/densityRatio;
                        momentum2 *= density2/density;
                    }
                }
                T volumeFraction = 1.0;
                if ( isInterface )
                {
                    param.flag ( iX,iY ) = interface;
                    param.volumeFraction ( iX,iY ) = volumeFraction;
                    param.mass ( iX,iY ) = volumeFraction*density;
                    param.setDensity ( iX,iY, density );
                    param.outsideDensity ( iX,iY ) = referenceDensity;

                    if ( model!=freeSurface )
                    {
                        param.mass2 ( iX,iY ) = ( 1.-volumeFraction ) *density2;
                    }
                }
                else
                {
                    param.flag ( iX,iY ) = empty;
                    param.attributeDynamics ( iX,iY, new NoDynamics<T,Descriptor>() );
                    param.mass ( iX,iY ) = T();
                    param.volumeFraction ( iX,iY ) = T();
                    param.setDensity ( iX,iY, rhoEmpty );
                    param.setForce ( iX,iY, Array<T,2> ( T(),T() ) );
                    param.setMomentum ( iX,iY, Array<T,2> ( T(),T() ) );
                    param.outsideDensity ( iX,iY ) = referenceDensity;

                    if ( model!=freeSurface )
                    {
                        param.mass2 ( iX,iY ) = density2;
                    }
                }

                if ( model!=freeSurface )
                {
                    param.attributeDynamics2 ( iX,iY, dynamicsTemplate2->clone() );
                    iniCellAtEquilibrium ( param.cell2 ( iX,iY ), density2, momentum2 );
                    param.setDensity2 ( iX,iY, density2 );
                    param.setMomentum2 ( iX,iY, momentum2 );
                    param.setForce2 ( iX,iY, Array<T,2> ( T(),T() ) );
                }
            }
        }
    }
}


/* *************** Class TwoPhasePunchRectangle2D ******************************************* */

template< typename T,template<typename U> class Descriptor>
void TwoPhasePunchRectangle2D<T,Descriptor>
::processGenericBlocks ( Box2D domain, std::vector<AtomicBlock2D*> atomicBlocks )
{
    using namespace twoPhaseFlag;
    TwoPhaseProcessorParam2D<T,Descriptor> param ( atomicBlocks );
    typedef Descriptor<T> D;

    Dot2D offset = param.absOffset();

    for ( plint iX=domain.x0; iX<=domain.x1; ++iX )
    {
        for ( plint iY=domain.y0; iY<=domain.y1; ++iY )
        {
            plint px = iX+offset.x;
            plint py = iY+offset.y;
            if ( contained ( px,py,rectangle ) )
            {
                plint numContacts = 0;
                if ( px==rectangle.x0 )
                {
                    ++numContacts;
                }
                else if ( px==rectangle.x1 )
                {
                    ++numContacts;
                }
                if ( py==rectangle.y0 )
                {
                    ++numContacts;
                }
                else if ( py==rectangle.y1 )
                {
                    ++numContacts;
                }
                bool isInterface = false;
                T volumeFraction = 0.;
                if ( numContacts>0 )
                {
                    isInterface = true;
                    switch ( numContacts )
                    {
                    case 1:
                        volumeFraction=0.5;
                        break;
                    case 2:
                        volumeFraction=3./4.;
                        break;
                    case 3:
                        volumeFraction=7./8.;
                        break;
                    default:
                        PLB_ASSERT ( false );
                    }
                }
                T density = param.getDensity ( iX,iY );
                Array<T,2> momentum ( param.getMomentum ( iX,iY ) );
                if ( isEmpty ( param.flag ( iX,iY ) ) )
                {
                    density = param.getDensity2 ( iX,iY );
                    momentum = param.getMomentum2 ( iX,iY );
                    if ( model==dynamic )
                    {
                        T deltaP = density-referenceDensity;
                        density = referenceDensity + deltaP*densityRatio;
                    }
                }
                if ( model==constRho )
                {
                    momentum *= referenceDensity/density;
                    density = referenceDensity;
                }
                T density2 ( density );
                Array<T,2> momentum2 ( momentum );
                if ( model==dynamic )
                {
                    T deltaP = density-referenceDensity;
                    density2 = referenceDensity + deltaP/densityRatio;
                    momentum2 *= density2/density;
                }
                if ( isInterface )
                {
                    param.flag ( iX,iY ) = interface;
                    param.volumeFraction ( iX,iY ) = volumeFraction;
                    param.mass ( iX,iY ) = volumeFraction*density;
                    param.setDensity ( iX,iY, density );
                    param.outsideDensity ( iX,iY ) = referenceDensity;

                    // twophase
                    param.mass2 ( iX,iY ) = ( 1.-volumeFraction ) *density2;
                }
                else
                {
                    param.flag ( iX,iY ) = empty;
                    param.attributeDynamics ( iX,iY, new NoDynamics<T,Descriptor>() );
                    param.mass ( iX,iY ) = T();
                    param.volumeFraction ( iX,iY ) = T();
                    param.setDensity ( iX,iY, rhoEmpty );
                    param.setForce ( iX,iY, Array<T,2> ( T(),T() ) );
                    param.setMomentum ( iX,iY, Array<T,2> ( T(),T() ) );
                    param.outsideDensity ( iX,iY ) = referenceDensity;

                    // twophase
                    param.mass2 ( iX,iY ) = density2;
                }
                // twophase
                param.attributeDynamics2 ( iX,iY, dynamicsTemplate2->clone() );
                iniCellAtEquilibrium ( param.cell2 ( iX,iY ), density2, momentum2 );
                param.setDensity2 ( iX,iY, density2 );
                param.setMomentum2 ( iX,iY, momentum2 );
                param.setForce2 ( iX,iY, Array<T,2> ( T(),T() ) );
            }
        }
    }
}


/* *************** Class DefaultInitializeTwoPhase2D ******************************************* */

template< typename T,template<typename U> class Descriptor>
void DefaultInitializeTwoPhase2D<T,Descriptor>
::processGenericBlocks ( Box2D domain, std::vector<AtomicBlock2D*> atomicBlocks )
{
    typedef typename TwoPhaseInterfaceLists<T,Descriptor>::Node Node;
    using namespace twoPhaseFlag;
    TwoPhaseProcessorParam2D<T,Descriptor> param ( atomicBlocks );
    typedef Descriptor<T> D;

    // In the following, spot the interface cells and tag them.
    for ( plint iX=domain.x0; iX<=domain.x1; ++iX )
    {
        for ( plint iY=domain.y0; iY<=domain.y1; ++iY )
        {
            if ( isEmpty ( param.flag ( iX,iY ) ) )
            {
                for ( plint iPop=1; iPop<D::q; iPop++ )
                {
                    plint nextX = iX + D::c[iPop][0];
                    plint nextY = iY + D::c[iPop][1];
                    // Note: there is no conflict of concurrent read/write on param.flag, because
                    // read tests is a cell is fluid, and write only converts empty cells to interface.
                    if ( isFullWet ( param.flag ( nextX,nextY ) ) )
                    {
                        param.flag ( iX,iY ) = interface;
                    }
                }
            }
        }
    }
    for ( plint iX=domain.x0; iX<=domain.x1; ++iX )
    {
        for ( plint iY=domain.y0; iY<=domain.y1; ++iY )
        {
            Array<T,2> j ( ( T ) 0., ( T ) 0. );
            iniCellAtEquilibrium ( param.cell ( iX,iY ), rhoIni, j/rhoIni );
            Array<T,2> j2 ( ( T ) 0., ( T ) 0. );
            if ( !useFreeSurfaceLimit )
            {
                iniCellAtEquilibrium ( param.cell2 ( iX,iY ), rhoIni, j2/rhoIni );
            }
            param.setDensity ( iX,iY, rhoIni );
            param.setMomentum ( iX,iY, j );
            param.setForce ( iX,iY, g );

            switch ( param.flag ( iX,iY ) )
            {
            case fluid:
            case protect:
                param.attributeDynamics ( iX,iY, dynamicsTemplate->clone() );
                param.mass ( iX,iY ) = rhoIni;
                param.volumeFraction ( iX,iY ) = ( T ) 1.;
                break;
            case interface:
                param.attributeDynamics ( iX,iY, dynamicsTemplate->clone() );
                param.mass ( iX,iY ) = 0.5 * rhoIni;
                param.volumeFraction ( iX,iY ) = ( T ) 0.5;
                break;
            case empty:
            case protectEmpty:
                param.attributeDynamics ( iX,iY, new NoDynamics<T,Descriptor>() );
                param.setForce ( iX,iY, Array<T,2> ( 0.,0. ) );
                param.mass ( iX,iY ) = ( T ) 0.;
                param.volumeFraction ( iX,iY ) = ( T ) 0.;
                break;
            case wall:
                param.attributeDynamics ( iX,iY, new BounceBack<T,Descriptor> ( rhoIni ) );
                param.mass ( iX,iY ) = ( T ) 0.;
                param.volumeFraction ( iX,iY ) = ( T ) 0.;
                break;
            default:
                // Invalid free-surface flag.
                PLB_ASSERT ( false );
            }
            if ( !useFreeSurfaceLimit )
            {
                // twophase
                param.setDensity2 ( iX,iY, rhoIni );
                param.setMomentum2 ( iX,iY, j2 );
                param.setForce2 ( iX,iY, g2 );
                switch ( param.flag ( iX,iY ) )
                {
                case fluid:
                case protect:
                    param.attributeDynamics2 ( iX,iY, new NoDynamics<T,Descriptor>() );
                    param.setForce2 ( iX,iY, Array<T,2> ( 0.,0. ) );
                    param.mass2 ( iX,iY ) = ( T ) 0.;
                    break;
                case interface:
                    param.attributeDynamics2 ( iX,iY, dynamicsTemplate2->clone() );
                    param.mass2 ( iX,iY ) = 0.5 * rhoIni;
                    break;
                case empty:
                case protectEmpty:
                    param.attributeDynamics2 ( iX,iY, dynamicsTemplate2->clone() );
                    param.mass2 ( iX,iY ) = rhoIni;
                    break;
                case wall:
                    param.attributeDynamics2 ( iX,iY, new BounceBack<T,Descriptor> ( rhoIni ) );
                    param.mass2 ( iX,iY ) = ( T ) 0.;
                    break;
                default:
                    // Invalid free-surface flag.
                    PLB_ASSERT ( false );
                }
            }
        }
    }
}


/* *************** Class PartiallyDefaultInitializeTwoPhase2D ******************************************* */

template< typename T,template<typename U> class Descriptor>
void PartiallyDefaultInitializeTwoPhase2D<T,Descriptor>
::processGenericBlocks ( Box2D domain, std::vector<AtomicBlock2D*> atomicBlocks )
{
    typedef typename TwoPhaseInterfaceLists<T,Descriptor>::Node Node;
    using namespace twoPhaseFlag;
    TwoPhaseProcessorParam2D<T,Descriptor> param ( atomicBlocks );
    typedef Descriptor<T> D;

    // In the following, spot the interface cells and tag them.
    for ( plint iX=domain.x0; iX<=domain.x1; ++iX )
    {
        for ( plint iY=domain.y0; iY<=domain.y1; ++iY )
        {
            if ( isEmpty ( param.flag ( iX,iY ) ) )
            {
                for ( plint iPop=1; iPop<D::q; iPop++ )
                {
                    plint nextX = iX + D::c[iPop][0];
                    plint nextY = iY + D::c[iPop][1];
                    // Note: there is no conflict of concurrent read/write on param.flag, because
                    // read tests is a cell is fluid, and write only converts empty cells to interface.
                    if ( isFullWet ( param.flag ( nextX,nextY ) ) )
                    {
                        param.flag ( iX,iY ) = interface;
                        param.volumeFraction ( iX,iY ) = ( T ) 0.;
                    }
                }
            }
        }
    }
    for ( plint iX=domain.x0; iX<=domain.x1; ++iX )
    {
        for ( plint iY=domain.y0; iY<=domain.y1; ++iY )
        {
            Array<T,2> j ( ( T ) 0., ( T ) 0. );
            iniCellAtEquilibrium ( param.cell ( iX,iY ), rhoIni, j/rhoIni );
            param.setDensity ( iX,iY, rhoIni );
            param.setMomentum ( iX,iY, j );
            param.setForce ( iX,iY, g );

            switch ( param.flag ( iX,iY ) )
            {
            case fluid:
            case protect:
                param.attributeDynamics ( iX,iY, dynamicsTemplate->clone() );
                param.mass ( iX,iY ) = rhoIni;
                param.volumeFraction ( iX,iY ) = ( T ) 1.;
                break;
            case interface:
                param.attributeDynamics ( iX,iY, dynamicsTemplate->clone() );
                param.mass ( iX,iY ) = rhoIni * param.volumeFraction ( iX,iY );
                break;
            case empty:
            case protectEmpty:
                param.attributeDynamics ( iX,iY, new NoDynamics<T,Descriptor>() );
                param.setForce ( iX,iY, Array<T,2> ( 0.,0. ) );
                param.mass ( iX,iY ) = ( T ) 0.;
                break;
            case wall:
                param.attributeDynamics ( iX,iY, new BounceBack<T,Descriptor> ( rhoIni ) );
                param.mass ( iX,iY ) = ( T ) 0.;
                break;
            default:
                // Invalid free-surface flag.
                PLB_ASSERT ( false );
            }
            if ( !useFreeSurfaceLimit )
            {
                Array<T,2> j2 ( ( T ) 0., ( T ) 0. );
                param.setDensity2 ( iX,iY, rhoIni );
                param.setMomentum2 ( iX,iY, j2 );
                param.setForce2 ( iX,iY, g2 );
                iniCellAtEquilibrium ( param.cell2 ( iX,iY ), rhoIni, j2/rhoIni );
                switch ( param.flag ( iX,iY ) )
                {
                case fluid:
                case protect:
                    param.attributeDynamics2 ( iX,iY, new NoDynamics<T,Descriptor>() );
                    param.setForce2 ( iX,iY, Array<T,2> ( 0.,0. ) );
                    param.mass2 ( iX,iY ) = ( T ) 0.;
                    break;
                case interface:
                    param.attributeDynamics2 ( iX,iY, dynamicsTemplate2->clone() );
                    param.mass2 ( iX,iY ) = rhoIni * ( 1.-param.volumeFraction ( iX,iY ) );
                    break;
                case empty:
                case protectEmpty:
                    param.attributeDynamics2 ( iX,iY, dynamicsTemplate2->clone() );
                    param.mass2 ( iX,iY ) = rhoIni;
                    break;
                case wall:
                    param.attributeDynamics2 ( iX,iY, new BounceBack<T,Descriptor> ( rhoIni ) );
                    param.mass2 ( iX,iY ) = ( T ) 0.;
                    break;
                default:
                    // Invalid free-surface flag.
                    PLB_ASSERT ( false );
                }
            }
        }
    }
}


/* *************** Class ConstantIniVelocityTwoPhase2D ******************************************* */

template<typename T,template<typename U> class Descriptor, class Function>
void ConstantIniVelocityTwoPhase2D<T,Descriptor,Function>::processGenericBlocks (
    Box2D domain,std::vector<AtomicBlock2D*> atomicBlocks )
{
    using namespace twoPhaseFlag;
    TwoPhaseProcessorParam2D<T,Descriptor> param ( atomicBlocks );

    Dot2D absOfs = param.absOffset();

    for ( plint iX=domain.x0; iX<=domain.x1; ++iX )
    {
        plint posX = iX + absOfs.x;
        for ( plint iY=domain.y0; iY<=domain.y1; ++iY )
        {
            plint posY = iY + absOfs.y;
            if ( f ( posX, posY ) )
            {
                if ( imposeToFluid1 && isWet ( param.flag ( iX,iY ) ) )
                {
                    T rho1 = param.getDensity ( iX,iY );
                    Array<T,2> j1 ( velocity*rho1 );
                    iniCellAtEquilibrium ( param.cell ( iX,iY ), rho1, velocity );
                    param.setMomentum ( iX,iY, j1 );
                }

                if ( imposeToFluid2 && !useFreeSurfaceLimit &&
                        ( param.flag ( iX,iY ) == empty || param.flag ( iX,iY ) == interface ) )
                {
                    T rho2 = param.getDensity2 ( iX,iY );
                    Array<T,2> j2 ( velocity*rho2 );
                    iniCellAtEquilibrium ( param.cell2 ( iX,iY ), rho2, velocity );
                    param.setMomentum ( iX,iY, j2 );
                }
            }
        }
    }
}


/* *************** Class TwoPhaseMassChange2D ******************************************* */

template< typename T,template<typename U> class Descriptor>
void TwoPhaseMassChange2D<T,Descriptor>::processGenericBlocks (
    Box2D domain,std::vector<AtomicBlock2D*> atomicBlocks )
{
    typedef Descriptor<T> D;
    using namespace twoPhaseFlag;
    TwoPhaseProcessorParam2D<T,Descriptor> param ( atomicBlocks );

    // This loop updates the mass, summarizing  Eq. 6/7, and Eq.8, in
    // the N. Thuerey e.a. technical report "Interactive Free Surface Fluids
    // with the Lattice Boltzmann Method".
    for ( plint iX=domain.x0; iX<=domain.x1; ++iX )
    {
        for ( plint iY=domain.y0; iY<=domain.y1; ++iY )
        {
            Cell<T,Descriptor>& cell = param.cell ( iX,iY );
            int flag = param.flag ( iX,iY );
            if ( isFullWet ( flag ) )
            {
                twoPhaseTemplates<T,Descriptor>::massExchangeFluidCell1 ( param, iX,iY );
            }
            else if ( flag==interface )
            {
                for ( plint iPop=0; iPop < D::q; ++iPop )
                {
                    plint nextX = iX + D::c[iPop][0];
                    plint nextY = iY + D::c[iPop][1];

                    int nextFlag = param.flag ( nextX,nextY );
                    plint opp = indexTemplates::opposite<D> ( iPop );
                    // Calculate mass at time t+1 on interface cell --> eq 7 Thurey's paper.
                    if ( isFullWet ( nextFlag ) )
                    {
                        param.mass ( iX,iY ) +=
                            ( cell[opp] - param.cell ( nextX,nextY ) [iPop] );
                    }
                    else if ( nextFlag==interface )
                    {
                        param.mass ( iX,iY ) +=
                            ( cell[opp] - param.cell ( nextX,nextY ) [iPop] ) *
                            0.5* ( param.volumeFraction ( nextX,nextY ) + param.volumeFraction ( iX,iY ) );
                    }
                }
            }
            if ( !useFreeSurfaceLimit )
            {
                if ( isEmpty ( flag ) )
                {
                    twoPhaseTemplates<T,Descriptor>::massExchangeFluidCell2 ( param, iX,iY );
                }
                else if ( flag==interface )
                {
                    for ( plint iPop=0; iPop < D::q; ++iPop )
                    {
                        plint nextX = iX + D::c[iPop][0];
                        plint nextY = iY + D::c[iPop][1];

                        int nextFlag = param.flag ( nextX,nextY );
                        plint opp = indexTemplates::opposite<D> ( iPop );
                        // Calculate mass at time t+1 on interface cell --> eq 7 Thurey's paper.
                        Cell<T,Descriptor>& cell2 = param.cell2 ( iX,iY );
                        if ( isEmpty ( nextFlag ) )
                        {
                            param.mass2 ( iX,iY ) +=
                                ( cell2[opp] - param.cell2 ( nextX,nextY ) [iPop] );
                        }
                        else if ( nextFlag==interface )
                        {
                            param.mass2 ( iX,iY ) +=
                                ( cell2[opp] - param.cell2 ( nextX,nextY ) [iPop] ) *
                                0.5* ( 1.-param.volumeFraction ( nextX,nextY ) + 1.-param.volumeFraction ( iX,iY ) );
                        }
                    }
                }
            }
        }
    }
}


/* *************** Class TwoPhaseComputeInterfaceLists2D ******************************************* */

template< typename T, template<typename> class Descriptor>
T TwoPhaseComputeInterfaceLists2D<T,Descriptor>::kappa = 1.e-3;

template< typename T,template<typename U> class Descriptor>
void TwoPhaseComputeInterfaceLists2D<T,Descriptor>::processGenericBlocks (
    Box2D domain,std::vector<AtomicBlock2D*> atomicBlocks )
{
    typedef Descriptor<T> D;
    typedef typename TwoPhaseInterfaceLists<T,Descriptor>::Node Node;
    typedef typename TwoPhaseInterfaceLists<T,Descriptor>::ExtrapolInfo ExtrapolInfo;
    TwoPhaseProcessorParam2D<T,Descriptor> param ( atomicBlocks );
    using namespace twoPhaseFlag;

    param.massExcess().clear();
    param.massExcess2().clear();
    param.interfaceToFluid().clear();
    param.interfaceToEmpty().clear();
    param.emptyToInterface().clear();
    param.fluidToInterface().clear();

    // interfaceToFluid needs to be computed in bulk+2.
    for ( plint iX=domain.x0-2; iX<=domain.x1+2; ++iX )
    {
        for ( plint iY=domain.y0-2; iY<=domain.y1+2; ++iY )
        {
            Node node ( iX,iY );
            // Eq. 11 in Thuerey's technical report.
            if ( param.flag ( iX,iY ) == interface ) // Interface cell.
            {
                if ( param.volumeFraction ( iX,iY ) > T ( 1 ) +kappa ) // Interface cell is filled.
                {
                    // Elements are added even if they belong to the envelope, because they may be
                    //   needed further down in the same data processor.
                    bool isAdjacentToProtected = false;
                    for ( plint iPop=1; iPop < D::q; ++iPop )
                    {
                        plint nextX = iX+D::c[iPop][0];
                        plint nextY = iY+D::c[iPop][1];

                        if ( param.flag ( nextX,nextY ) ==protectEmpty )
                        {
                            isAdjacentToProtected = true;
                            break;
                        }
                    }
                    if ( !isAdjacentToProtected )
                    {
                        param.interfaceToFluid().insert ( node );
                    }
                }
                else if ( param.volumeFraction ( iX,iY ) < kappa
                          && contained ( iX,iY,domain.enlarge ( 1 ) ) )
                {
                    // Avoid the case where an empty cell has a fluid neighbor without
                    // interface cell between them.
                    bool isAdjacentToProtected = false;
                    for ( plint iPop=1; iPop < D::q; ++iPop )
                    {
                        plint nextX = iX+D::c[iPop][0];
                        plint nextY = iY+D::c[iPop][1];
                        if ( param.flag ( nextX,nextY ) ==protect )
                        {
                            isAdjacentToProtected = true;
                            break;
                        }
                    }
                    if ( !isAdjacentToProtected )
                    {
                        PLB_ASSERT ( param.flag ( iX,iY ) != protectEmpty );
                        // Elements are added even if they belong to the envelope, because they may be
                        //   needed further down in the same data processor.
                        param.interfaceToEmpty().insert ( node );
                    }
                }
            }
        }
    }


    // Where interface cells have become fluid, neighboring cells must be prevented from
    //   being empty, because otherwise there's no interface cell between empty and fluid.
    typename std::set<Node>::iterator iToFl = param.interfaceToFluid().begin();
    for ( ; iToFl != param.interfaceToFluid().end(); ++iToFl )
    {
        // The node here may belong to the 1st envelope.
        Node node = *iToFl;
        plint iX=node[0];
        plint iY=node[1];


        for ( plint iPop=1; iPop < D::q; ++iPop )
        {
            plint nextX = iX+D::c[iPop][0];
            plint nextY = iY+D::c[iPop][1];
            Node nextNode ( nextX,nextY );

            // If one of my neighbors switches interface->fluid, then I shall be prevented
            //     from switching interface->empty at the same time step.
            if ( contained ( nextX,nextY,domain.enlarge ( 1 ) ) && param.flag ( nextX,nextY ) == interface )
            {
                param.interfaceToEmpty().erase ( nextNode );
            }
            // If one of my neighbors switches interface->fluid and I am empty I shall become
            //   interface.
            else if ( contained ( nextX,nextY,domain ) && isEmpty ( param.flag ( nextX,nextY ) ) )
            {
                param.emptyToInterface().insert ( std::make_pair ( nextNode, ExtrapolInfo() ) );
            }
        }
    }

    typename std::set<Node>::iterator iToE = param.interfaceToEmpty().begin();
    for ( ; iToE != param.interfaceToEmpty().end(); ++iToE )
    {
        Node node = *iToE;
        plint iX = node[0];
        plint iY = node[1];

        for ( plint iPop=1; iPop < D::q; ++iPop )
        {
            plint nextX = iX+D::c[iPop][0];
            plint nextY = iY+D::c[iPop][1];
            Node nextNode ( nextX,nextY );

            if ( param.flag ( nextX,nextY ) ==fluid )
            {
                PLB_ASSERT ( param.flag ( nextX,nextY ) != protectEmpty );
                param.fluidToInterface().insert ( std::make_pair ( nextNode, ExtrapolInfo() ) );
            }
        }
    }

    typename std::map<Node,ExtrapolInfo>::iterator flToI = param.fluidToInterface().begin();
    for ( ; flToI != param.fluidToInterface().end(); ++flToI )
    {
        Node node = flToI->first;
        plint iX = node[0];
        plint iY = node[1];


        if ( model != freeSurface )
        {
            // fluid->interface for phase 1 means, empty->interface for phase 2.
            // If non-bulk elements are left in the list, disregard to avoid accessing undefined neighbors.
            if ( contained ( iX,iY,domain ) )
            {
                // In case of dynamic/constRho model:
                // For initialization of the new cell, compute average density
                //   and momentum on neighbors.
                T averageDensity2 = T ( 0 );
                Array<T,2> averageMomentum2 ( T ( 0 ),T ( 0 ) );
                //Array<T,6> averagePiNeq2; averagePiNeq2.resetToZero();
                T sumWeights = ( T ) 0;
                for ( plint iPop=1; iPop<D::q; ++iPop )
                {
                    plint nextX = iX+D::c[iPop][0];
                    plint nextY = iY+D::c[iPop][1];

                    if ( isEmpty ( param.flag ( nextX,nextY ) ) ||
                            param.flag ( nextX,nextY ) ==interface )
                    {
                        //T weight = D::t[iPop]*param.volumeFraction(nextX,nextY,nextZ);
                        T weight = D::t[iPop];
                        sumWeights += weight;
                        averageDensity2 += weight * param.getDensity2 ( nextX,nextY );
                        averageMomentum2 += weight * param.getMomentum2 ( nextX,nextY );
                        //Array<T,6> nextPiNeq2;
                        //Cell<T,Descriptor>& cell2 = param.cell2(nextX,nextY,nextZ);
                        //cell2.getDynamics().computePiNeq(cell2, nextPiNeq2);
                        //averagePiNeq2 += weight*nextPiNeq2;
                    }
                }
                if ( sumWeights<1.e-6 )
                {
                    averageDensity2 = rhoDefault;
                }
                else
                {
                    T invSum = T ( 1 ) /sumWeights;
                    averageDensity2 *= invSum;
                    averageMomentum2 *= invSum;
                    //averagePiNeq2 *= invSum;
                }
                T newDensity = rhoDefault;
                if ( model==dynamic || model==bubblePressure || model==constRho )
                {
                    newDensity = averageDensity2;
                }
                // In case of kinetic model: take the density and momentum from the
                // partner field (but remember to subtract the surface tension term).
                else
                {
                    T rho1 = param.getDensity ( iX,iY );
                    newDensity = rho1 - surfaceTension * param.curvature ( iX,iY ) * D::invCs2;
                }
                // In both cases, the new velocity is initialized in the same way as it is computed in
                // TwoPhaseMacroscopic2D.
                Array<T,2> velocity ( param.getMomentum ( iX,iY ) /param.getDensity ( iX,iY ) );
                Array<T,2> velocity2 ( averageMomentum2/averageDensity2 );
                Array<T,2> velAverage ( ( velocity+densityRatio*velocity2 ) / ( 1.+densityRatio ) );
                Array<T,2> newMomentum ( velAverage*newDensity );

                flToI->second.density = newDensity;
                flToI->second.j = newMomentum;
                //flToI->second.PiNeq = averagePiNeq2;
            }
        }
    }


    // Compute density and momentum for cells that will switch state empty->interface.
    //   It is sufficient to do this is bulk+0.
    //   This loop performs read-only access to the lattice.
    typename std::map<Node,ExtrapolInfo>::iterator iEtoI = param.emptyToInterface().begin();
    for ( ; iEtoI != param.emptyToInterface().end(); ++iEtoI )
    {
        Node node = iEtoI->first;
        plint iX=node[0];
        plint iY=node[1];


        // If non-bulk elements are left in the list, disregard to avoid accessing undefined neighbors.
        // In case of dynamic/constRho model:
        // For initialization of the new cell, compute average density
        //   and momentum on neighbors.
        T averageDensity = T ( 0 );
        Array<T,2> averageMomentum ( T ( 0 ),T ( 0 ) );
        //Array<T,6> averagePiNeq; averagePiNeq.resetToZero();
        T sumWeights = ( T ) 0;
        for ( plint iPop=1; iPop<D::q; ++iPop )
        {
            plint nextX = iX+D::c[iPop][0];
            plint nextY = iY+D::c[iPop][1];

            // Warning: it is not accounted for the fact that neighbors can have excess mass. It
            //   might be good to account for this in the future.
            if ( isWet ( param.flag ( nextX,nextY ) ) )
            {
                //T weight = D::t[iPop]*param.volumeFraction(nextX,nextY,nextZ);
                T weight = D::t[iPop];
                sumWeights += weight;
                averageDensity += weight * param.getDensity ( nextX,nextY );
                averageMomentum += weight * param.getMomentum ( nextX,nextY );
                //Array<T,6> nextPiNeq;
                //Cell<T,Descriptor>& cell = param.cell(nextX,nextY,nextZ);
                //cell.getDynamics().computePiNeq(cell, nextPiNeq);
                //averagePiNeq += weight*nextPiNeq;
            }
        }
        if ( sumWeights<1.e-6 )
        {
            averageDensity = rhoDefault;
        }
        else
        {
            T invSum = T ( 1 ) /sumWeights;
            averageDensity *= invSum;
            averageMomentum *= invSum;
            //averagePiNeq *= invSum;
        }

        T newDensity = T();
        if ( model==dynamic || model==bubblePressure || dynamic==constRho )
        {
            newDensity = averageDensity;
        }
        // In case of kinetic model: take density and momentum from partner field.
        else if ( model==kinetic )
        {
            newDensity = param.getDensity2 ( iX,iY );
        }

        if ( model==freeSurface )
        {
            iEtoI->second.density = averageDensity;
            iEtoI->second.j = averageMomentum;
        }
        else
        {
            // The new velocity is initialized in the same way as it is computed in
            // TwoPhaseMacroscopic2D.
            Array<T,2> velocity ( averageMomentum/averageDensity );
            Array<T,2> velocity2 ( param.getMomentum2 ( iX,iY ) /param.getDensity2 ( iX,iY ) );
            Array<T,2> velAverage ( ( velocity+densityRatio*velocity2 ) / ( 1.+densityRatio ) );
            Array<T,2> newMomentum = velAverage*averageDensity;
            iEtoI->second.density = newDensity;
            iEtoI->second.j = newMomentum;
        }
        //iEtoI->second.PiNeq = averagePiNeq;
    }
}

template< typename T,template<typename U> class Descriptor>
void TwoPhaseEqualMassExcessReDistribution2D<T,Descriptor>
::processGenericBlocks ( Box2D domain,std::vector<AtomicBlock2D*> atomicBlocks )
{
    typedef Descriptor<T> D;
    typedef typename TwoPhaseInterfaceLists<T,Descriptor>::Node Node;
    using namespace twoPhaseFlag;
    TwoPhaseProcessorParam2D<T,Descriptor> param ( atomicBlocks );

    Box2D originalDomain ( domain );

    typename std::map<Node,T>::iterator iEle = param.massExcess().begin();
    for ( ; iEle != param.massExcess().end(); ++iEle )
    {
        Array<plint,2> node = iEle->first;
        plint iX = node[0];
        plint iY = node[1];

        // Check for valid interface neighbors to re-distribute mass
        if ( contained ( iX,iY,domain.enlarge ( 1 ) ) )
        {
            std::vector<int> indX, indY, indZ;
            plint numValidNeighbors = 0;

            // Check for interface neighbors in the LB directions.
            for ( plint iPop=1; iPop<D::q; iPop++ )
            {
                plint nextX = iX + D::c[iPop][0];
                plint nextY = iY + D::c[iPop][1];

                if ( param.flag ( nextX,nextY ) == interface )
                {
                    if ( contained ( nextX,nextY,domain ) )
                    {
                        indX.push_back ( nextX );
                        indY.push_back ( nextY );
                    }
                    numValidNeighbors++;
                }
            }

            // Mass re-distribution
            if ( numValidNeighbors != 0 )
            {
                int indSize = ( int ) indX.size();
                T massToRedistribute = iEle->second/ ( T ) numValidNeighbors;

                for ( int i = 0; i < indSize; i++ )
                {
                    int nextX = indX[i];
                    int nextY = indY[i];

                    param.mass ( nextX,nextY ) += massToRedistribute;
                    param.volumeFraction ( nextX,nextY ) =
                        param.mass ( nextX,nextY ) /
                        param.getDensity ( nextX,nextY );
                }
            }
            else
            {
                if ( contained ( iX,iY,originalDomain ) )
                {
                    param.addToLostMass ( iEle->second );
                }
            }
        }
    }

    // twophase
    typename std::map<Node,T>::iterator iEle2 = param.massExcess2().begin();
    for ( ; iEle2 != param.massExcess2().end(); ++iEle2 )
    {
        Array<plint,2> node = iEle2->first;
        plint iX = node[0];
        plint iY = node[1];


        // Check for valid interface neighbors to re-distribute mass
        if ( contained ( iX,iY,domain.enlarge ( 1 ) ) )
        {
            std::vector<int> indX, indY, indZ;
            plint numValidNeighbors = 0;

            // Check for interface neighbors in the LB directions.
            for ( plint iPop=1; iPop<D::q; iPop++ )
            {
                plint nextX = iX + D::c[iPop][0];
                plint nextY = iY + D::c[iPop][1];

                // for phase 2, "empty" means "fluid".
                if ( isEmpty ( param.flag ( nextX,nextY ) ) || param.flag ( nextX,nextY ) == interface )
                {
                    if ( contained ( nextX,nextY,domain ) )
                    {
                        indX.push_back ( nextX );
                        indY.push_back ( nextY );
                    }
                    numValidNeighbors++;
                }
            }

            // Mass re-distribution
            if ( numValidNeighbors != 0 )
            {
                int indSize = ( int ) indX.size();
                T massToRedistribute = iEle2->second/ ( T ) numValidNeighbors;

                for ( int i = 0; i < indSize; i++ )
                {
                    int nextX = indX[i];
                    int nextY = indY[i];

                    param.mass2 ( nextX,nextY ) += massToRedistribute;
                    if ( isEmpty ( param.flag ( nextX,nextY ) ) )
                    {
                        param.setDensity2 ( nextX,nextY,
                                            param.getDensity2 ( nextX,nextY ) + massToRedistribute );
                        for ( plint iPop=0; iPop<D::q; ++iPop )
                        {
                            param.cell2 ( nextX,nextY ) [iPop] += D::t[iPop]*massToRedistribute;
                        }
                    }
                }
            }
            else
            {
                if ( contained ( iX,iY,originalDomain ) )
                {
                    param.addToLostMass2 ( iEle2->second );
                }
            }
        }
    }
}


/* *************** Class TwoPhaseCompletion2D ******************************************* */

template< typename T,template<typename U> class Descriptor>
void TwoPhaseCompletion2D<T,Descriptor>
::processGenericBlocks ( Box2D domain,std::vector<AtomicBlock2D*> atomicBlocks )
{
    typedef Descriptor<T> D;
    using namespace twoPhaseFlag;
    typedef typename TwoPhaseInterfaceLists<T,Descriptor>::Node Node;

    TwoPhaseProcessorParam2D<T,Descriptor> param ( atomicBlocks );

    // In this data processor, populations are both written locally and read non-locally.
    // To guarantee data consistency, a first loop makes only read accesses and stores
    // the necessary information into the list neighborOppositePop. A second loop reads
    // from this list and assigns values to populations.
    std::map<Node, Array<T,D::q> > neighborOppositePop, neighborOppositePop2;
    for ( plint iX=domain.x0; iX<=domain.x1; ++iX )
    {
        for ( plint iY=domain.y0; iY<=domain.y1; ++iY )
        {
            if ( param.flag ( iX,iY ) == interface )
            {
                // Here we are on an interface node. The entire set of fi's is reconstructed.
                bool needsModification = false, needsModification2 = false;
                Array<T,D::q> savedPop, savedPop2;
                savedPop[0] = -2.;
                savedPop2[0] = -2.;
                for ( plint iPop=1; iPop < D::q; ++iPop )
                {
                    // This is one of the tricky points of the code
                    // we have to decide if the f_is from the neighborhood
                    // have to be re-update by using the Thurey's rule, which
                    // states that f_i's coming from nearest neighs. that are empty cells,
                    // have to be re-updated.
                    // I like the eq.   f^{in}_i(x,t+dt) = f^{out}_i(x-e_i,t);
                    // This eq. makes me think that the neigh. that I have to check
                    // (to control is status e.g. empty or fluid ?) has to be pos-c_i
                    plint prevX = iX-D::c[iPop][0];
                    plint prevY = iY-D::c[iPop][1];


                    plint opp = indexTemplates::opposite<D> ( iPop );

                    // Should I also change particle distribution function coming from
                    // bounceBack nodes? Well ideally no ... but there is for sure some
                    // cell configuration where these f_is are not well defined because
                    // they are probably coming from empty cells

                    // If the f_i[iPop] would be streamed from an empty cell
                    if ( isEmpty ( param.flag ( prevX,prevY ) ) ||
                            param.flag ( prevX,prevY ) == wall )
                    {
                        savedPop[iPop] = param.cell ( prevX,prevY ) [opp];
                        needsModification = true;
                    }
                    else
                    {
                        savedPop[iPop] = ( T )-2.;
                    }

                    if ( !useFreeSurfaceLimit )
                    {
                        if ( param.flag ( prevX,prevY ) == fluid ||
                                param.flag ( prevX,prevY ) == wall )
                        {
                            savedPop2[iPop] = param.cell2 ( prevX,prevY ) [opp];
                            needsModification2 = true;
                        }
                        else
                        {
                            savedPop2[iPop] = ( T )-2.;
                        }
                    }
                }
                if ( needsModification )
                {
                    neighborOppositePop.insert ( std::pair<Node,Array<T,D::q> > ( Node ( iX,iY ), savedPop ) );
                }
                if ( !useFreeSurfaceLimit )
                {
                    if ( needsModification2 )
                    {
                        neighborOppositePop2.insert ( std::pair<Node,Array<T,D::q> > ( Node ( iX,iY ), savedPop2 ) );
                    }
                }
            }
        }
    }

    typename std::map<Node, Array<T,D::q> >::const_iterator nodes = neighborOppositePop.begin();
    for ( ; nodes != neighborOppositePop.end(); ++nodes )
    {
        Node node = nodes->first;
        plint iX = node[0];
        plint iY = node[1];

        Array<T,D::q> neighborOppPop = nodes->second;
        for ( plint iPop=1; iPop < D::q; ++iPop )
        {
            if ( neighborOppPop[iPop]> ( T )-1. )
            {
                // Velocity is simply taken from the previous time step.
                Array<T,2> j = param.getMomentum ( iX,iY );
                T jSqr = VectorTemplate<T,Descriptor>::normSqr ( j );
                // Remember: the value of pressure on an interface node has been set in
                // TwoPhaseMacroscopic2D, and is equal to the ambient pressure for a
                // single free-surface fluid, or in the case of a binary pressure, an
                // averaged value.
                T rhoBar = Descriptor<T>::rhoBar ( param.getDensity ( iX,iY ) );
                T feq_i = param.cell ( iX,iY ).computeEquilibrium ( iPop, rhoBar, j, jSqr );
                plint opp = indexTemplates::opposite<D> ( iPop );
                T feq_opp_i = param.cell ( iX,iY ).computeEquilibrium ( opp, rhoBar, j, jSqr );
                param.cell ( iX,iY ) [iPop] = feq_i + feq_opp_i - neighborOppPop[iPop];
            }
        }
    }

    // twophase
    typename std::map<Node, Array<T,D::q> >::const_iterator nodes2 = neighborOppositePop2.begin();
    for ( ; nodes2 != neighborOppositePop2.end(); ++nodes2 )
    {
        Node node = nodes2->first;
        plint iX = node[0];
        plint iY = node[1];


        Array<T,D::q> neighborOppPop2 = nodes2->second;
        for ( plint iPop=1; iPop < D::q; ++iPop )
        {
            if ( neighborOppPop2[iPop]> ( T )-1. )
            {
                // Velocity is simply taken from the previous time step.
                Array<T,2> j = param.getMomentum2 ( iX,iY );
                T jSqr = VectorTemplate<T,Descriptor>::normSqr ( j );
                // Remember: the value of pressure on an interface node has been set in
                // F, and is equal to the ambient pressure for a
                // single free-surface fluid, or in the case of a binary pressure, an
                // averaged value.
                T rhoBar = Descriptor<T>::rhoBar ( param.getDensity2 ( iX,iY ) );
                T feq_i = param.cell2 ( iX,iY ).computeEquilibrium ( iPop, rhoBar, j, jSqr );
                plint opp = indexTemplates::opposite<D> ( iPop );
                T feq_opp_i = param.cell2 ( iX,iY ).computeEquilibrium ( opp, rhoBar, j, jSqr );
                param.cell2 ( iX,iY ) [iPop] = feq_i + feq_opp_i - neighborOppPop2[iPop];
            }
        }
    }
}


template< typename T,template<typename U> class Descriptor>
void VerifyTwoPhase2D<T,Descriptor>
::processGenericBlocks ( Box2D domain,std::vector<AtomicBlock2D*> atomicBlocks )
{
    typedef Descriptor<T> D;
    using namespace twoPhaseFlag;
    TwoPhaseProcessorParam2D<T,Descriptor> param ( atomicBlocks );

    // Save macroscopic fields in external scalars and update the mass-fraction.
    for ( plint iX=domain.x0-1; iX<=domain.x1+1; ++iX )
    {
        for ( plint iY=domain.y0-1; iY<=domain.y1+1; ++iY )
        {
            if ( isFullWet ( param.flag ( iX,iY ) ) )
            {
                for ( plint iPop=1; iPop<D::q; iPop++ )
                {
                    plint nextX = iX+D::c[iPop][0];
                    plint nextY = iY+D::c[iPop][1];
                    if ( isEmpty ( param.flag ( nextX,nextY ) ) )
                    {
                        PLB_ASSERT ( false );
                    }
                }
            }
        }
    }
}

template< typename T,template<typename U> class Descriptor>
void TwoPhaseMacroscopic2D<T,Descriptor>
::processGenericBlocks ( Box2D domain,std::vector<AtomicBlock2D*> atomicBlocks )
{
    typedef Descriptor<T> D;
    using namespace twoPhaseFlag;
    TwoPhaseProcessorParam2D<T,Descriptor> param ( atomicBlocks );

    T lostMass = param.getSumLostMass();
    plint numInterfaceCells = param.getNumInterfaceCells();
    T massPerCell = T();
    if ( numInterfaceCells>0 )
    {
        massPerCell = lostMass / ( T ) numInterfaceCells;
    }
    T massPerCell2 = T();
    if ( model!=freeSurface )
    {
        T lostMass2 = param.getSumLostMass2();
        if ( numInterfaceCells>0 )
        {
            massPerCell2 = lostMass2 / ( T ) numInterfaceCells;
        }
    }

    // Save macroscopic fields in external scalars and update the mass-fraction.
    for ( plint iX=domain.x0; iX<=domain.x1; ++iX )
    {
        for ( plint iY=domain.y0; iY<=domain.y1; ++iY )
        {
            T rhoBar, rhoBar2;
            Array<T,2> j ( 0.,0. ), j2 ( 0.,0. );
            T density=param.outsideDensity ( iX,iY ), density2=param.outsideDensity ( iX,iY );
            if ( isWet ( param.flag ( iX,iY ) ) )
            {
                momentTemplates<T,Descriptor>::get_rhoBar_j ( param.cell ( iX,iY ), rhoBar, j );
                density = Descriptor<T>::fullRho ( rhoBar );
            }
            if ( model!=freeSurface )
            {
                if ( isEmpty ( param.flag ( iX,iY ) ) || param.flag ( iX,iY ) ==interface )
                {
                    momentTemplates<T,Descriptor>::get_rhoBar_j ( param.cell2 ( iX,iY ), rhoBar2, j2 );
                    density2 = Descriptor<T>::fullRho ( rhoBar2 );
                }
            }
            if ( param.flag ( iX,iY ) ==interface )
            {
                param.mass ( iX,iY ) += massPerCell;

                T newDensity = T();
                if ( model==constRho || model==freeSurface )
                {
                    newDensity = param.outsideDensity ( iX,iY );
                }
                else
                {
                    param.mass2 ( iX,iY ) += massPerCell2;
                    T totalMass = param.mass ( iX,iY ) + param.mass2 ( iX,iY );
                    newDensity = totalMass; // So that totalVolume=1.
                }

                // A safety check to avoid division-by-zero situations
                // in subsequent code, in the V = m/rho formula.
                if ( newDensity<0.1*rhoDefault )
                {
                    newDensity = 0.1*rhoDefault;
                }
                // Default behavior: no rescaling.
                T newDensity1 = newDensity;
                T newDensity2 = newDensity;

                // With the two dynamic models, the pressure is computed in a more
                // sophisticated way, with proper rescaling when the density ratio
                // is other than 1.
                if ( model==dynamic || model==bubblePressure )
                {
                    PLB_ASSERT ( model!=freeSurface );
                    plint numP=0, numP2=0;
                    T deltaP = T(), deltaP2 = T();
                    T referenceDensity = param.outsideDensity ( iX,iY );
                    for ( plint dx=-1; dx<=1; ++dx )
                    {
                        for ( plint dy=-1; dy<=1; ++dy )
                        {
                            plint px = iX+dx;
                            plint py = iY+dy;
                            int fl = param.flag ( px,py );
                            if ( isEmpty ( fl ) )
                            {
                                T rhoBar_;
                                Array<T,2> j_;
                                momentTemplates<T,Descriptor>::get_rhoBar_j ( param.cell2 ( px,py ), rhoBar_, j_ );
                                deltaP2 += Descriptor<T>::fullRho ( rhoBar_ )-referenceDensity;
                                ++numP2;
                            }
                            else if ( isFullWet ( fl ) )
                            {
                                T rhoBar_;
                                Array<T,2> j_;
                                momentTemplates<T,Descriptor>::get_rhoBar_j ( param.cell ( px,py ), rhoBar_, j_ );
                                // The pressure offset for the surface tension is always added to fluid 1. Here,
                                // we subtract it from fluid 1's pressure so the pressure can be coupled to fluid 2.
                                deltaP += Descriptor<T>::fullRho ( rhoBar_ )-referenceDensity
                                          - surfaceTension * param.curvature ( px,py ) * D::invCs2;
                                ++numP;
                            }
                        }
                    }
                    // For algorithmic reasons, an interface cell has always one fluid and one empty
                    // neighbor, at least. However, it can be that the user has manually interfered
                    // with the fluid setup (for example by punching a hole into the domain), and
                    // violated this condition. We therefore address this possibility here explicitly
                    // in order to avoid divisions by zero: in such a case we simply keep the pressure
                    // from the previous time step.
                    if ( numP==0 )
                    {
                        deltaP = param.getDensity ( iX,iY )-referenceDensity
                                 - surfaceTension * param.curvature ( iX,iY ) * D::invCs2;
                    }
                    else
                    {
                        deltaP /= numP;
                    }
                    if ( numP2==0 )
                    {
                        deltaP2 = param.getDensity2 ( iX,iY )-referenceDensity;
                    }
                    else
                    {
                        deltaP2 /= numP2;
                    }
                    // In the dynamic model there's a two-way coupling. It is formulated in such a way
                    // that it is equal to the free-surface pressure-correction model when the density
                    // ratio is zero.
                    if ( model==dynamic )
                    {
                        newDensity1 = referenceDensity+0.5* ( deltaP+deltaP2*densityRatio );
                        newDensity2 = referenceDensity+0.5* ( densityRatio*deltaP+ ( 2.0-densityRatio ) *deltaP2 );
                    }
                    // In the bubble-pressure model, there's only a one-way coupling from fluid2 to fluid1.
                    else   // model=bubblePressure
                    {
                        newDensity1 = referenceDensity+deltaP2*densityRatio;
                        newDensity2 = referenceDensity+deltaP2;
                    }
                }
                // In the constRho model, the pressure is "constant", i.e. imposed from outside by the user.
                else if ( model==constRho )
                {
                    newDensity1 = param.outsideDensity ( iX,iY );
                    newDensity2 = param.outsideDensity ( iX,iY );
                }

                param.volumeFraction ( iX,iY ) = param.mass ( iX,iY ) /newDensity1;
                // On interface cells, adjust the pressure to the ambient pressure.
                j *= newDensity1/density;
                density = newDensity1;
                if ( model!=freeSurface )
                {
                    j2 *= newDensity2/density2;
                    density2 = newDensity2;
                }
            }

            if ( isFullWet ( param.flag ( iX,iY ) ) )
            {
                param.volumeFraction ( iX,iY ) = T ( 1 );
            }

            if ( model!=freeSurface )
            {
                // The following is the velocity coupling for all two-phase models. Its validity has
                // been verified in two-phase Couette flows.
                if ( param.flag ( iX,iY ) ==interface )
                {
                    T rho = param.getDensity ( iX,iY );
                    T rho2 = param.getDensity2 ( iX,iY );
                    Array<T,2> velocity ( j/rho );
                    Array<T,2> velocity2 ( j2/rho2 );
                    Array<T,2> velAverage ( ( velocity+densityRatio*velocity2 ) / ( 1.+densityRatio ) );
                    j = rho*velAverage;
                    j2 = rho2*velAverage;
                }
            }

            if ( isWet ( param.flag ( iX,iY ) ) )
            {
                Array<T,2> force = param.getForce ( iX,iY );
                T tau = T ( 1 ) /param.cell ( iX,iY ).getDynamics().getOmega();
                // Two comments:
                // - Here the force is multiplied by rho0 and not rho so that, under
                //   gravity, a linear pressure profile is obtained.
                // - The force is not multiplied by the volume fraction (some authors
                //   do multiply it by the volumeFraction), because there is a
                //   point-wise interpretation of quantities like momentum.
                j += rhoDefault*tau*force;
            }
            param.setDensity ( iX,iY, density );
            param.setMomentum ( iX,iY, j );

            if ( model!=freeSurface )
            {
                if ( isEmpty ( param.flag ( iX,iY ) ) || param.flag ( iX,iY ) ==interface )
                {
                    Array<T,2> force = param.getForce2 ( iX,iY );
                    T tau = T ( 1 ) /param.cell2 ( iX,iY ).getDynamics().getOmega();
                    j2 += rhoDefault*tau*force;
                }

                param.setDensity2 ( iX,iY, density2 );
                param.setMomentum2 ( iX,iY, j2 );
            }
        }
    }
}

template< typename T,template<typename U> class Descriptor>
void TwoPhaseInterfaceFilter2D<T,Descriptor>
::processGenericBlocks ( Box2D domain,std::vector<AtomicBlock2D*> atomicBlocks )
{
    typedef Descriptor<T> D;
    using namespace twoPhaseFlag;
    TwoPhaseProcessorParam2D<T,Descriptor> param ( atomicBlocks );

    std::vector< Array<plint,2> > interfaceNodes;

    // Save macroscopic fields in external scalars and update the mass-fraction.
    for ( plint iX=domain.x0; iX<=domain.x1; ++iX )
    {
        for ( plint iY=domain.y0; iY<=domain.y1; ++iY )
        {
            if ( param.flag ( iX,iY ) ==interface )
            {
                if ( contained ( iX,iY, domain ) )
                {
                    interfaceNodes.push_back ( Array<plint,2> ( iX,iY ) );
                }
            }
        }
    }

    std::vector<T> newRho, newRho2;
    std::vector<Array<T,2> > newJ, newJ2;
    for ( pluint i=0; i<interfaceNodes.size(); ++i )
    {
        Array<plint,2> node = interfaceNodes[i];
        plint iX=node[0];
        plint iY=node[1];

        T rho = T(), rho2 = T();
        Array<T,2> j, j2;
        j.resetToZero();
        j2.resetToZero();
        plint numNodes=0;
        T weight = 0.;
        for ( plint dx=-2; dx<=2; ++dx )
        {
            for ( plint dy=-2; dy<=2; ++dy )
            {
                plint px = iX+dx;
                plint py = iY+dy;
                if ( param.flag ( px,py ) ==interface )
                {
                    T newWeight = T();
                    if ( dx==0 && dy==0 )
                    {
                        newWeight = 2.0;
                    }
                    else
                    {
                        newWeight = 1./norm ( Array<T,2> ( dx,dy ) );
                    }
                    rho += newWeight*param.getDensity ( px,py );
                    rho2 += newWeight*param.getDensity2 ( px,py );
                    j += newWeight*param.getMomentum ( px,py );
                    j2 += newWeight*param.getMomentum2 ( px,py );
                    weight += newWeight;
                    ++numNodes;
                }
            }
        }
        PLB_ASSERT ( numNodes>0 );
        rho /= weight;
        rho2 /= weight;
        j /= weight;
        j2 /= weight;
        newRho.push_back ( rho );
        newRho2.push_back ( rho2 );
        newJ.push_back ( j );
        newJ2.push_back ( j2 );
    }
    for ( pluint i=0; i<interfaceNodes.size(); ++i )
    {
        Array<plint,2> node = interfaceNodes[i];
        plint iX=node[0];
        plint iY=node[1];


        param.setDensity ( iX,iY, newRho[i] );
        param.setDensity2 ( iX,iY, newRho2[i] );
        param.setMomentum ( iX,iY, newJ[i] );
        param.setMomentum2 ( iX,iY, newJ2[i] );
    }
}

/* *************** Class TwoPhaseIniInterfaceToAnyNodes2D ******************************************* */

template< typename T,template<typename U> class Descriptor>
void TwoPhaseIniInterfaceToAnyNodes2D<T,Descriptor>
::processGenericBlocks ( Box2D domain,std::vector<AtomicBlock2D*> atomicBlocks )
{
    typedef Descriptor<T> D;
    typedef typename TwoPhaseInterfaceLists<T,Descriptor>::Node Node;
    typedef typename TwoPhaseInterfaceLists<T,Descriptor>::ExtrapolInfo ExtrapolInfo;
    using namespace twoPhaseFlag;

    TwoPhaseProcessorParam2D<T,Descriptor> param ( atomicBlocks );

    // 1. For interface->fluid nodes, update in the flag matrix,
    //   and compute and store mass excess from these cells.
    typename std::set<Node>::iterator iEle = param.interfaceToFluid().begin();
    for ( ; iEle != param.interfaceToFluid().end(); ++iEle )
    {
        Node node = *iEle;

        plint iX = node[0];
        plint iY = node[1];

        if ( contained ( iX,iY,domain.enlarge ( 1 ) ) )
        {
            bool isAdjacentToProtected = false;
            for ( plint iPop=1; iPop < D::q; ++iPop )
            {
                plint nextX = iX+D::c[iPop][0];
                plint nextY = iY+D::c[iPop][1];
                if ( param.flag ( nextX,nextY ) ==protectEmpty )
                {
                    isAdjacentToProtected = true;
                    break;
                }
            }
            if ( !isAdjacentToProtected )
            {
                T massExcess = param.mass ( iX,iY ) - param.getDensity ( iX,iY );
                param.massExcess().insert ( std::pair<Node,T> ( node,massExcess ) );
                param.mass ( iX,iY ) = param.getDensity ( iX,iY );
                param.volumeFraction ( iX,iY ) = ( T ) 1;
                PLB_ASSERT ( param.flag ( iX,iY ) != protectEmpty );
                param.flag ( iX,iY ) = fluid;

                if ( model!=freeSurface )
                {
                    // interface->fluid for phase 1 means interface->empty for phase 2.
                    param.attributeDynamics2 ( iX,iY,new NoDynamics<T,Descriptor>() );
                    param.setForce2 ( iX,iY, Array<T,2> ( T(),T() ) );
                    param.setDensity2 ( iX,iY, rhoDefault );
                    param.setMomentum2 ( iX,iY, Array<T,2> ( T(),T() ) );
                    T massExcess2 = param.mass2 ( iX,iY );
                    param.massExcess2().insert ( std::pair<Node,T> ( node,massExcess2 ) );
                    param.mass2 ( iX,iY ) = T();
                }
            }
        }
    }

    // 2. For interface->empty nodes, update in the flag matrix,
    //   and compute and store mass excess from these cells.
    iEle = param.interfaceToEmpty().begin();
    for ( ; iEle != param.interfaceToEmpty().end(); ++iEle )
    {
        Node node = *iEle;
        plint iX = node[0];
        plint iY = node[1];


        param.flag ( iX,iY ) = empty;
        param.attributeDynamics ( iX,iY, new NoDynamics<T,Descriptor>() );

        T massExcess = param.mass ( iX,iY );
        param.massExcess().insert ( std::pair<Node,T> ( node,massExcess ) );

        param.mass ( iX,iY ) = T();
        param.volumeFraction ( iX,iY ) = T();
        param.setDensity ( iX,iY, rhoDefault );
        param.setMomentum ( iX,iY, Array<T,2> ( T(),T() ) );
        if ( model!=freeSurface )
        {
            // interface->empty for phase 1 means, interface->fluid for phase 2.
            T massExcess2 = param.mass2 ( iX,iY ) - param.getDensity2 ( iX,iY );
            param.massExcess2().insert ( std::pair<Node,T> ( node,massExcess2 ) );
            param.mass2 ( iX,iY ) = param.getDensity2 ( iX,iY );
        }
    }


    // Execute the empty->interface steps for phase 2.
    typename std::map<Node,ExtrapolInfo>::iterator flToI = param.fluidToInterface().begin();
    for ( ; flToI != param.fluidToInterface().end(); ++flToI )
    {
        Node node = flToI->first;
        plint iX = node[0];
        plint iY = node[1];


        param.flag ( iX,iY ) = interface;
        if ( model != freeSurface )
        {
            if ( contained ( iX,iY,domain ) )
            {
                param.attributeDynamics2 (
                    iX,iY, dynamicsTemplate2->clone() );

                T density = flToI->second.density;
                Array<T,2> momentum ( flToI->second.j );
                //Array<T,6> PiNeq(flToI->second.PiNeq);
                //T rhoBar = D::rhoBar(density);
                //T jSqr = normSqr(momentum);
                Cell<T,Descriptor>& cell2 = param.cell2 ( iX,iY );
                //cell2.getDynamics().regularize(cell2, rhoBar, momentum, jSqr, PiNeq);
                iniCellAtEquilibrium ( cell2, density, momentum/density );
                param.setForce2 ( iX,iY, force2 );
                // Change density, but leave mass and volumeFraction at 0, as they are later
                //   recomputed (Warning: this is probably correct, but there remains a small doubt).
                param.setMomentum2 ( iX,iY, momentum );
                param.setDensity2 ( iX,iY, density );
                param.mass2 ( iX,iY ) = T();
            }
        }
    }
}

/* *************** Class TwoPhaseIniEmptyToInterfaceNodes2D ******************************************* */

template< typename T,template<typename U> class Descriptor>
void TwoPhaseIniEmptyToInterfaceNodes2D<T,Descriptor>
::processGenericBlocks ( Box2D domain,std::vector<AtomicBlock2D*> atomicBlocks )
{
    typedef Descriptor<T> D;
    typedef typename TwoPhaseInterfaceLists<T,Descriptor>::Node Node;
    typedef typename TwoPhaseInterfaceLists<T,Descriptor>::ExtrapolInfo ExtrapolInfo;
    using namespace twoPhaseFlag;
    TwoPhaseProcessorParam2D<T,Descriptor> param ( atomicBlocks );


    // Elements that have switched state empty->interface are initialized at equilibrium.
    //   It is sufficient to initialize them in bulk+0.
    //   This loop performs write-only access on the lattice.
    typename std::map<Node,ExtrapolInfo>::iterator iEtoI = param.emptyToInterface().begin();
    for ( ; iEtoI != param.emptyToInterface().end(); ++iEtoI )
    {
        Node node = iEtoI->first;

        plint iX=node[0];
        plint iY=node[1];


        // If non-bulk elements are left in the list, disregard to avoid accessing undefined neighbors.
        T newDensity = iEtoI->second.density;
        Array<T,2> newMomentum = iEtoI->second.j;

        param.attributeDynamics (
            iX,iY, dynamicsTemplate->clone() );
        //Array<T,6> PiNeq(iEtoI->second.PiNeq);
        //T rhoBar = D::rhoBar(newDensity);
        //T jSqr = normSqr(newMomentum);
        Cell<T,Descriptor>& cell = param.cell ( iX,iY );
        //cell.getDynamics().regularize(cell, rhoBar, newMomentum, jSqr, PiNeq);
        iniCellAtEquilibrium ( cell, newDensity, newMomentum/newDensity );
        param.setForce ( iX,iY, force );
        // Change density, but leave mass and volumeFraction at 0, as they are later
        //   recomputed (Warning: this is probably correct, but there remains a small doubt).
        param.setMomentum ( iX,iY, newMomentum );
        param.setDensity ( iX,iY, newDensity );
        param.mass ( iX,iY ) = T();
        param.volumeFraction ( iX,iY ) = T();
        PLB_ASSERT ( param.flag ( iX,iY ) != protectEmpty );
        param.flag ( iX,iY ) = interface;

        // empty->interface for phase 1 means fluid->interface for phase 2.
        // Nothing needs to be done here, because the fluid node has already proper values.
    }
}

/* *************** Class TwoPhaseRemoveFalseInterfaceCells2D ******************************************* */

template< typename T,template<typename U> class Descriptor>
void TwoPhaseRemoveFalseInterfaceCells2D<T,Descriptor>
::processGenericBlocks ( Box2D domain,std::vector<AtomicBlock2D*> atomicBlocks )
{
    typedef Descriptor<T> D;
    typedef typename TwoPhaseInterfaceLists<T,Descriptor>::Node Node;
    using namespace twoPhaseFlag;

    TwoPhaseProcessorParam2D<T,Descriptor> param ( atomicBlocks );

    /// In the following, the flag status of cells is read (non-locally) and
    /// modified (locally). To avoid conflict, two loops are made, the first
    /// of which reads only, and the second writes. The vectors "interfaceToFluidNodes"
    /// and "interfaceToEmptyNodes" store coordinates of nodes that will switch
    /// status.
    std::vector<Node> interfaceToFluidNodes, interfaceToEmptyNodes;
    for ( plint iX=domain.x0-1; iX<=domain.x1+1; ++iX )
    {
        for ( plint iY=domain.y0-1; iY<=domain.y1+1; ++iY )
        {
            Node node ( iX,iY );
            if ( param.flag ( iX,iY ) == interface )
            {
                bool noFluidNeighbor = true;

                for ( plint iPop=1; iPop<D::q; iPop++ )
                {
                    plint nextX = iX+D::c[iPop][0];
                    plint nextY = iY+D::c[iPop][1];

                    if ( isFullWet ( param.flag ( nextX,nextY ) ) ) noFluidNeighbor = false;
                }
                if ( noFluidNeighbor )
                {
                    bool allInterface = true;
                    for ( plint iPop=1; iPop<D::q; iPop++ )
                    {
                        plint nextX = iX+D::c[iPop][0];
                        plint nextY = iY+D::c[iPop][1];
                        int fl = param.flag ( nextX,nextY );
                        if ( fl!=interface && fl!=wall )
                        {
                            allInterface = false;
                        }
                    }
                    // By default (if it's not the case that all
                    // neighbors are interface), the interface cell is
                    // converted to empty (because it has no fluid neighbor).
                    bool convertToFluid = false;
                    if ( allInterface )
                    {
                        convertToFluid = param.volumeFraction ( iX,iY ) >=0.5;
                    }
                    if ( convertToFluid )
                    {
                        interfaceToFluidNodes.push_back ( Node ( iX,iY ) );
                        // Store the coordinates, so flag on this node
                        // can be changed in a loop outside the current one.

                        T massExcess = param.mass ( iX,iY ) - param.getDensity ( iX,iY );
                        param.massExcess().insert ( std::pair<Node,T> ( node,massExcess ) );
                        param.mass ( iX,iY ) = param.getDensity ( iX,iY );
                        param.volumeFraction ( iX,iY ) = T ( 1 );

                        if ( model!=freeSurface )
                        {
                            // interface->fluid for phase 1 means interface->empty for phase 2.
                            T massExcess2 = param.mass2 ( iX,iY );
                            param.massExcess2().insert ( std::pair<Node,T> ( node,massExcess2 ) );
                            param.mass2 ( iX,iY ) = T();
                            param.setDensity2 ( iX,iY, rhoDefault );
                            param.setMomentum2 ( iX,iY, Array<T,2> ( T(),T() ) );
                            param.attributeDynamics2 ( iX,iY,new NoDynamics<T,Descriptor>() );
                            param.setForce2 ( iX,iY, Array<T,2> ( T(),T() ) );
                        }
                    }
                    else   // convert to empty
                    {
                        interfaceToEmptyNodes.push_back ( Node ( iX,iY ) );
                        // Store the coordinates, so flag on this node
                        // can be changed in a loop outside the current one.

                        T massExcess = param.mass ( iX,iY );
                        param.massExcess().insert ( std::pair<Node,T> ( node,massExcess ) );

                        param.attributeDynamics ( iX,iY,new NoDynamics<T,Descriptor>() );
                        param.mass ( iX,iY ) = T();
                        param.volumeFraction ( iX,iY ) = T();
                        param.setForce ( iX,iY, Array<T,2> ( T(),T() ) );
                        // Don't modify density and momentum, because they are needed by the second phase.
                        param.setDensity ( iX,iY, rhoDefault );
                        param.setMomentum ( iX,iY, Array<T,2> ( T(),T() ) );

                        if ( model != freeSurface )
                        {
                            // interface->empty for phase 1 means interface->fluid for phase 2.
                            T massExcess2 = param.mass2 ( iX,iY ) - param.getDensity2 ( iX,iY );
                            param.massExcess2().insert ( std::pair<Node,T> ( node,massExcess2 ) );
                            param.mass2 ( iX,iY ) = param.getDensity2 ( iX,iY );
                        }
                    }
                }
            }
        }
    }

    for ( pluint i=0; i<interfaceToFluidNodes.size(); ++i )
    {
        Node const& pos = interfaceToFluidNodes[i];
        PLB_ASSERT ( param.flag ( pos[0],pos[1] ) != protectEmpty );
        param.flag ( pos[0],pos[1] ) = fluid;
    }
    for ( pluint i=0; i<interfaceToEmptyNodes.size(); ++i )
    {
        Node const& pos = interfaceToEmptyNodes[i];
        PLB_ASSERT ( param.flag ( pos[0],pos[1] ) != protectEmpty );
        param.flag ( pos[0],pos[1] ) = empty;
    }
}

/* *************** Class TwoPhaseOutletMaximumVolumeFraction2D ******************************************* */

template< typename T,template<typename U> class Descriptor>
void TwoPhaseOutletMaximumVolumeFraction2D<T,Descriptor>
::processGenericBlocks ( Box2D domain, std::vector<AtomicBlock2D*> atomicBlocks )
{
    using namespace twoPhaseFlag;
    TwoPhaseProcessorParam2D<T,Descriptor> param ( atomicBlocks );
    typedef Descriptor<T> D;

    for ( plint iX=domain.x0; iX<=domain.x1; ++iX )
    {
        for ( plint iY=domain.y0; iY<=domain.y1; ++iY )
        {
            // ====================
            // First part: Fluid 1.
            // ====================
            // The pressure is imposed for fluid1. The value is taken from outside-density.
            T rho0 = param.outsideDensity ( iX,iY );
            // The philosophy is that the volume-fraction for fluid 1 takes a maximum
            // value, usually 0.5. In this way, fluid 1 is prevented from filling up the
            // outflow cells and putting the algorithm into unsolvable dilemmas.
            T maximumMass = rho0*volumeFraction;
            // The maximum volume fraction in translated into a maximum mass, given that
            // the density is known.
            if ( param.mass ( iX,iY ) > maximumMass )
            {
                param.mass ( iX,iY ) = maximumMass;
            }
            param.setDensity ( iX,iY, rho0 );
            // The idea is that the unknowns due to the boundary condition are the same
            // as the unknowns due to free surface. Therefore, the boundary condition can
            // be handled by the free-surface algorithm.

            // To end with, modify the populations to impose exactly the wished values for
            // density and velocity.
            T oldRhoBar;
            Array<T,2> oldJ;
            momentTemplates<T,Descriptor>::get_rhoBar_j ( param.cell ( iX,iY ), oldRhoBar, oldJ );
            T oldJsqr = normSqr ( oldJ );

            T newRhoBar = D::rhoBar ( rho0 );
            Array<T,2> newJ = param.getMomentum ( iX,iY );
            T newJsqr = normSqr ( newJ );

            for ( plint iPop=0; iPop < D::q; ++iPop )
            {
                T feq_old_i = param.cell ( iX,iY ).computeEquilibrium ( iPop, oldRhoBar, oldJ, oldJsqr );
                T feq_new_i = param.cell ( iX,iY ).computeEquilibrium ( iPop, newRhoBar, newJ, newJsqr );
                param.cell ( iX,iY ) [iPop] += feq_new_i - feq_old_i;
            }


            // ====================
            // Second part: Fluid 2.
            // ====================
            if ( model!=freeSurface )
            {
                // Remember that the maximum volume to be used is
                // imposed by the variable volumeFraction. Fluid 1 already occupies
                // part of the volume, and Fluid 2 occupies the remaining part, as
                // computed in the following.
                T volume2 = 1.- ( volumeFraction+param.volumeFraction ( iX,iY ) );
                T maximumMass2 = rho0*volume2;
                if ( maximumMass2<0. ) maximumMass2 = 0.;
                if ( param.mass2 ( iX,iY ) > maximumMass2 )
                {
                    param.mass2 ( iX,iY ) = maximumMass2;
                }
                param.setDensity2 ( iX,iY, rho0 );
                // For fluid 2, the free surface is not equal to the interface of the
                // boundary condition. As a matter of fact, we need a boundary condition
                // for fluid 2 even when the  free surface is very far away. Therefore, in
                // the following we manually implement a boundary condition equivalent
                // to the closure scheme of the free-surface algorithm.
                Cell<T,Descriptor>& cell2 = param.cell2 ( iX,iY );
                for ( plint iPop=1; iPop < D::q; ++iPop )
                {
                    plint nextX = iX + D::c[iPop][0];
                    plint nextY = iY + D::c[iPop][1];

                    int nextFlag = param.flag ( nextX,nextY );
                    if ( nextFlag==wall )
                    {
                        plint opp = indexTemplates::opposite<D> ( iPop );
                        Array<T,2> j2 = param.getMomentum2 ( iX,iY );
                        T j2Sqr = normSqr ( j2 );
                        T rhoBar2 = D::rhoBar ( rho0 );
                        T feq_i = cell2.computeEquilibrium ( iPop, rhoBar2, j2, j2Sqr );
                        T feq_opp_i = cell2.computeEquilibrium ( opp, rhoBar2, j2, j2Sqr );
                        cell2[opp] = -param.cell2 ( nextX,nextY ) [iPop]+feq_i+feq_opp_i;
                    }
                }
                Array<T,2> j2;
                T rhoBar2;
                momentTemplates<T,Descriptor>::get_rhoBar_j ( cell2, rhoBar2, j2 );
                param.setMomentum2 ( iX,iY, j2 );
            }
        }
    }
}

/* *************** Class TwoPhaseComputePressure2D ******************************************* */

template< typename T,template<typename U> class Descriptor>
void TwoPhaseComputePressure2D<T,Descriptor>
::processGenericBlocks ( Box2D domain, std::vector<AtomicBlock2D*> atomicBlocks )
{
    using namespace twoPhaseFlag;
    ScalarField2D<T>* pressure =
        ( model==freeSurface ) ?
        dynamic_cast<ScalarField2D<T>*> ( atomicBlocks[10] ) :
        dynamic_cast<ScalarField2D<T>*> ( atomicBlocks[14] );
    PLB_ASSERT ( pressure );
    TwoPhaseProcessorParam2D<T,Descriptor> param ( atomicBlocks );
    typedef Descriptor<T> D;
    Dot2D offset = computeRelativeDisplacement ( *atomicBlocks[0], *pressure );

    for ( plint iX=domain.x0; iX<=domain.x1; ++iX )
    {
        for ( plint iY=domain.y0; iY<=domain.y1; ++iY )
        {
            if ( param.flag ( iX,iY ) ==empty )
            {
                if ( computeFluid2 && model!=freeSurface )
                {
                    pressure->get ( iX+offset.x,iY+offset.y ) =
                        ( param.getDensity2 ( iX,iY )-param.outsideDensity ( iX,iY ) ) *densityRatio*D::cs2;
                }
                else
                {
                    pressure->get ( iX+offset.x,iY+offset.y ) = T();
                }
            }
            else
            {
                if ( computeFluid1 )
                {
                    pressure->get ( iX+offset.x,iY+offset.y ) =
                        ( param.getDensity ( iX,iY )-param.outsideDensity ( iX,iY ) ) *D::cs2;
                }
                else
                {
                    pressure->get ( iX+offset.x,iY+offset.y ) = T();
                }
            }
        }
    }
}

/* *************** Class TwoPhaseComputeVelocity2D ******************************************* */

template< typename T,template<typename U> class Descriptor>
void TwoPhaseComputeVelocity2D<T,Descriptor>
::processGenericBlocks ( Box2D domain, std::vector<AtomicBlock2D*> atomicBlocks )
{
    using namespace twoPhaseFlag;
    TensorField2D<T,2>* velocity =
        useFreeSurfaceLimit ?
        dynamic_cast<TensorField2D<T,2>*> ( atomicBlocks[10] ) :
        dynamic_cast<TensorField2D<T,2>*> ( atomicBlocks[14] );
    PLB_ASSERT ( velocity );
    TwoPhaseProcessorParam2D<T,Descriptor> param ( atomicBlocks );
    typedef Descriptor<T> D;
    Dot2D offset = computeRelativeDisplacement ( *atomicBlocks[0], *velocity );

    for ( plint iX=domain.x0; iX<=domain.x1; ++iX )
    {
        for ( plint iY=domain.y0; iY<=domain.y1; ++iY )
        {
            if ( param.flag ( iX,iY ) ==empty )
            {
                if ( computeFluid2 && !useFreeSurfaceLimit )
                {
                    velocity->get ( iX+offset.x,iY+offset.y ) =
                        param.getMomentum2 ( iX,iY ) /param.getDensity2 ( iX,iY );
                }
                else
                {
                    velocity->get ( iX+offset.x,iY+offset.y ).resetToZero();
                }
            }
            else
            {
                if ( computeFluid1 )
                {
                    velocity->get ( iX+offset.x,iY+offset.y ) =
                        param.getMomentum ( iX,iY ) /param.getDensity ( iX,iY );
                }
                else
                {
                    velocity->get ( iX+offset.x,iY+offset.y ).resetToZero();
                }
            }
        }
    }
}

/* *************** Class TwoPhaseAverageVelocity2D ******************************************* */

template< typename T,template<typename U> class Descriptor>
TwoPhaseAverageVelocity2D<T,Descriptor>::TwoPhaseAverageVelocity2D ( TwoPhaseModel model_ )
    : sumVel1_ID (
        this->getStatistics().subscribeSum(),
        this->getStatistics().subscribeSum() ),
    sumVel2_ID (
        this->getStatistics().subscribeSum(),
        this->getStatistics().subscribeSum() ),
    weight1_ID ( this->getStatistics().subscribeSum() ),
    weight2_ID ( this->getStatistics().subscribeSum() ),
    model ( model_ )
{ }

template< typename T,template<typename U> class Descriptor>
void TwoPhaseAverageVelocity2D<T,Descriptor>
::processGenericBlocks ( Box2D domain, std::vector<AtomicBlock2D*> atomicBlocks )
{
    using namespace twoPhaseFlag;
    TwoPhaseProcessorParam2D<T,Descriptor> param ( atomicBlocks );
    typedef Descriptor<T> D;
    BlockStatistics& statistics = this->getStatistics();

    for ( plint iX=domain.x0; iX<=domain.x1; ++iX )
    {
        for ( plint iY=domain.y0; iY<=domain.y1; ++iY )
        {
            int fl = param.flag ( iX,iY );
            T vf = param.volumeFraction ( iX,iY );
            if ( isWet ( fl ) )
            {
                Array<T,2> vel1 =
                    param.getMomentum ( iX,iY ) /param.getDensity ( iX,iY );
                statistics.gatherSum ( sumVel1_ID[0], vel1[0]*vf );
                statistics.gatherSum ( sumVel1_ID[1], vel1[1]*vf );
                statistics.gatherSum ( weight1_ID, vf );
            }
            if ( model!=freeSurface )
            {
                if ( isEmpty ( fl ) || fl==interface )
                {
                    Array<T,2> vel2 =
                        param.getMomentum2 ( iX,iY ) /param.getDensity2 ( iX,iY );
                    statistics.gatherSum ( sumVel2_ID[0], vel2[0]*vf );
                    statistics.gatherSum ( sumVel2_ID[1], vel2[1]*vf );
                    statistics.gatherSum ( weight2_ID, vf );
                }
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
Array<T,2> TwoPhaseAverageVelocity2D<T,Descriptor>::getAverageVelocity() const
{
    T weight = this->getStatistics().getSum ( weight1_ID );
    Array<T,2> vel (
        this->getStatistics().getSum ( sumVel1_ID[0] ) /weight,
        this->getStatistics().getSum ( sumVel1_ID[1] ) /weight );
    return vel;
}

template<typename T, template<typename U> class Descriptor>
Array<T,2> TwoPhaseAverageVelocity2D<T,Descriptor>::getAverageVelocity2() const
{
    T weight = this->getStatistics().getSum ( weight2_ID );
    Array<T,2> vel (
        this->getStatistics().getSum ( sumVel2_ID[0] ) /weight,
        this->getStatistics().getSum ( sumVel2_ID[1] ) /weight );
    return vel;
}


/* *************** Class TwoPhaseAveragePressure2D ******************************************* */

template< typename T,template<typename U> class Descriptor>
TwoPhaseAveragePressure2D<T,Descriptor>::TwoPhaseAveragePressure2D ( T densityRatio_, T rhoDefault_, TwoPhaseModel model_ )
    : sumRho1_ID ( this->getStatistics().subscribeSum() ),
      sumRho2_ID ( this->getStatistics().subscribeSum() ),
      weight1_ID ( this->getStatistics().subscribeSum() ),
      weight2_ID ( this->getStatistics().subscribeSum() ),
      densityRatio ( densityRatio_ ),
      rhoDefault ( rhoDefault_ ),
      model ( model_ )
{ }

template< typename T,template<typename U> class Descriptor>
void TwoPhaseAveragePressure2D<T,Descriptor>
::processGenericBlocks ( Box2D domain, std::vector<AtomicBlock2D*> atomicBlocks )
{
    using namespace twoPhaseFlag;
    TwoPhaseProcessorParam2D<T,Descriptor> param ( atomicBlocks );
    typedef Descriptor<T> D;
    BlockStatistics& statistics = this->getStatistics();

    for ( plint iX=domain.x0; iX<=domain.x1; ++iX )
    {
        for ( plint iY=domain.y0; iY<=domain.y1; ++iY )
        {
            int fl = param.flag ( iX,iY );
            T vf = param.volumeFraction ( iX,iY );
            if ( isWet ( fl ) )
            {
                T rho1 = param.getDensity ( iX,iY );
                statistics.gatherSum ( sumRho1_ID, rho1*vf );
                statistics.gatherSum ( weight1_ID, vf );
            }
            if ( model!=freeSurface )
            {
                if ( isEmpty ( fl ) || fl==interface )
                {
                    T rho2 = param.getDensity2 ( iX,iY );
                    if ( model==dynamic || model==bubblePressure || model==constRho )
                        // In the dynamic models, only the pressure fluctuation (rho-referenceDensity)
                        // is rescaled to the pressure units of fluid 1. The static component
                        // ("referenceDensity") is already expressed in units of fluid 1. It is for
                        // example equal to the pressure correction computed by pattern matching to
                        // cope with varying bubble volumes.
                    {
                        T referenceDensity = param.outsideDensity ( iX,iY );
                        rho2 = referenceDensity + densityRatio* ( rho2-referenceDensity );
                    }
                    statistics.gatherSum ( sumRho2_ID, rho2* ( 1.-vf ) );
                    statistics.gatherSum ( weight2_ID, 1.-vf );
                }
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
T TwoPhaseAveragePressure2D<T,Descriptor>::getAveragePressure() const
{
    typedef Descriptor<T> D;
    T rho = this->getStatistics().getSum ( sumRho1_ID ) /this->getStatistics().getSum ( weight1_ID );
    return D::cs2* ( rho-rhoDefault );
}

template<typename T, template<typename U> class Descriptor>
T TwoPhaseAveragePressure2D<T,Descriptor>::getAveragePressure2() const
{
    typedef Descriptor<T> D;
    T rho = this->getStatistics().getSum ( sumRho2_ID ) /this->getStatistics().getSum ( weight2_ID );
    bool dynamicModel = ( model==dynamic || model==bubblePressure || model==constRho );
    if ( dynamicModel )
    {
        return D::cs2* ( rho-rhoDefault );
    }
    else
    {
        // In the quasi-static model, the full pressure term is rescaled to the pressure
        // units of fluid 1. In case of a non-unity density ratio, there will be a
        // pressure jump through the interface: this is the major deficiency of the
        // quasi-static model.
        return D::cs2*densityRatio* ( rho-rhoDefault );
    }
}

}  // namespace plb

#endif  // TWO_PHASE_MODEL_2D_HH

