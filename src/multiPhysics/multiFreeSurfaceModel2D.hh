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

#ifndef MULTI_FREE_SURFACE_MODEL_2D_HH
#define MULTI_FREE_SURFACE_MODEL_2D_HH

#include "core/globalDefs.h"
#include "core/block2D.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include "atomicBlock/dataProcessor2D.h"
#include "atomicBlock/blockLattice2D.h"
#include "atomicBlock/atomicContainerBlock2D.h"
#include "multiPhysics/freeSurfaceModel2D.h"
#include "multiPhysics/freeSurfaceTemplates.h"
#include "multiPhysics/multiFreeSurfaceModel2D.h"

#include <cmath>

namespace plb
{

/* *************** Class MultiFreeSurfaceOneWayCoupling2D ******************************************* */

template<typename T,template<typename U> class Descriptor>
void MultiFreeSurfaceOneWayCoupling2D<T,Descriptor>
::processGenericBlocks ( Box2D domain,std::vector<AtomicBlock2D*> atomicBlocks )
{
    PLB_ASSERT ( atomicBlocks.size() == 20 );

    using namespace twoPhaseFlag;

    std::vector<AtomicBlock2D*> atomicBlocks1, atomicBlocks2;

    for ( int i = 0; i < 10; i++ )
        atomicBlocks1.push_back ( atomicBlocks[i] );

    for ( int i = 0; i < 10; i++ )
        atomicBlocks2.push_back ( atomicBlocks[10 + i] );

    FreeSurfaceProcessorParam2D<T,Descriptor> param1 ( atomicBlocks1 );
    FreeSurfaceProcessorParam2D<T,Descriptor> param2 ( atomicBlocks2 );

    for ( plint iX=domain.x0; iX<=domain.x1; ++iX )
    {
        for ( plint iY=domain.y0; iY<=domain.y1; ++iY )
        {
            if ( param1.flag ( iX,iY ) ==interface && isWet ( param2.flag ( iX,iY ) ) )
            {
                // Velocity coupling
                T rhoBar1;
                Array<T,2> j1;
                momentTemplates<T,Descriptor>::get_rhoBar_j ( param1.cell ( iX,iY ), rhoBar1, j1 );
                T density1 = Descriptor<T>::fullRho ( rhoBar1 );

                T rhoBar2;
                Array<T,2> j2;
                momentTemplates<T,Descriptor>::get_rhoBar_j ( param2.cell ( iX,iY ), rhoBar2, j2 );
                T density2 = Descriptor<T>::fullRho ( rhoBar2 );

                T tau = T ( 1 ) /param1.cell ( iX,iY ).getDynamics().getOmega();
                for ( plint iD=0; iD<3; ++iD )
                {
                    //j1[iD] += interactionStrength*density1*(j2[iD]/density2 - j1[iD]/density1)*tau;
                    j1[iD] += interactionStrength*rhoDefault1* ( j2[iD]/density2 - j1[iD]/density1 ) *tau;
                }

                param1.setMomentum ( iX,iY, j1 );
            }
        }
    }
}

/* *************** Class MultiFreeSurfaceVelocityContinuityCoupling2D ******************************************* */

template<typename T,template<typename U> class Descriptor>
void MultiFreeSurfaceVelocityContinuityCoupling2D<T,Descriptor>
::processGenericBlocks ( Box2D domain,std::vector<AtomicBlock2D*> atomicBlocks )
{
    PLB_ASSERT ( atomicBlocks.size() == 20 );

    using namespace twoPhaseFlag;

    std::vector<AtomicBlock2D*> atomicBlocks1, atomicBlocks2;

    for ( int i = 0; i < 10; i++ )
        atomicBlocks1.push_back ( atomicBlocks[i] );

    for ( int i = 0; i < 10; i++ )
        atomicBlocks2.push_back ( atomicBlocks[10 + i] );

    FreeSurfaceProcessorParam2D<T,Descriptor> param1 ( atomicBlocks1 );
    FreeSurfaceProcessorParam2D<T,Descriptor> param2 ( atomicBlocks2 );

    T convertRho1to2 = rhoDefault2/rhoDefault1;
    T convertRho2to1 = rhoDefault1/rhoDefault2;

    for ( plint iX=domain.x0; iX<=domain.x1; ++iX )
    {
        for ( plint iY=domain.y0; iY<=domain.y1; ++iY )
        {
            Array<T,2> velocity ( ( T ) 0.0, ( T ) 0.0 );

            // Enforce continuity of velocity only when both fluid cells belong to the interface.
            if ( param1.flag ( iX,iY ) ==interface && param2.flag ( iX,iY ) ==interface )
            {
                //if ((param2.flag(iX,iY)==interface && isWet(param1.flag(iX,iY))) ||
                //    (param1.flag(iX,iY)==interface && isWet(param2.flag(iX,iY)))) {
                T rhoBar1;
                Array<T,2> j1;
                momentTemplates<T,Descriptor>::get_rhoBar_j ( param1.cell ( iX,iY ), rhoBar1, j1 );
                T density1 = Descriptor<T>::fullRho ( rhoBar1 );

                T rhoBar2;
                Array<T,2> j2;
                momentTemplates<T,Descriptor>::get_rhoBar_j ( param2.cell ( iX,iY ), rhoBar2, j2 );
                T density2 = Descriptor<T>::fullRho ( rhoBar2 );

                velocity = 0.5 * ( j1/density1 + j2/density2 );
            }

            if ( param1.flag ( iX,iY ) ==interface && isWet ( param2.flag ( iX,iY ) ) )
            {
                // Compute the average pressure from fluid cells of both fluids.
                T averageDensity = T();
                plint numNeighbors = 0;
                for ( plint i=-1; i<2; i++ )
                {
                    for ( plint j=-1; j<2; j++ )
                    {
                        if ( isFullWet ( param1.flag ( iX+i,iY+j ) ) )
                        {
                            ++numNeighbors;
                            averageDensity += param1.getDensity ( iX+i,iY+j );
                        }
                        if ( isFullWet ( param2.flag ( iX+i,iY+j ) ) )
                        {
                            ++numNeighbors;
                            averageDensity += convertRho2to1 * param2.getDensity ( iX+i,iY+j );
                        }
                    }
                }
                if ( numNeighbors==0 )
                {
                    averageDensity = rhoDefault1;
                }
                else
                {
                    averageDensity /= ( T ) numNeighbors;
                }
                param1.setDensity ( iX,iY, averageDensity );

                param1.volumeFraction ( iX,iY ) = param1.mass ( iX,iY ) / averageDensity;
                if ( param2.flag ( iX,iY ) == interface )
                {
                    param1.setMomentum ( iX,iY, param1.getDensity ( iX,iY ) * velocity );
                }
            }

            if ( param2.flag ( iX,iY ) ==interface && isWet ( param1.flag ( iX,iY ) ) )
            {
                // Compute the average pressure from fluid cells of both fluids.
                T averageDensity = T();
                plint numNeighbors = 0;
                for ( plint i=-1; i<2; i++ )
                {
                    for ( plint j=-1; j<2; j++ )
                    {
                        if ( isFullWet ( param1.flag ( iX+i,iY+j ) ) )
                        {
                            ++numNeighbors;
                            averageDensity += convertRho1to2 * param1.getDensity ( iX+i,iY+j );
                        }
                        if ( isFullWet ( param2.flag ( iX+i,iY+j ) ) )
                        {
                            ++numNeighbors;
                            averageDensity += param2.getDensity ( iX+i,iY+j );
                        }
                    }
                }
                if ( numNeighbors==0 )
                {
                    averageDensity = rhoDefault2;
                }
                else
                {
                    averageDensity /= ( T ) numNeighbors;
                }
                param2.setDensity ( iX,iY, averageDensity );

                param2.volumeFraction ( iX,iY ) = param2.mass ( iX,iY ) / averageDensity;
                if ( param1.flag ( iX,iY ) == interface )
                {
                    param2.setMomentum ( iX,iY, param2.getDensity ( iX,iY ) * velocity );
                }
            }
        }
    }
}

/* *************** Class MultiFreeSurfaceRepellingForceCoupling2D ******************************************* */

template<typename T,template<typename U> class Descriptor>
T MultiFreeSurfaceRepellingForceCoupling2D<T,Descriptor>::deltaFunction ( T r, T h )
{
    Precision precision;
    if ( sizeof ( T ) == sizeof ( float ) )
        precision = FLT;
    else if ( sizeof ( T ) == sizeof ( double ) )
        precision = DBL;
    else if ( sizeof ( T ) == sizeof ( long double ) )
        precision = LDBL;
    else
        PLB_ASSERT ( false );

    T eps = getEpsilon<T> ( precision );

    PLB_ASSERT ( r > ( T ) 0 || fabs ( r ) <= eps );
    PLB_ASSERT ( h > ( T ) 0 && fabs ( h ) >  eps );

    static T pi = acos ( ( T ) -1.0 );

    T delta = 0.0;
    if ( r < h || fabs ( r - h ) <= eps )
    {
        delta = 0.5 * ( 1.0 + cos ( pi * r / h ) );
    }

    return delta;
}

template<typename T,template<typename U> class Descriptor>
void MultiFreeSurfaceRepellingForceCoupling2D<T,Descriptor>
::processGenericBlocks ( Box2D domain,std::vector<AtomicBlock2D*> atomicBlocks )
{
    PLB_ASSERT ( atomicBlocks.size() == 20 );

    using namespace twoPhaseFlag;

    std::vector<AtomicBlock2D*> atomicBlocks1, atomicBlocks2;

    for ( int i = 0; i < 10; i++ )
        atomicBlocks1.push_back ( atomicBlocks[i] );

    for ( int i = 0; i < 10; i++ )
        atomicBlocks2.push_back ( atomicBlocks[10 + i] );

    FreeSurfaceProcessorParam2D<T,Descriptor> param1 ( atomicBlocks1 );
    FreeSurfaceProcessorParam2D<T,Descriptor> param2 ( atomicBlocks2 );

    T convertRho1to2 = rhoDefault2/rhoDefault1;
    T convertRho2to1 = rhoDefault1/rhoDefault2;

    for ( plint iX=domain.x0; iX<=domain.x1; ++iX )
    {
        for ( plint iY=domain.y0; iY<=domain.y1; ++iY )
        {
            // Implement the repelling force in terms of a momentum correction. Remember that fluids can have strong
            // density differences. In order to have the units of a momentum, this correction term must be
            // multiplied by the density.
            if ( isWet ( param1.flag ( iX,iY ) ) && isWet ( param2.flag ( iX,iY ) ) )
            {
                // Add a force to both fluids. Remember the envelope is 3.
                Array<T,2> force1 ( ( T ) 0, ( T ) 0, ( T ) 0 );
                for ( plint i=-2; i<3; i++ )
                {
                    plint I = iX + i;
                    for ( plint j=-2; j<3; j++ )
                    {
                        plint J = iY + j;
                        if ( isWet ( param1.flag ( I,J ) ) && isWet ( param2.flag ( I,J ) ) )
                        {
                            T VFdiff = param1.volumeFraction ( I,J ) + param2.volumeFraction ( I,J ) - 1.0;
                            if ( VFdiff> ( T ) 0.0 )
                            {
                                Array<T,2> outwardNormalToFluid1 = param1.getNormal ( I,J ) - param2.getNormal ( I,J );
                                outwardNormalToFluid1 /= norm ( outwardNormalToFluid1 );
                                T r = sqrt ( i*i + j*j );
                                T delta = deltaFunction ( r, ( T ) 3 ); // h is 3 because envelope is 3.
                                force1 += -VFdiff*delta*outwardNormalToFluid1;
                            }
                        }
                    }
                }
                force1 *= interactionStrength;

                // force1 is an acceleration.
                T tau1 = T ( 1 ) /param1.cell ( iX,iY ).getDynamics().getOmega();
                param1.setMomentum ( iX,iY, param1.getMomentum ( iX,iY ) + rhoDefault1*tau1*force1 );
                T tau2 = T ( 1 ) /param2.cell ( iX,iY ).getDynamics().getOmega();
                param2.setMomentum ( iX,iY, param2.getMomentum ( iX,iY ) + rhoDefault2*tau2* ( -force1 ) );

                /*
                // force1 is a force.
                T tau1 = T(1)/param1.cell(iX,iY).getDynamics().getOmega();
                param1.setMomentum(iX,iY, param1.getMomentum(iX,iY) + tau1*force1);
                T tau2 = T(1)/param2.cell(iX,iY).getDynamics().getOmega();
                param2.setMomentum(iX,iY, param2.getMomentum(iX,iY) + tau2*(-force1));
                */
            }

            Array<T,2> j;

            if ( param1.flag ( iX,iY ) ==interface && isWet ( param2.flag ( iX,iY ) ) )
            {
                // Compute the average pressure from fluid cells of both fluids.
                T averageDensity = T();
                plint numNeighbors = 0;
                for ( plint i=-1; i<2; i++ )
                {
                    for ( plint j=-1; j<2; j++ )
                    {
                        if ( isFullWet ( param1.flag ( iX+i,iY+j ) ) )
                        {
                            ++numNeighbors;
                            averageDensity += param1.getDensity ( iX+i,iY+j );
                        }
                        if ( isFullWet ( param2.flag ( iX+i,iY+j ) ) )
                        {
                            ++numNeighbors;
                            averageDensity += convertRho2to1 * param2.getDensity ( iX+i,iY+j );
                        }
                    }
                }
                if ( numNeighbors==0 )
                {
                    averageDensity = rhoDefault1;
                }
                else
                {
                    averageDensity /= ( T ) numNeighbors;
                }
                T oldDensity = param1.getDensity ( iX,iY );
                param1.setDensity ( iX,iY, averageDensity );
                param1.volumeFraction ( iX,iY ) = param1.mass ( iX,iY ) / averageDensity;
                j = param1.getMomentum ( iX,iY );
                j *= averageDensity / oldDensity;
                param1.setMomentum ( iX,iY, j );


                /*
                // Implement the repelling force in terms of a momentum correction. Remember that fluids can have strong
                // density differences. In order to have the units of a momentum, this correction term must be
                // multiplied by the density.
                T VFdiff = param1.volumeFraction(iX,iY) + param2.volumeFraction(iX,iY) - 1.0;
                Array<T,2> normalToInterface = param1.getNormal(iX, iY, iZ);
                //Array<T,2> normalToInterface = -param2.getNormal(iX, iY, iZ);
                if (VFdiff>(T)0.0) {
                    T tau = T(1)/param2.cell(iX,iY).getDynamics().getOmega();
                    param2.setMomentum(iX,iY, param2.getMomentum(iX,iY) +
                            interactionStrength*rhoDefault2*VFdiff*tau*normalToInterface);
                }
                */
            }

            if ( param2.flag ( iX,iY ) ==interface && isWet ( param1.flag ( iX,iY ) ) )
            {
                // Compute the average pressure from fluid cells of both fluids.
                T averageDensity = T();
                plint numNeighbors = 0;
                for ( plint i=-1; i<2; i++ )
                {
                    for ( plint j=-1; j<2; j++ )
                    {
                        for ( plint k=-1; k<2; k++ )
                        {
                            if ( isFullWet ( param1.flag ( iX+i,iY+j ) ) )
                            {
                                ++numNeighbors;
                                averageDensity += convertRho1to2 * param1.getDensity ( iX+i,iY+j );
                            }
                            if ( isFullWet ( param2.flag ( iX+i,iY+j ) ) )
                            {
                                ++numNeighbors;
                                averageDensity += param2.getDensity ( iX+i,iY+j );
                            }
                        }
                    }
                }
                if ( numNeighbors==0 )
                {
                    averageDensity = rhoDefault2;
                }
                else
                {
                    averageDensity /= ( T ) numNeighbors;
                }

                T oldDensity = param2.getDensity ( iX,iY );
                param2.setDensity ( iX,iY, averageDensity );
                param2.volumeFraction ( iX,iY ) = param2.mass ( iX,iY ) / averageDensity;
                j = param2.getMomentum ( iX,iY );
                j *= averageDensity / oldDensity;
                param2.setMomentum ( iX,iY, j );

                /*
                // Implement the repelling force in terms of a momentum correction. Remember that fluids can have strong
                // density differences. In order to have the units of a momentum, this correction term must be
                // multiplied by the density.
                T VFdiff = param1.volumeFraction(iX,iY) + param2.volumeFraction(iX,iY) - 1.0;
                Array<T,2> normalToInterface = param2.getNormal(iX,iY);
                //Array<T,2> normalToInterface = -param1.getNormal(iX,iY);
                if (VFdiff>(T)0.0) {
                    T tau = T(1)/param1.cell(iX,iY).getDynamics().getOmega();
                    param1.setMomentum(iX,iY, param1.getMomentum(iX,iY) +
                            interactionStrength*rhoDefault1*VFdiff*tau*normalToInterface);
                }
                */
            }
        }
    }
}

/* *************** Class MultiFreeSurfaceComplexCoupling2D ******************************************* */

template<typename T,template<typename U> class Descriptor>
T MultiFreeSurfaceComplexCoupling2D<T,Descriptor>::deltaFunction ( T r, T h )
{
    Precision precision;
    if ( sizeof ( T ) == sizeof ( float ) )
        precision = FLT;
    else if ( sizeof ( T ) == sizeof ( double ) )
        precision = DBL;
    else if ( sizeof ( T ) == sizeof ( long double ) )
        precision = LDBL;
    else
        PLB_ASSERT ( false );

    T eps = getEpsilon<T> ( precision );

    PLB_ASSERT ( r > ( T ) 0 || fabs ( r ) <= eps );
    PLB_ASSERT ( h > ( T ) 0 && fabs ( h ) >  eps );

    static T pi = acos ( ( T ) -1.0 );

    T delta = 0.0;
    if ( r < h || fabs ( r - h ) <= eps )
    {
        delta = 0.5 * ( 1.0 + cos ( pi * r / h ) );
    }

    return delta;
}

template<typename T,template<typename U> class Descriptor>
void MultiFreeSurfaceComplexCoupling2D<T,Descriptor>
::processGenericBlocks ( Box2D domain,std::vector<AtomicBlock2D*> atomicBlocks )
{
    PLB_ASSERT ( atomicBlocks.size() == 20 );

    using namespace twoPhaseFlag;

    std::vector<AtomicBlock2D*> atomicBlocks1, atomicBlocks2;

    for ( int i = 0; i < 10; i++ )
        atomicBlocks1.push_back ( atomicBlocks[i] );

    for ( int i = 0; i < 10; i++ )
        atomicBlocks2.push_back ( atomicBlocks[10 + i] );

    FreeSurfaceProcessorParam2D<T,Descriptor> param1 ( atomicBlocks1 );
    FreeSurfaceProcessorParam2D<T,Descriptor> param2 ( atomicBlocks2 );

    T convertRho1to2 = rhoDefault2/rhoDefault1;
    T convertRho2to1 = rhoDefault1/rhoDefault2;

    for ( plint iX=domain.x0; iX<=domain.x1; ++iX )
    {
        for ( plint iY=domain.y0; iY<=domain.y1; ++iY )
        {
            Array<T,2> velocity ( ( T ) 0.0, ( T ) 0.0 );

            if ( param1.flag ( iX,iY ) ==interface && param2.flag ( iX,iY ) ==interface )
            {
                T rhoBar1;
                Array<T,2> j1;
                momentTemplates<T,Descriptor>::get_rhoBar_j ( param1.cell ( iX,iY ), rhoBar1, j1 );
                T density1 = Descriptor<T>::fullRho ( rhoBar1 );

                T rhoBar2;
                Array<T,2> j2;
                momentTemplates<T,Descriptor>::get_rhoBar_j ( param2.cell ( iX,iY ), rhoBar2, j2 );
                T density2 = Descriptor<T>::fullRho ( rhoBar2 );

                velocity = 0.5 * ( j1/density1 + j2/density2 );
            }

            // Compute the repelling force.

            Array<T,2> force1 ( ( T ) 0.0, ( T ) 0.0, ( T ) 0.0 );

            if ( isWet ( param1.flag ( iX,iY ) ) && isWet ( param2.flag ( iX,iY ) ) )
            {
                // Remember the envelope is 3.
                for ( plint i=-2; i<3; i++ )
                {
                    plint I = iX + i;
                    for ( plint j=-2; j<3; j++ )
                    {
                        plint J = iY + j;
                        if ( isWet ( param1.flag ( I,J ) ) && isWet ( param2.flag ( I,J ) ) )
                        {
                            T VFdiff = param1.volumeFraction ( I,J ) + param2.volumeFraction ( I,J ) - 1.0;
                            if ( VFdiff> ( T ) 0.0 )
                            {
                                Array<T,2> outwardNormalToFluid1 = param1.getNormal ( I,J ) - param2.getNormal ( I,J );
                                outwardNormalToFluid1 /= norm ( outwardNormalToFluid1 );
                                T r = sqrt ( i*i + j*j );
                                T delta = deltaFunction ( r, ( T ) 3 ); // h is 3 because envelope is 3.
                                force1 += -VFdiff*delta*outwardNormalToFluid1;
                            }
                        }
                    }
                }
                force1 *= interactionStrength;
            }

            if ( param1.flag ( iX,iY ) ==interface && isWet ( param2.flag ( iX,iY ) ) )
            {
                // Compute the average pressure from fluid cells of both fluids.
                T averageDensity = T();
                plint numNeighbors = 0;
                for ( plint i=-1; i<2; i++ )
                {
                    for ( plint j=-1; j<2; j++ )
                    {
                        for ( plint k=-1; k<2; k++ )
                        {
                            if ( isFullWet ( param1.flag ( iX+i,iY+j ) ) )
                            {
                                ++numNeighbors;
                                averageDensity += param1.getDensity ( iX+i,iY+j );
                            }
                            if ( isFullWet ( param2.flag ( iX+i,iY+j ) ) )
                            {
                                ++numNeighbors;
                                averageDensity += convertRho2to1 * param2.getDensity ( iX+i,iY+j );
                            }
                        }
                    }
                }
                if ( numNeighbors==0 )
                {
                    averageDensity = rhoDefault1;
                }
                else
                {
                    averageDensity /= ( T ) numNeighbors;
                }
                param1.setDensity ( iX,iY, averageDensity );
                param1.volumeFraction ( iX,iY ) = param1.mass ( iX,iY ) / averageDensity;
                if ( param2.flag ( iX,iY ) == interface )
                {
                    param1.setMomentum ( iX,iY, param1.getDensity ( iX,iY ) * velocity );
                }
            }

            if ( param2.flag ( iX,iY ) ==interface && isWet ( param1.flag ( iX,iY ) ) )
            {
                // Compute the average pressure from fluid cells of both fluids.
                T averageDensity = T();
                plint numNeighbors = 0;
                for ( plint i=-1; i<2; i++ )
                {
                    for ( plint j=-1; j<2; j++ )
                    {
                        for ( plint k=-1; k<2; k++ )
                        {
                            if ( isFullWet ( param1.flag ( iX+i,iY+j ) ) )
                            {
                                ++numNeighbors;
                                averageDensity += convertRho1to2 * param1.getDensity ( iX+i,iY+j );
                            }
                            if ( isFullWet ( param2.flag ( iX+i,iY+j ) ) )
                            {
                                ++numNeighbors;
                                averageDensity += param2.getDensity ( iX+i,iY+j );
                            }
                        }
                    }
                }
                if ( numNeighbors==0 )
                {
                    averageDensity = rhoDefault2;
                }
                else
                {
                    averageDensity /= ( T ) numNeighbors;
                }
                param2.setDensity ( iX,iY, averageDensity );
                param2.volumeFraction ( iX,iY ) = param2.mass ( iX,iY ) / averageDensity;
                if ( param1.flag ( iX,iY ) == interface )
                {
                    param2.setMomentum ( iX,iY, param2.getDensity ( iX,iY ) * velocity );
                }
            }

            // Implement the repelling force in terms of a momentum correction. Remember that fluids can have strong
            // density differences. In order to have the units of a momentum, this correction term must be
            // multiplied by the density.

            if ( isWet ( param1.flag ( iX,iY ) ) && isWet ( param2.flag ( iX,iY ) ) )
            {
                // Add a force to both fluids. Remember the envelope is 3.

                // force1 is an acceleration.
                T tau1 = T ( 1 ) /param1.cell ( iX,iY ).getDynamics().getOmega();
                param1.setMomentum ( iX,iY, param1.getMomentum ( iX,iY ) + rhoDefault1*tau1*force1 );
                T tau2 = T ( 1 ) /param2.cell ( iX,iY ).getDynamics().getOmega();
                param2.setMomentum ( iX,iY, param2.getMomentum ( iX,iY ) + rhoDefault2*tau2* ( -force1 ) );

                /*
                // force1 is a force.
                T tau1 = T(1)/param1.cell(iX,iY).getDynamics().getOmega();
                param1.setMomentum(iX,iY, param1.getMomentum(iX,iY) + tau1*force1);
                T tau2 = T(1)/param2.cell(iX,iY).getDynamics().getOmega();
                param2.setMomentum(iX,iY, param2.getMomentum(iX,iY) + tau2*(-force1));
                */
            }
        }
    }
}

}  // namespace plb

#endif  // MULTI_FREE_SURFACE_MODEL_2D_HH
