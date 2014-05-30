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

#ifndef FREE_SURFACE_MODEL_2D_HH
#define FREE_SURFACE_MODEL_2D_HH

#include "core/globalDefs.h"
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

/* *************** Class TwoPhaseComputeNormals2D ******************************************* */

template< typename T,template<typename U> class Descriptor>
void TwoPhaseComputeNormals2D<T,Descriptor>::processGenericBlocks (
    Box2D domain,std::vector<AtomicBlock2D*> atomicBlocks )
{
    using namespace twoPhaseFlag;
    FreeSurfaceProcessorParam2D<T,Descriptor> param ( atomicBlocks );

    // Smooth the volume fraction twice. (At the end include also a 1-cell layer around "domain".)
    plint nx = domain.getNx() + 4;
    plint ny = domain.getNy() + 4;
    ScalarField2D<T> smoothVolumeFractionTmp ( nx, ny );
    for ( plint iX=domain.x0-2; iX<=domain.x1+2; ++iX )
    {
        plint i = iX - domain.x0 + 2;
        for ( plint iY=domain.y0-2; iY<=domain.y1+2; ++iY )
        {
            plint j = iY - domain.y0 + 2;
            smoothVolumeFractionTmp.get ( i, j ) = param.smoothVolumeFraction ( iX, iY );
        }
    }

    nx = domain.getNx() + 2;
    ny = domain.getNy() + 2;
    ScalarField2D<T> smoothVolumeFraction ( nx, ny );
    for ( plint iX=domain.x0-1; iX<=domain.x1+1; ++iX )
    {
        plint i    = iX - domain.x0 + 1;
        plint iTmp = iX - domain.x0 + 2;
        for ( plint iY=domain.y0-1; iY<=domain.y1+1; ++iY )
        {
            plint j    = iY - domain.y0 + 1;
            plint jTmp = iY - domain.y0 + 2;

            if ( param.flag ( iX, iY ) == wall )
            {
                smoothVolumeFraction.get ( i, j ) = smoothVolumeFractionTmp.get ( iTmp, jTmp );
                continue;
            }

            T val = 0.0;
            int n = 0;
            for ( int dx = -1; dx < 2; dx++ )
            {
                plint nextX    = iX   + dx;
                plint nextXTmp = iTmp + dx;
                for ( int dy = -1; dy < 2; dy++ )
                {
                    plint nextY    = iY   + dy;
                    plint nextYTmp = jTmp + dy;
                    if ( ! ( dx == 0 && dy == 0 ) && param.flag ( nextX, nextY ) != wall )
                    {
                        n++;
                        val += smoothVolumeFractionTmp.get ( nextXTmp, nextYTmp );
                    }
                }
                if ( n != 0 )
                {
                    val /= ( T ) n;
                }
                else
                {
                    val = smoothVolumeFractionTmp.get ( iTmp, jTmp );
                }

                smoothVolumeFraction.get ( i, j ) = val;
            }
        }
    }

    T eps = getEpsilon<T> ( precision );

    // The outward pointing unit normal is: n = - grad(VOF) / ||grad(VOF)||.

    typedef Descriptor<T> D;
    for ( plint iX=domain.x0; iX<=domain.x1; ++iX )
    {
        for ( plint iY=domain.y0; iY<=domain.y1; ++iY )
        {
            Array<T,2> normal ( ( T ) 0.0, ( T ) 0.0 );

            if ( param.flag ( iX, iY ) == wall )
            {
                param.setNormal ( iX, iY, normal );
                continue;
            }

            int useLB = 1;
            for ( plint iPop = 1; iPop < D::q; ++iPop )
            {
                plint nextX = iX + D::c[iPop][0];
                plint nextY = iY + D::c[iPop][1];
                if ( param.flag ( nextX, nextY ) == wall )
                {
                    useLB = 0;
                    break;
                }
            }

            if ( useLB )
            {
                // Compute the gradient of the smoothed volume fraction "the lattice Boltzmann way".
                // With this method wall cells with a not well defined volume fraction cannot exist
                // at the neighborhood of the point under consideration. This is because wall
                // cells have to be excluded, and there exists no asymmetric lattice-Boltzmann
                // differencing scheme. One-sided first order finite differences have to be used
                // instead.
                for ( plint iPop = 1; iPop < D::q; ++iPop )
                {
                    plint nextX = iX + D::c[iPop][0];
                    plint nextY = iY + D::c[iPop][1];

                    plint i = nextX - domain.x0 + 1;
                    plint j = nextY - domain.y0 + 1;

                    T svf = smoothVolumeFraction.get ( i, j );

                    normal[0] += D::t[iPop]*D::c[iPop][0]*svf;
                    normal[1] += D::t[iPop]*D::c[iPop][1]*svf;
                }
                normal *= D::invCs2;
                T nn = norm ( normal );
                if ( nn <= eps )
                {
                    normal = Array<T,2> ( ( T ) 0.0, ( T ) 0.0 );
                }
                else
                {
                    normal /= -nn;
                }
            }
            else
            {
                // Compute the gradient of the smoothed volume fraction with finite differences
                // excluding the wall cells.
                int fx1 = param.flag ( iX - 1, iY );
                int fx2 = param.flag ( iX + 1, iY );

                int fy1 = param.flag ( iX, iY - 1 );
                int fy2 = param.flag ( iX, iY + 1 );

                plint i, j;
                T h;
                T v1, v2;

                i = iX - domain.x0 + 1;
                j = iY - domain.y0 + 1;

                h = ( fx1 == wall || fx2 == wall ) ? ( T ) 1.0 : ( T ) 2.0;

                v1 = ( fx1 == wall ) ? smoothVolumeFraction.get ( i, j ) :
                     smoothVolumeFraction.get ( i - 1, j );
                v2 = ( fx2 == wall ) ? smoothVolumeFraction.get ( i, j ) :
                     smoothVolumeFraction.get ( i + 1, j );

                normal[0] = ( v2 - v1 ) / h;

                h = ( fy1 == wall || fy2 == wall ) ? ( T ) 1.0 : ( T ) 2.0;

                v1 = ( fy1 == wall ) ? smoothVolumeFraction.get ( i, j ) :
                     smoothVolumeFraction.get ( i, j - 1 );
                v2 = ( fy2 == wall ) ? smoothVolumeFraction.get ( i, j ) :
                     smoothVolumeFraction.get ( i, j + 1 );

                normal[1] = ( v2 - v1 ) / h;

                T nn = norm ( normal );
                if ( nn <= eps )
                {
                    normal = Array<T,2> ( ( T ) 0.0, ( T ) 0.0 );
                }
                else
                {
                    normal /= -nn;
                }
            }

            param.setNormal ( iX, iY, normal );
        }
    }
}

/* *************** Class FreeSurfaceGeometry2D ******************************** */

template<typename T,template<typename U> class Descriptor>
ScalarField2D<int> *FreeSurfaceGeometry2D<T,Descriptor>::getInterfaceFlags ( Box2D domain,
        FreeSurfaceProcessorParam2D<T,Descriptor>& param )
{
    using namespace twoPhaseFlag;

    // Define a temporary scalar field for local use in this function. This scalar field will contain 1 extra
    // layer of cells around "domain".
    plint nx = domain.x1 - domain.x0 + 1;
    plint ny = domain.y1 - domain.y0 + 1;
    ScalarField2D<int> *tmp = new ScalarField2D<int> ( nx + 2, ny + 2, ( int ) unTagged );
    PLB_ASSERT ( tmp );

    // First tag all regular and contact line interface cells. (Loop along 1 envelope cell as well).
    // All interface tags are stored in the temporary storage.
    for ( plint iX=domain.x0-1; iX<=domain.x1+1; ++iX )
    {
        plint indX = iX - domain.x0 + 1;
        for ( plint iY=domain.y0-1; iY<=domain.y1+1; ++iY )
        {
            plint indY = iY - domain.y0 + 1;
            {
                if ( param.flag ( iX, iY ) != interface )
                {
                    tmp->get ( indX, indY ) = notInterface;
                    continue;
                }

                // Find all wall neighbors and store their indices.
                int numWallNeighbors = 0;
                std::vector<Array<plint,2> > wallNeighborIndex;
                for ( int dx = -1; dx < 2; dx++ )
                {
                    plint i = iX + dx;
                    for ( int dy = -1; dy < 2; dy++ )
                    {
                        plint j = iY + dy;
                        if ( ! ( dx == 0 && dy == 0 ) )
                        {
                            if ( param.flag ( i, j ) == wall )
                            {
                                numWallNeighbors++;
                                wallNeighborIndex.push_back ( Array<plint,2> ( i, j ) );
                            }
                        }
                    }
                }

                if ( numWallNeighbors == 0 )
                {
                    tmp->get ( indX, indY ) = regular;
                    continue;
                }

                for ( int dx = -1; dx < 2; dx++ )
                {
                    plint i = iX + dx;
                    for ( int dy = -1; dy < 2; dy++ )
                    {
                        plint j = iY + dy;
                        if ( ! ( dx == 0 && dy == 0 ) )
                        {
                            if ( ( contactAngle >  90.0 && isFullWet ( param.flag ( i, j ) ) ) ||
                                    ( contactAngle <= 90.0 && isEmpty ( param.flag ( i, j ) ) ) )
                            {
                                for ( int dxx = -1; dxx < 2; dxx++ )
                                {
                                    plint ii = i + dxx;
                                    for ( int dyy = -1; dyy < 2; dyy++ )
                                    {
                                        plint jj = j + dyy;
                                        if ( ! ( dxx == 0 && dyy == 0 ) )
                                        {
                                            if ( param.flag ( ii, jj ) == wall )
                                            {
                                                for ( int iWall = 0; iWall < numWallNeighbors; iWall++ )
                                                {
                                                    if ( ii == wallNeighborIndex[iWall][0] &&
                                                            jj == wallNeighborIndex[iWall][1] )
                                                    {
                                                        tmp->get ( indX, indY ) = contactLine;
                                                        goto label0;
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            label0:
                continue;
            }
        }
    }

// Define a scalar field with the interface flags that will be returned from this function.
    ScalarField2D<int> *interfaceFlag = new ScalarField2D<int> ( nx, ny, ( int ) unTagged );
    PLB_ASSERT ( interfaceFlag );

// Now tag all adjacent interface cells and copy all information to the scalar field to be returned.
// At this point all cells that have the flag "unTagged" are non-contact-line interface cells with wall neighbors,
// so they are either "regular" or "adjacent" interface cells.
    for ( plint iX=domain.x0; iX<=domain.x1; ++iX )
    {
        plint indXtmp = iX - domain.x0 + 1;
        plint indX    = iX - domain.x0;
        for ( plint iY=domain.y0; iY<=domain.y1; ++iY )
        {
            plint indYtmp = iY - domain.y0 + 1;
            plint indY    = iY - domain.y0;
            if ( tmp->get ( indXtmp, indYtmp ) != unTagged )
            {
                interfaceFlag->get ( indX, indY ) = tmp->get ( indXtmp, indYtmp );
            }
            else
            {
                int isAdjacent = 0;
                for ( int dx = -1; dx < 2; dx++ )
                {
                    plint i = indXtmp + dx;
                    for ( int dy = -1; dy < 2; dy++ )
                    {
                        plint j = indYtmp + dy;
                        if ( ! ( dx == 0 && dy == 0 ) )
                        {
                            if ( tmp->get ( i, j ) == contactLine )
                            {
                                isAdjacent = 1;
                                interfaceFlag->get ( indX, indY ) = adjacent;
                                goto label1;
                            }
                        }
                    }
                }
            label1:
                if ( !isAdjacent )
                {
                    interfaceFlag->get ( indX, indY ) = regular;
                }
            }
        }
    }

// Check for untagged cells
#ifdef PLB_DEBUG
    for ( plint i = 0; i < nx; i++ )
    {
        for ( plint j = 0; j < ny; j++ )
        {
            PLB_ASSERT ( interfaceFlag->get ( i, j ) != unTagged );
        }
    }
#endif

    delete tmp;

    return interfaceFlag;
}

template<typename T,template<typename U> class Descriptor>
void FreeSurfaceGeometry2D<T,Descriptor>::computeHeights2D ( FreeSurfaceProcessorParam2D<T,Descriptor>& param,
        int integrationDirection, plint iX, plint iY, T h[3][3] )
{
    using namespace twoPhaseFlag;

    // Compute the vector parallel to the integration direction.
    Array<int,2> integrationVector;
    integrationVector[0] = integrationDirection == 0 ? 1 : 0;
    integrationVector[1] = integrationDirection == 1 ? 1 : 0;

    // Compute the vectors tangent to the plane which is normal to the integration vector.
    int iTangentDirection0 = integrationDirection == 0 ? 1 : ( integrationDirection == 1 ) ? 2 : 0;
    int iTangentDirection1 = integrationDirection == 0 ? 2 : ( integrationDirection == 1 ) ? 0 : 1;
    Array<int,2> tangent0;
    tangent0[0] = iTangentDirection0 == 0 ? 1 : 0;
    tangent0[1] = iTangentDirection0 == 1 ? 1 : 0;
    Array<int,2> tangent1;
    tangent1[0] = iTangentDirection1 == 0 ? 1 : 0;
    tangent1[1] = iTangentDirection1 == 1 ? 1 : 0;

    // Calculate the integration stencil width.
    int maxLim = 3;
    for ( int d0 = -1; d0 <= 1; d0++ )
    {
        for ( int d1 = -1; d1 <= 1; d1++ )
        {
            plint posX = iX + d0 * tangent0[0] + d1 * tangent1[0];
            plint posY = iY + d0 * tangent0[1] + d1 * tangent1[1];
            if ( param.flag ( posX, posY ) == wall )
            {
                continue;
            }
            for ( int d = 1; d <= maxLim; d++ )
            {
                plint nextX = posX + d * integrationVector[0];
                plint nextY = posY + d * integrationVector[1];
                if ( param.flag ( nextX, nextY ) == wall )
                {
                    maxLim = std::min ( maxLim, d - 1 );
                    break;
                }
            }
        }
    }

    int minLim = 3;
    for ( int d0 = -1; d0 <= 1; d0++ )
    {
        for ( int d1 = -1; d1 <= 1; d1++ )
        {
            plint posX = iX + d0 * tangent0[0] + d1 * tangent1[0];
            plint posY = iY + d0 * tangent0[1] + d1 * tangent1[1];
            if ( param.flag ( posX, posY ) == wall )
            {
                continue;
            }
            for ( int d = 1; d <= minLim; d++ )
            {
                plint nextX = posX - d * integrationVector[0];
                plint nextY = posY - d * integrationVector[1];
                if ( param.flag ( nextX, nextY ) == wall )
                {
                    minLim = std::min ( minLim, d - 1 );
                    break;
                }
            }
        }
    }

    // Properly initialize heights to -1.
    for ( int i = 0; i < 3; i++ )
        for ( int j = 0; j < 3; j++ )
            h[i][j] = -1.0;

    // Integrate.
    for ( int d0 = -1; d0 <= 1; d0++ )
    {
        int i = d0 + 1;
        for ( int d1 = -1; d1 <= 1; d1++ )
        {
            int j = d1 + 1;

            plint posX = iX + d0 * tangent0[0] + d1 * tangent1[0];
            plint posY = iY + d0 * tangent0[1] + d1 * tangent1[1];
            if ( param.flag ( posX, posY ) == wall )
            {
                continue;
            }
            h[i][j] = 0.0;
            for ( int d = -minLim; d <= maxLim; d++ )
            {
                plint nextX = posX + d * integrationVector[0];
                plint nextY = posY + d * integrationVector[1];
                h[i][j] += param.volumeFraction ( nextX, nextY );
            }
        }
    }

    // Extrapolate on walls. (No contact angle algorithm).
    for ( int i = 0; i < 3; i++ )
    {
        for ( int j = 0; j < 3; j++ )
        {
            if ( fabs ( h[i][j] + 1.0 ) <= eps )
            {
                h[i][j] = h[1][1];
            }
        }
    }
}

template<typename T,template<typename U> class Descriptor>
void FreeSurfaceGeometry2D<T,Descriptor>::computeHeights2D ( FreeSurfaceProcessorParam2D<T,Descriptor>& param,
        Array<int,2>& wallTangent0, Array<int,2>& wallTangent1, int integrationDirection, plint iX, plint iY,T h[3] )
{
    using namespace twoPhaseFlag;

    // Compute the vector parallel to the integration direction.
    Array<int,2> integrationVector;
    integrationVector[0] = integrationDirection == 0 ? wallTangent0[0] : wallTangent1[0];
    integrationVector[1] = integrationDirection == 0 ? wallTangent0[1] : wallTangent1[1];


    // Compute the vector tangent to the line which is normal to the integration vector.
    Array<int,2> tangent;
    tangent[0] = integrationDirection == 0 ? wallTangent1[0] : wallTangent0[0];
    tangent[1] = integrationDirection == 0 ? wallTangent1[1] : wallTangent0[1];


    // Calculate the integration stencil width.
    int maxLim = 3;
    for ( int d0 = -1; d0 <= 1; d0++ )
    {
        plint posX = iX + d0 * tangent[0];
        plint posY = iY + d0 * tangent[1];

        if ( param.flag ( posX, posY ) == wall )
        {
            continue;
        }
        for ( int d = 1; d <= maxLim; d++ )
        {
            plint nextX = posX + d * integrationVector[0];
            plint nextY = posY + d * integrationVector[1];

            if ( param.flag ( nextX, nextY ) == wall )
            {
                maxLim = std::min ( maxLim, d - 1 );
                break;
            }
        }
    }

    int minLim = 3;
    for ( int d0 = -1; d0 <= 1; d0++ )
    {
        plint posX = iX + d0 * tangent[0];
        plint posY = iY + d0 * tangent[1];

        if ( param.flag ( posX, posY ) == wall )
        {
            continue;
        }
        for ( int d = 1; d <= minLim; d++ )
        {
            plint nextX = posX - d * integrationVector[0];
            plint nextY = posY - d * integrationVector[1];

            if ( param.flag ( nextX, nextY ) == wall )
            {
                minLim = std::min ( minLim, d - 1 );
                break;
            }
        }
    }

    // Properly initialize heights to -1.
    h[0] = h[1] = h[2] = -1.0;

    // Integrate.
    for ( int d0 = -1; d0 <= 1; d0++ )
    {
        int i = d0 + 1;
        plint posX = iX + d0 * tangent[0];
        plint posY = iY + d0 * tangent[1];
        if ( param.flag ( posX, posY ) == wall )
        {
            continue;
        }
        h[i] = 0.0;
        for ( int d = -minLim; d <= maxLim; d++ )
        {
            plint nextX = posX + d * integrationVector[0];
            plint nextY = posY + d * integrationVector[1];
            h[i] += param.volumeFraction ( nextX, nextY );
        }
    }

    // Extrapolate on walls. (No contact angle algorithm).
    for ( int i = 0; i < 3; i++ )
    {
        if ( fabs ( h[i] + 1.0 ) <= eps )
        {
            h[i] = h[1];
        }
    }
}

template<typename T,template<typename U> class Descriptor>
void FreeSurfaceGeometry2D<T,Descriptor>::processGenericBlocks ( Box2D domain, std::vector<AtomicBlock2D*> atomicBlocks )
{
    using namespace twoPhaseFlag;
    FreeSurfaceProcessorParam2D<T,Descriptor> param ( atomicBlocks );

    Array<T,2> zeroVector ( ( T ) 0, ( T ) 0, ( T ) 0 );

    if ( !useContactAngle )
    {
        for ( plint iX=domain.x0; iX<=domain.x1; ++iX )
        {
            for ( plint iY=domain.y0; iY<=domain.y1; ++iY )
                if ( param.flag ( iX, iY ) == interface )
                {
                    if ( param.isBoundary ( iX, iY ) )
                    {
                        param.curvature ( iX, iY ) = 0.0;
                        param.setNormal ( iX, iY, zeroVector );
                        continue;
                    }
                    // Locally smooth the volume fraction to compute an estimate of the normal.
                    T svfcp = param.smoothVolumeFraction ( iX, iY );
                    T svfx0 = param.flag ( iX - 1, iY ) != wall ? param.smoothVolumeFraction ( iX - 1, iY ) : svfcp;
                    T svfx1 = param.flag ( iX + 1, iY ) != wall ? param.smoothVolumeFraction ( iX + 1, iY ) : svfcp;
                    T svfy0 = param.flag ( iX, iY - 1 ) != wall ? param.smoothVolumeFraction ( iX, iY - 1 ) : svfcp;
                    T svfy1 = param.flag ( iX, iY + 1 ) != wall ? param.smoothVolumeFraction ( iX, iY + 1 ) : svfcp;

                    // Compute a normalized grad(VF) (inward-pointing normal).
                    Array<T,2> gradVF;
                    gradVF[0] = 0.5 * ( svfx1 - svfx0 );
                    gradVF[1] = 0.5 * ( svfy1 - svfy0 );
                    T norm_gradVF = norm ( gradVF );
                    if ( norm_gradVF <= eps )
                    {
                        param.curvature ( iX, iY ) = 0.0;
                        param.setNormal ( iX, iY, zeroVector );
                        continue;
                    }
                    gradVF /= norm_gradVF;

                    T abs0 = fabs ( gradVF[0] );
                    T abs1 = fabs ( gradVF[1] );

                    int integrationDirection=2;

                    if ( abs0 > abs1 )
                        integrationDirection = 0;
                    else
                        integrationDirection = 1;


                    T h[3][3];
                    computeHeights2D ( param, integrationDirection, iX, iY, h );

                    T dh0 = 0.5 * ( h[2][1] - h[0][1] );
                    T dh1 = 0.5 * ( h[1][2] - h[1][0] );

                    T dh00 = h[2][1] - 2.0 * h[1][1] + h[0][1];
                    T dh11 = h[1][2] - 2.0 * h[1][1] + h[1][0];

                    T dh01 = 0.25 * ( h[2][2] - h[2][0] - h[0][2] + h[0][0] );

                    T value = - ( dh00 + dh11 + dh00 * dh1 * dh1 + dh11 * dh0 * dh0 - 2.0 * dh01 * dh0 * dh1 ) /
                              pow ( 1.0 + dh0 * dh0 + dh1 * dh1, 1.5 );

                    param.curvature ( iX, iY ) = value;

                    T sgn = -gradVF[integrationDirection] < 0.0 ? -1.0 : 1.0;
                    Array<T,2> normal;
                    if ( integrationDirection == 0 )
                    {
                        normal = Array<T,2> ( sgn, -dh0, -dh1 );
                    }
                    else if ( integrationDirection == 1 )
                    {
                        normal = Array<T,2> ( -dh1, sgn, -dh0 );
                    }
                    else
                    {
                        normal = Array<T,2> ( -dh0, -dh1, sgn );
                    }
                    T norm_normal = norm ( normal );
                    if ( norm_normal <= eps )
                    {
                        param.setNormal ( iX, iY, zeroVector );
                    }
                    else
                    {
                        param.setNormal ( iX, iY, normal / norm_normal );
                    }
                }
                else
                {
                    param.curvature ( iX, iY ) = 0.0;
                    param.setNormal ( iX, iY, zeroVector );
                }
        }
    }
    else     // Use contact angles.
    {
        // First compute the flags.
        ScalarField2D<int> *interfaceFlag = getInterfaceFlags ( domain, param );

        /* New contact angle algorithm. This algorithm still does not properly treat the adjacent cells. */

        // First loop over all the regular and adjacent interface cells and calculate the curvature and the normal vectors.
        // When the appropriate algorithm for the adjacent cells is implemented, they must be removed from these loops.
        for ( plint iX=domain.x0; iX<=domain.x1; ++iX )
        {
            plint i = iX - domain.x0;
            for ( plint iY=domain.y0; iY<=domain.y1; ++iY )
            {
                plint j = iY - domain.y0;

                if ( interfaceFlag->get ( i, j ) == regular || interfaceFlag->get ( i, j ) == adjacent )
                {
                    if ( param.isBoundary ( iX, iY ) )
                    {
                        param.curvature ( iX, iY ) = 0.0;
                        param.setNormal ( iX, iY, zeroVector );
                        continue;
                    }
                    // Locally smooth the volume fraction to compute an estimate of the normal.
                    T svfcp = param.smoothVolumeFraction ( iX, iY );
                    T svfx0 = param.flag ( iX - 1, iY ) != wall ? param.smoothVolumeFraction ( iX - 1, iY ) : svfcp;
                    T svfx1 = param.flag ( iX + 1, iY ) != wall ? param.smoothVolumeFraction ( iX + 1, iY ) : svfcp;
                    T svfy0 = param.flag ( iX, iY - 1 ) != wall ? param.smoothVolumeFraction ( iX, iY - 1 ) : svfcp;
                    T svfy1 = param.flag ( iX, iY + 1 ) != wall ? param.smoothVolumeFraction ( iX, iY + 1 ) : svfcp;

                    // Compute a normalized grad(VF) (inward-pointing normal).
                    Array<T,2> gradVF;
                    gradVF[0] = 0.5 * ( svfx1 - svfx0 );
                    gradVF[1] = 0.5 * ( svfy1 - svfy0 );
                    T norm_gradVF = norm ( gradVF );
                    if ( norm_gradVF <= eps )
                    {
                        param.curvature ( iX, iY ) = 0.0;
                        param.setNormal ( iX, iY, zeroVector );
                        continue;
                    }
                    gradVF /= norm_gradVF;

                    T abs0 = fabs ( gradVF[0] );
                    T abs1 = fabs ( gradVF[1] );

                    int integrationDirection=2;

                    if ( abs0 > abs1 )
                        integrationDirection = 0;
                    else
                        integrationDirection = 1;


                    T h[3][3];
                    computeHeights2D ( param, integrationDirection, iX, iY, h );

                    T dh0 = 0.5 * ( h[2][1] - h[0][1] );
                    T dh1 = 0.5 * ( h[1][2] - h[1][0] );

                    T dh00 = h[2][1] - 2.0 * h[1][1] + h[0][1];
                    T dh11 = h[1][2] - 2.0 * h[1][1] + h[1][0];

                    T dh01 = 0.25 * ( h[2][2] - h[2][0] - h[0][2] + h[0][0] );

                    T value = - ( dh00 + dh11 + dh00 * dh1 * dh1 + dh11 * dh0 * dh0 - 2.0 * dh01 * dh0 * dh1 ) /
                              pow ( 1.0 + dh0 * dh0 + dh1 * dh1, 1.5 );

                    param.curvature ( iX, iY ) = value;

                    T sgn = -gradVF[integrationDirection] < 0.0 ? -1.0 : 1.0;
                    Array<T,2> normal;
                    if ( integrationDirection == 0 )
                    {
                        normal = Array<T,2> ( sgn, -dh0, -dh1 );
                    }
                    else if ( integrationDirection == 1 )
                    {
                        normal = Array<T,2> ( -dh1, sgn, -dh0 );
                    }
                    else
                    {
                        normal = Array<T,2> ( -dh0, -dh1, sgn );
                    }
                    T norm_normal = norm ( normal );
                    if ( norm_normal <= eps )
                    {
                        param.setNormal ( iX, iY, zeroVector );
                    }
                    else
                    {
                        param.setNormal ( iX, iY, normal / norm_normal );
                    }
                }
                else
                {
                    param.curvature ( iX, iY ) = 0.0;
                    param.setNormal ( iX, iY, zeroVector );
                }
            }
        }

        // Then loop over all the contact-line interface cells and calculate the curvature and the
        // normal vectors according to the specified contact angle.
        T pi = acos ( -1.0 );
        T contactAngleRad = contactAngle * pi / 180.0;
        T tanContactAngle = tan ( contactAngleRad );

        for ( plint iX=domain.x0; iX<=domain.x1; ++iX )
        {
            plint i = iX - domain.x0;

            for ( plint iY=domain.y0; iY<=domain.y1; ++iY )
            {
                plint j = iY - domain.y0;
                if ( interfaceFlag->get ( i, j ) == contactLine )
                {
                    if ( param.isBoundary ( iX, iY ) )
                    {
                        param.curvature ( iX, iY ) = 0.0;
                        param.setNormal ( iX, iY, zeroVector );
                        continue;
                    }

                    // First decide where is the wall.
                    int numWallCells = 0;
                    // Computation of the inward-pointing wall normal (not unitary).
                    Array<int,2> inwardWallNormal ( 0, 0 );
                    for ( int dx = -1; dx < 2; dx++ )
                    {
                        for ( int dy = -1; dy < 2; dy++ )
                        {
                            if ( param.flag ( iX + dx, iY + dy ) == wall )
                            {
                                inwardWallNormal += Array<int,2> ( -dx, -dy );
                                numWallCells++;
                            }
                        }
                    }
                    PLB_ASSERT ( numWallCells != 0 );
                    int norm2_inwardWallNormal = inwardWallNormal[0] * inwardWallNormal[0] +
                                                 inwardWallNormal[1] * inwardWallNormal[1];
                    PLB_ASSERT ( norm2_inwardWallNormal != 0 );

                    int iWallNormalDirection;
                    // The inwardWallNormal is aligned with one axis.
                    if ( inwardWallNormal[0] != 0 && inwardWallNormal[1] == 0 )
                    {
                        iWallNormalDirection = 0;
                    }
                    else if ( inwardWallNormal[1] != 0 && inwardWallNormal[0] == 0 )
                    {
                        iWallNormalDirection = 1;
                    }
                    else if ( inwardWallNormal[0] == 0 && inwardWallNormal[1] == 0 )
                    {
                        iWallNormalDirection = 2;
                    }
                    else
                    {
                        // The inwardWallNormal is not aligned with one axis.
                        Array<int,2> sumDirection[3];
                        sumDirection[0] = inwardWallNormal[0] == 0 ? Array<int,2> ( 0,  0,  0 ) :
                                          ( inwardWallNormal[0] >  0 ? Array<int,2> ( 1,  0,  0 ) :
                                            Array<int,2> ( -1,  0,  0 ) );
                        sumDirection[1] = inwardWallNormal[1] == 0 ? Array<int,2> ( 0,  0,  0 ) :
                                          ( inwardWallNormal[1] >  0 ? Array<int,2> ( 0,  1,  0 ) :
                                            Array<int,2> ( 0, -1,  0 ) );
                        sumDirection[2] = inwardWallNormal[2] == 0 ? Array<int,2> ( 0,  0,  0 ) :
                                          ( inwardWallNormal[2] >  0 ? Array<int,2> ( 0,  0,  1 ) :
                                            Array<int,2> ( 0,  0, -1 ) );
                        T sum[3] = { std::numeric_limits<T>::max(),
                                     std::numeric_limits<T>::max(),
                                     std::numeric_limits<T>::max()
                                   };
                        for ( int iSum = 0; iSum < 3; iSum++ )
                        {
                            if ( sumDirection[iSum][0] + sumDirection[iSum][1] != 0 )
                            {
                                sum[iSum] = 0.0;
                                for ( int d = 0; d <= 3; d++ )
                                {
                                    plint posX = iX + d * sumDirection[iSum][0];
                                    plint posY = iY + d * sumDirection[iSum][1];
                                    if ( param.flag ( posX, posY ) != wall )
                                    {
                                        sum[iSum] += param.volumeFraction ( posX, posY );
                                    }
                                }
                            }
                        }

                        // The wall normal direction is the direction of the smallest sum.
                        if ( sum[0] < sum[1] )
                        {
                            if ( sum[0] < sum[2] )
                            {
                                iWallNormalDirection = 0;
                            }
                            // sum[0]<sum[1] && sum[0] >= sum[2]
                            else
                            {
                                iWallNormalDirection = 2;
                            }
                        }
                        // sum[0] >= sum[1]
                        else
                        {
                            if ( sum[1] < sum[2] )
                            {
                                iWallNormalDirection = 1;
                            }
                            // sum[0] >= sum[1] && sum[1] >= sum[2]
                            else
                            {
                                iWallNormalDirection = 2;
                            }
                        }
                    }

                    // Reset the inward wall normal to be unitary and to contain information on the direction.
                    inwardWallNormal[0] = iWallNormalDirection != 0 ? 0 : ( inwardWallNormal[0] > 0 ? 1 : -1 );
                    inwardWallNormal[1] = iWallNormalDirection != 1 ? 0 : ( inwardWallNormal[1] > 0 ? 1 : -1 );

                    // Define a wall normal that shows only the wall normal axis.
                    Array<int,2> wallNormal;
                    wallNormal[0] = iWallNormalDirection == 0 ? 1 : 0;
                    wallNormal[1] = iWallNormalDirection == 1 ? 1 : 0;

                    // Compute the wall tangent vectors.
                    int iWallTangentDirection0 = iWallNormalDirection == 0 ? 1 : ( iWallNormalDirection == 1 ) ? 2 : 0;
                    int iWallTangentDirection1 = iWallNormalDirection == 0 ? 2 : ( iWallNormalDirection == 1 ) ? 0 : 1;
                    Array<int,2> wallTangent0;
                    wallTangent0[0] = iWallTangentDirection0 == 0 ? 1 : 0;
                    wallTangent0[1] = iWallTangentDirection0 == 1 ? 1 : 0;
                    Array<int,2> wallTangent1;
                    wallTangent1[0] = iWallTangentDirection1 == 0 ? 1 : 0;
                    wallTangent1[1] = iWallTangentDirection1 == 1 ? 1 : 0;

                    // Locally smooth the volume fraction to compute an estimate of the 2D normal.
                    T svfcp = param.smoothVolumeFraction ( iX, iY );
                    plint posX, posY;
                    posX = iX - wallTangent0[0];
                    posY = iY - wallTangent0[1];
                    T svf00 = param.flag ( posX, posY ) != wall ? param.smoothVolumeFraction ( posX, posY ) : svfcp;
                    posX = iX + wallTangent0[0];
                    posY = iY + wallTangent0[1];

                    T svf01 = param.flag ( posX, posY ) != wall ? param.smoothVolumeFraction ( posX, posY ) : svfcp;
                    posX = iX - wallTangent1[0];
                    posY = iY - wallTangent1[1];
                    T svf10 = param.flag ( posX, posY ) != wall ? param.smoothVolumeFraction ( posX, posY ) : svfcp;
                    posX = iX + wallTangent1[0];
                    posY = iY + wallTangent1[1];
                    T svf11 = param.flag ( posX, posY ) != wall ? param.smoothVolumeFraction ( posX, posY ) : svfcp;

                    // Compute a normalized 2D grad(VF) (inward-pointing 2D normal).
                    Array<T,2> gradVF2D;
                    gradVF2D[0] = 0.5 * ( svf01 - svf00 );
                    gradVF2D[1] = 0.5 * ( svf11 - svf10 );
                    T norm_gradVF2D = norm ( gradVF2D );
                    if ( norm_gradVF2D <= eps )
                    {
                        param.curvature ( iX, iY ) = 0.0;
                        param.setNormal ( iX, iY, zeroVector );
                        continue;
                    }
                    gradVF2D /= norm_gradVF2D;

                    T abs02D = fabs ( gradVF2D[0] );
                    T abs12D = fabs ( gradVF2D[1] );

                    int integrationDirection2D = 1; // wallTangent1.
                    if ( abs02D > abs12D )
                    {
                        integrationDirection2D = 0; // wallTangent0.
                    }

                    T h2D[3];
                    computeHeights2D ( param, wallTangent0, wallTangent1, integrationDirection2D, iX, iY, h2D );

                    T dh2D = 0.5 * ( h2D[2] - h2D[0] );

                    T sgn2D = -gradVF2D[integrationDirection2D] < 0.0 ? -1.0 : 1.0;
                    Array<T,2> normal2D;
                    if ( integrationDirection2D == 0 ) // With respect to wallTangent0 and wallTangent1.
                    {
                        normal2D = Array<T,2> ( sgn2D, -dh2D );
                    }
                    else
                    {
                        normal2D = Array<T,2> ( -dh2D, sgn2D );
                    }
                    T norm_normal2D = norm ( normal2D );
                    if ( norm_normal2D <= eps )
                    {
                        param.curvature ( iX, iY ) = 0.0;
                        param.setNormal ( iX, iY, zeroVector );
                        continue;
                    }

                    Array<T,2> normal; // 2D outward unit normal.
                    T wallNormalComponent = norm_normal2D / tanContactAngle;
                    normal[0] = normal2D[0] * wallTangent0[0] + normal2D[1] * wallTangent1[0] + wallNormalComponent * wallNormal[0];
                    normal[1] = normal2D[0] * wallTangent0[1] + normal2D[1] * wallTangent1[1] + wallNormalComponent * wallNormal[1];
                    T norm_normal = norm ( normal );
                    if ( norm_normal <= eps )
                    {
                        param.setNormal ( iX, iY, zeroVector );
                    }
                    else
                    {
                        param.setNormal ( iX, iY, normal / norm_normal );
                    }

                    // Now compute the curvature.
                    // First compute the 2D height functions.
                    int integrationDirection;
                    if ( integrationDirection2D == 0 )
                    {
                        integrationDirection = wallTangent0[0] != 0 ? 0 : ( wallTangent0[1] != 0 ? 1 : 2 );
                    }
                    else
                    {
                        integrationDirection = wallTangent1[0] != 0 ? 0 : ( wallTangent1[1] != 0 ? 1 : 2 );
                    }
                    T h[3][3];
                    computeHeights2D ( param, integrationDirection, iX, iY, h );

                    // Determine the orientation of the elements of h.
                    int iTangentDirection0 = integrationDirection == 0 ? 1 : ( integrationDirection == 1 ) ? 2 : 0;
                    int iTangentDirection1 = integrationDirection == 0 ? 2 : ( integrationDirection == 1 ) ? 0 : 1;
                    Array<int,2> tangent0;
                    tangent0[0] = iTangentDirection0 == 0 ? 1 : 0;
                    tangent0[1] = iTangentDirection0 == 1 ? 1 : 0;
                    Array<int,2> tangent1;
                    tangent1[0] = iTangentDirection1 == 0 ? 1 : 0;
                    tangent1[1] = iTangentDirection1 == 1 ? 1 : 0;

                    int i0 = -1;
                    int j0 = -1;
                    if ( inwardWallNormal[0] == tangent0[0] &&
                            inwardWallNormal[1] == tangent0[1] )
                    {
                        i0 = 0;
                    }
                    else if ( inwardWallNormal[0] == -tangent0[0] &&
                              inwardWallNormal[1] == -tangent0[1] )
                    {
                        i0 = 2;
                    }
                    else if ( inwardWallNormal[0] == tangent1[0] &&
                              inwardWallNormal[1] == tangent1[1] )
                    {
                        j0 = 0;
                    }
                    else if ( inwardWallNormal[0] == -tangent1[0] &&
                              inwardWallNormal[1] == -tangent1[1] )
                    {
                        j0 = 2;
                    }
                    else
                    {
                        PLB_ASSERT ( false );
                    }

                    Array<T,2> v1, v2; // In the wallTangent0, wallTangent1 base.
                    v1[0] = fabs ( normal2D[0] );
                    v1[1] = fabs ( normal2D[1] );
                    if ( integrationDirection2D == 0 )
                    {
                        v2[0] = 1.0;
                        v2[1] = 0.0;
                    }
                    else
                    {
                        v2[0] = 0.0;
                        v2[1] = 1.0;
                    }
                    T cosAlpha = cos ( angleBetweenVectors ( v1, v2 ) );
                    T correction = 1.0 / ( tanContactAngle * cosAlpha );

                    if ( i0 != -1 && j0 == -1 )
                    {
                        for ( int d = 0; d < 3; d++ )
                        {
                            h[i0][d] = h[1][d] + correction;
                        }
                    }
                    else if ( i0 == -1 && j0 != -1 )
                    {
                        for ( int d = 0; d < 3; d++ )
                        {
                            h[d][j0] = h[d][1] + correction;
                        }
                    }
                    else
                    {
                        PLB_ASSERT ( false );
                    }

                    T dh0 = 0.5 * ( h[2][1] - h[0][1] );
                    T dh1 = 0.5 * ( h[1][2] - h[1][0] );

                    T dh00 = h[2][1] - 2.0 * h[1][1] + h[0][1];
                    T dh11 = h[1][2] - 2.0 * h[1][1] + h[1][0];

                    T dh01 = 0.25 * ( h[2][2] - h[2][0] - h[0][2] + h[0][0] );

                    T value = - ( dh00 + dh11 + dh00 * dh1 * dh1 + dh11 * dh0 * dh0 - 2.0 * dh01 * dh0 * dh1 ) /
                              pow ( 1.0 + dh0 * dh0 + dh1 * dh1, 1.5 );

                    param.curvature ( iX, iY ) = value;
                }
            }
        }

        delete interfaceFlag;
    }
}

/* *************** Class TwoPhaseComputeCurvature2D ******************************** */

template<typename T,template<typename U> class Descriptor>
void TwoPhaseComputeCurvature2D<T,Descriptor>::processGenericBlocks ( Box2D domain, std::vector<AtomicBlock2D*> atomicBlocks )
{
    using namespace twoPhaseFlag;
    FreeSurfaceProcessorParam2D<T,Descriptor> param ( atomicBlocks );

    // Tensor field to hold a temporary vector field of unit normals. (Include also a 1-cell layer around "domain".)
    plint nx = domain.getNx() + 2;
    plint ny = domain.getNy() + 2;
    TensorField2D<T,2> tmpNormal ( nx, ny );
    for ( plint iX=domain.x0-1; iX<=domain.x1+1; ++iX )
    {
        plint i = iX - domain.x0 + 1;
        for ( plint iY=domain.y0-1; iY<=domain.y1+1; ++iY )
        {
            plint j = iY - domain.y0 + 1;
            tmpNormal.get ( i, j ) = param.getNormal ( iX, iY );
        }
    }

    // Enforce contact angles.
    if ( useContactAngle )
    {
        Dot2D absOffset = param.absOffset();
        T eps = getEpsilon<T> ( precision );

        for ( plint iX=domain.x0-1; iX<=domain.x1+1; ++iX )
        {
            for ( plint iY=domain.y0-1; iY<=domain.y1+1; ++iY )
            {
                if ( contained ( iX+absOffset.x, iY+absOffset.y, globalBoundingBox ) )
                {
                    if ( param.flag ( iX, iY ) == interface )
                    {
                        int isaContactAngleCell = 0;
                        int numWallCells = 0;
                        int numEmptyCells = 0;
                        // Computation of the inward-pointing wall normal.
                        Array<int,2> tmpWallNormal ( 0, 0 );
                        for ( int i = -1; i < 2; i++ )
                        {
                            for ( int j = -1; j < 2; j++ )
                            {
                                if ( contained ( iX+i+absOffset.x, iY+j+absOffset.y,globalBoundingBox ) )
                                {
                                    int flg = param.flag ( iX+i, iY+j );
                                    if ( flg == wall )
                                    {
                                        tmpWallNormal += Array<int,2> ( -i, -j );
                                        numWallCells++;
                                    }
                                    else if ( isEmpty ( flg ) )
                                    {
                                        numEmptyCells++;
                                    }
                                }
                            }
                        }
                        Array<T,2> wallNormal;
                        if ( numWallCells != 0 && numEmptyCells != 0 )
                        {
                            int norm2tmpWallNormal = tmpWallNormal[0] * tmpWallNormal[0] +
                                                     tmpWallNormal[1] * tmpWallNormal[1];
                            if ( norm2tmpWallNormal != 0 )
                            {
                                T tmpNormWallNormal = ( T ) sqrt ( ( T ) norm2tmpWallNormal );
                                wallNormal[0] = ( T ) tmpWallNormal[0] / tmpNormWallNormal;
                                wallNormal[1] = ( T ) tmpWallNormal[1] / tmpNormWallNormal;
                                isaContactAngleCell = 1;
                            }
                        }
                        if ( isaContactAngleCell )
                        {
                            // Construction of a new orthonormal basis.
//                             Array<T,2> wallTangent0 ( -wallNormal[1], wallNormal[0] );
//                             Array<T,2> wallTangent1 ( ( T ) 0.0, ( T ) 0.0 );
                            //FIXME
//                             gramSchmidt ( wallNormal, wallTangent0, wallTangent1 );

                            // Transformation matrix from the standard Euclidean basis to the basis
                            // (wallTangent0, wallTangent1, wallNormal).
                            T a[2][2];
                            a[0][0] = wallNormal[0];
                            a[0][1] = wallNormal[1];
                            a[1][0] = -wallNormal[1];
                            a[1][1] = wallNormal[0];
//TODO FIX THIS
                            T det = a[0][0]*a[1][1] - a[0][1]*a[1][0];


                            PLB_ASSERT ( fabs ( det ) > eps );

                            // Make sure that the new basis is counter-clockwise oriented.
#if 0
                            if ( det < ( T ) 0.0 )
                            {
                                Array<T,2> tmp ( wallTangent0 );
                                wallTangent0 = wallTangent1;
                                wallTangent1 = tmp;

                                a[0][0] = wallTangent0[0];
                                a[0][1] = wallTangent0[1];
                                a[1][0] = wallTangent1[0];
                                a[1][1] = wallTangent1[1];

                                det = -det;
                            }
#endif
                            T inv_det = 1.0 / fabs ( det );

                            // Inverse of the transformation matrix.
                            T inv[2][2];
                            inv[0][0] = inv_det *a[1][1];
                            inv[0][1] = -inv_det *a[0][1];
                            inv[1][0] = -inv_det *a[1][0];
                            inv[1][1] = inv_det *a[0][0];

                            // Express the outward pointing unit normal of the free surface
                            // in the new basis.
                            Array<T,2> normal = param.getNormal ( iX, iY );
                            Array<T,2> newNormal ( ( T ) 0.0, ( T ) 0.0 );
                            for ( int i = 0; i < 2; i++ )
                            {
                                for ( int j = 0; j < 2; j++ )
                                {
                                    newNormal[i] += normal[j] * inv[j][i];
                                }
                            }

                            // Spherical coordinates.
                            // The contact angle is the angle between the free surface normal vector,
                            // and the wall normal.
                            T phi = atan2 ( newNormal[1], newNormal[0] );
                            T theta = contactAngle; // In radians.

                            newNormal[0] = cos ( phi ) * sin ( theta );
                            newNormal[1] = sin ( phi ) * sin ( theta );

                            normal.resetToZero();
                            for ( int i = 0; i < 2; i++ )
                            {
                                for ( int j = 0; j < 2; j++ )
                                {
                                    normal[i] += newNormal[j] * a[j][i];
                                }
                            }

                            // Enforce the required free surface normal at the interface cell under
                            // consideration, so that the correct contact angle is implicitly imposed.

                            plint i = iX - domain.x0 + 1;
                            plint j = iY - domain.y0 + 1;

                            tmpNormal.get ( i, j ) = normal;
                        }
                    }
                }
            }
        }
    }

    // Compute the curvature as the divergence of the vector field of unit normals.
    typedef Descriptor<T> D;
    for ( plint iX=domain.x0; iX<=domain.x1; ++iX )
    {
        for ( plint iY=domain.y0; iY<=domain.y1; ++iY )
        {
            T curv = 0.0;

            if ( param.flag ( iX, iY ) != interface )
            {
                param.curvature ( iX, iY ) = curv;
                continue;
            }

            int useLB = 1;
            for ( plint iPop = 1; iPop < D::q; ++iPop )
            {
                plint nextX = iX + D::c[iPop][0];
                plint nextY = iY + D::c[iPop][1];

                if ( param.flag ( nextX, nextY ) == wall )
                {
                    useLB = 0;
                    break;
                }
            }

            if ( useLB )
            {
                // Compute the divergence of the normal vector field "the lattice Boltzmann way".
                curv = 0.0;
                for ( plint iPop=1; iPop < D::q; ++iPop )
                {
                    plint nextX = iX + D::c[iPop][0];
                    plint nextY = iY + D::c[iPop][1];

                    plint i = nextX - domain.x0 + 1;
                    plint j = nextY - domain.y0 + 1;

                    Array<T,2>& normal = tmpNormal.get ( i, j );

                    curv += D::t[iPop]* ( D::c[iPop][0]*normal[0] + D::c[iPop][1]*normal[1] );
                }
                curv *= D::invCs2;
            }
            else
            {
                // Compute the divergence with finite differences on the interface cells excluding wall cells.
                int fx1 = param.flag ( iX - 1, iY );
                int fx2 = param.flag ( iX + 1, iY );

                int fy1 = param.flag ( iX, iY - 1 );
                int fy2 = param.flag ( iX, iY + 1 );


                plint i, j;
                T h;
                T dnx_dx, dny_dy;
                T v1, v2;

                i = iX - domain.x0 + 1;
                j = iY - domain.y0 + 1;

                h = ( fx1 == wall || fx2 == wall ) ? ( T ) 1.0 : ( T ) 2.0;

                v1 = ( fx1 == wall ) ? tmpNormal.get ( i, j ) [0] : tmpNormal.get ( i - 1, j ) [0];
                v2 = ( fx2 == wall ) ? tmpNormal.get ( i, j ) [0] : tmpNormal.get ( i + 1, j ) [0];

                dnx_dx = ( v2 - v1 ) / h;

                h = ( fy1 == wall || fy2 == wall ) ? ( T ) 1.0 : ( T ) 2.0;

                v1 = ( fy1 == wall ) ? tmpNormal.get ( i, j ) [1] : tmpNormal.get ( i, j - 1 ) [1];
                v2 = ( fy2 == wall ) ? tmpNormal.get ( i, j ) [1] : tmpNormal.get ( i, j + 1 ) [1];

                dny_dy = ( v2 - v1 ) / h;

                curv = dnx_dx + dny_dy;
            }

            // We restrict the radius of curvature to be more always >=0.5, in lattice units.
            // A smaller radius makes no sense anyway, numerically speaking, and in this way
            // we avoid problems of the "division by zero" kind. (radius = 2/curvature)
            if ( fabs ( curv ) >4.0 )
            {
                if ( curv < 0. )
                {
                    curv = -4.0;
                }
                else
                {
                    curv = 4.0;
                }
            }
            param.curvature ( iX, iY ) = curv;
        }
    }
}

/* *************** Class FreeSurfaceMassChange2D ******************************************* */

template< typename T,template<typename U> class Descriptor>
void FreeSurfaceMassChange2D<T,Descriptor>::processGenericBlocks (
    Box2D domain,std::vector<AtomicBlock2D*> atomicBlocks )
{
    typedef Descriptor<T> D;
    using namespace twoPhaseFlag;
    FreeSurfaceProcessorParam2D<T,Descriptor> param ( atomicBlocks );

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
                freeSurfaceTemplates<T,Descriptor>::massExchangeFluidCell ( param, iX,iY );
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
        }
    }
}

/* *************** Class FreeSurfaceCompletion2D ******************************************* */

template< typename T,template<typename U> class Descriptor>
void FreeSurfaceCompletion2D<T,Descriptor>
::processGenericBlocks ( Box2D domain,std::vector<AtomicBlock2D*> atomicBlocks )
{
    typedef Descriptor<T> D;
    using namespace twoPhaseFlag;
    typedef typename InterfaceLists<T,Descriptor>::Node Node;

    FreeSurfaceProcessorParam2D<T,Descriptor> param ( atomicBlocks );

    // In this data processor, populations are both written locally and read non-locally.
    // To guarantee data consistency, a first loop makes only read accesses and stores
    // the necessary information into the list neighborOppositePop. A second loop reads
    // from this list and assigns values to populations.
    std::map<Node, Array<T,D::q> > neighborOppositePop;
    for ( plint iX=domain.x0; iX<=domain.x1; ++iX )
    {
        for ( plint iY=domain.y0; iY<=domain.y1; ++iY )
        {
            if ( param.flag ( iX,iY ) == interface )
            {
                // Here we are on an interface node. The entire set of fi's is reconstructed.
                bool needsModification = false;
                Array<T,D::q> savedPop;
                savedPop[0] = -2.;
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
                }
                if ( needsModification )
                {
                    neighborOppositePop.insert ( std::pair<Node,Array<T,D::q> > ( Node ( iX,iY ), savedPop ) );
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
                // F, and is equal to the ambient pressure for a
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
}

/* *************** Class FreeSurfaceMacroscopic2D ******************************** */

template< typename T,template<typename U> class Descriptor>
void FreeSurfaceMacroscopic2D<T,Descriptor>
::processGenericBlocks ( Box2D domain,std::vector<AtomicBlock2D*> atomicBlocks )
{
    typedef Descriptor<T> D;
    using namespace twoPhaseFlag;
    FreeSurfaceProcessorParam2D<T,Descriptor> param ( atomicBlocks );

    T lostMass = param.getSumLostMass();
    plint numInterfaceCells = param.getNumInterfaceCells();
    T massPerCell = T();
    if ( numInterfaceCells>0 )
    {
        massPerCell = lostMass / ( T ) numInterfaceCells;
    }

    // Save macroscopic fields in external scalars and update the mass-fraction.
    for ( plint iX=domain.x0; iX<=domain.x1; ++iX )
    {
        for ( plint iY=domain.y0; iY<=domain.y1; ++iY )
        {
            if ( isWet ( param.flag ( iX,iY ) ) )
            {
                T rhoBar;
                Array<T,2> j;
                momentTemplates<T,Descriptor>::get_rhoBar_j ( param.cell ( iX,iY ), rhoBar, j );
                T density = Descriptor<T>::fullRho ( rhoBar );
                param.setDensity ( iX,iY, density );

                if ( param.flag ( iX,iY ) ==interface )
                {
                    param.mass ( iX,iY ) += massPerCell;
                    T newDensity = param.outsideDensity ( iX,iY );
                    param.volumeFraction ( iX,iY ) = param.mass ( iX,iY ) /newDensity;
                    // On interface cells, adjust the pressure to the ambient pressure.
                    param.setDensity ( iX,iY, newDensity );
                    j *= newDensity/density;
                }
                else if ( isFullWet ( param.flag ( iX,iY ) ) )
                {
                    param.volumeFraction ( iX,iY ) = T ( 1 );
                }

                Array<T,2> force = param.getForce ( iX,iY );
                T tau = T ( 1 ) /param.cell ( iX,iY ).getDynamics().getOmega();
                // Two comments:
                // - Here the force is multiplied by rho0 and not rho so that, under
                //   gravity, a linear pressure profile is obtained.
                // - The force is not multiplied by the volume fraction (some authors
                //   do multiply it by the volumeFraction), because there is a
                //   point-wise interpretation of quantities like momentum.
                j += rhoDefault*tau*force;
                param.setMomentum ( iX,iY,j );
            }
        }
    }
}

/* *************** Class TwoPhaseAddSurfaceTension2D ******************************** */

template< typename T,template<typename U> class Descriptor>
void TwoPhaseAddSurfaceTension2D<T,Descriptor>
::processGenericBlocks ( Box2D domain,std::vector<AtomicBlock2D*> atomicBlocks )
{
    typedef Descriptor<T> D;
    using namespace twoPhaseFlag;
    FreeSurfaceProcessorParam2D<T,Descriptor> param ( atomicBlocks );

    // Save macroscopic fields in external scalars and add the surface tension effect.
    for ( plint iX=domain.x0; iX<=domain.x1; ++iX )
    {
        for ( plint iY=domain.y0; iY<=domain.y1; ++iY )
        {
            if ( param.flag ( iX,iY ) ==interface )
            {
                // This time I do not compute density and momentum from the populations...
                //T rhoBar;
                //Array<T,2> j;
                //momentTemplates<T,Descriptor>::get_rhoBar_j(param.cell(iX,iY,iZ), rhoBar, j);
                //T density = Descriptor<T>::fullRho(rhoBar);
                //param.setDensity(iX,iY, density);

                // ... I just read them from their matrices.
                T density = param.getDensity ( iX,iY );
                Array<T,2> j = param.getMomentum ( iX,iY );

                T newDensity = density;
                // Stored curvature is computed to be twice the mean curvature.
                newDensity += surfaceTension * param.curvature ( iX,iY ) * D::invCs2;
                param.volumeFraction ( iX,iY ) = param.mass ( iX,iY ) / newDensity;
                // On interface cells, adjust the pressure to incorporate surface tension.
                param.setDensity ( iX,iY, newDensity );
                Array<T,2> newJ = j*newDensity/density;
                param.setMomentum ( iX,iY, newJ );

                // TODO Are the following lines really necessary? To be tested.
                Cell<T,Descriptor>& cell = param.cell ( iX,iY );
                T oldRhoBar;
                Array<T,2> oldJ;
                momentTemplates<T,Descriptor>::get_rhoBar_j ( cell, oldRhoBar, oldJ );
                T oldJsqr = normSqr ( oldJ );
                T newRhoBar = Descriptor<T>::rhoBar ( newDensity );
                T newJsqr = normSqr ( newJ );
                for ( int iPop=0; iPop<Descriptor<T>::q; ++iPop )
                {
                    T oldEq = cell.getDynamics().computeEquilibrium ( iPop, oldRhoBar, oldJ, oldJsqr );
                    T newEq = cell.getDynamics().computeEquilibrium ( iPop, newRhoBar, newJ, newJsqr );
                    cell[iPop] += newEq - oldEq;
                }
            }
        }
    }
}

/* *************** Class FreeSurfaceComputeInterfaceLists2D ******************************************* */

template< typename T, template<typename> class Descriptor>
T FreeSurfaceComputeInterfaceLists2D<T,Descriptor>::kappa = 1.e-3;

template< typename T,template<typename U> class Descriptor>
void FreeSurfaceComputeInterfaceLists2D<T,Descriptor>
::processGenericBlocks ( Box2D domain,std::vector<AtomicBlock2D*> atomicBlocks )
{
    typedef Descriptor<T> D;
    typedef typename InterfaceLists<T,Descriptor>::Node Node;
    FreeSurfaceProcessorParam2D<T,Descriptor> param ( atomicBlocks );
    using namespace twoPhaseFlag;

    param.massExcess().clear();
    param.interfaceToFluid().clear();
    param.interfaceToEmpty().clear();
    param.emptyToInterface().clear();

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
                    param.interfaceToFluid().insert ( node );
                }
                else if ( param.volumeFraction ( iX,iY ) < kappa ) // Interface cell is empty.
                {
                    // Elements are added even if they belong to the envelope, because they may be
                    //   needed further down in the same data processor.
                    param.interfaceToEmpty().insert ( node );
                }
            }
        }
    }

    // Where interface cells have become fluid, neighboring cells must be prevented from
    //   being empty, because otherwise there's no interface cell between empty and fluid.
    typename std::set<Node>::iterator iEle = param.interfaceToFluid().begin();
    for ( ; iEle != param.interfaceToFluid().end(); ++iEle )
    {
        // The node here may belong to the 1st envelope.
        Node node = *iEle;
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
            else if ( contained ( nextX,nextY,domain.enlarge ( 1 ) ) && isEmpty ( param.flag ( nextX,nextY ) ) )
            {
                param.emptyToInterface().insert ( nextNode );
            }
        }
    }
}

/* *************** Class FreeSurfaceIniInterfaceToAnyNodes2D ******************************************* */

template< typename T,template<typename U> class Descriptor>
FreeSurfaceIniInterfaceToAnyNodes2D<T,Descriptor>::FreeSurfaceIniInterfaceToAnyNodes2D ( T rhoDefault_ )
    : rhoDefault ( rhoDefault_ )
{ }

template< typename T,template<typename U> class Descriptor>
void FreeSurfaceIniInterfaceToAnyNodes2D<T,Descriptor>
::processGenericBlocks ( Box2D domain,std::vector<AtomicBlock2D*> atomicBlocks )
{
    typedef Descriptor<T> D;
    typedef typename InterfaceLists<T,Descriptor>::Node Node;
    using namespace twoPhaseFlag;

    FreeSurfaceProcessorParam2D<T,Descriptor> param ( atomicBlocks );

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
            T saveMass = param.mass ( iX,iY );
            param.mass ( iX,iY ) = param.getDensity ( iX,iY );
            param.volumeFraction ( iX,iY ) = ( T ) 1;
            param.flag ( iX,iY ) = fluid;

            T massExcess = saveMass - param.getDensity ( iX,iY );
            param.massExcess().insert ( std::pair<Node,T> ( node,massExcess ) );
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

        if ( contained ( iX,iY,domain.enlarge ( 1 ) ) )
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
                param.flag ( iX,iY ) = empty;
                param.attributeDynamics ( iX,iY, new NoDynamics<T,Descriptor>() );

                T massExcess = param.mass ( iX,iY );
                param.massExcess().insert ( std::pair<Node,T> ( node,massExcess ) );

                param.mass ( iX,iY ) = T();
                param.volumeFraction ( iX,iY ) = T();
                param.setDensity ( iX,iY, rhoDefault );
                param.setForce ( iX,iY, Array<T,2> ( T(),T() ) );
                param.setMomentum ( iX,iY, Array<T,2> ( T(),T() ) );
                for ( plint iPop=1; iPop < D::q; ++iPop )
                {
                    plint nextX = iX+D::c[iPop][0];
                    plint nextY = iY+D::c[iPop][1];

                    // The concurrent read/write on param.flag is not an issue here, because the
                    // result in any case is that all adjacent fluid cells have become interface.
                    if ( param.flag ( nextX,nextY ) ==fluid )
                    {
                        param.flag ( nextX,nextY ) = interface;
                    }
                }
            }
        }
    }
}

/* *************** Class FreeSurfaceIniEmptyToInterfaceNodes2D ******************************************* */

template< typename T,template<typename U> class Descriptor>
void FreeSurfaceIniEmptyToInterfaceNodes2D<T,Descriptor>
::processGenericBlocks ( Box2D domain,std::vector<AtomicBlock2D*> atomicBlocks )
{
    typedef Descriptor<T> D;
    typedef typename InterfaceLists<T,Descriptor>::Node Node;
    using namespace twoPhaseFlag;
    FreeSurfaceProcessorParam2D<T,Descriptor> param ( atomicBlocks );

    // In this data processor, density and momentum are potentially read and written
    //   from the same node, because nodes can switch state. The following two vectors
    //   store temporary variables to avoid read/write in undefined order.
    std::vector<T> newDensity ( param.emptyToInterface().size() );
    std::vector<Array<T,2> > newMomentum ( param.emptyToInterface().size() );
    std::fill ( newDensity.begin(), newDensity.end(), T() );
    std::fill ( newMomentum.begin(), newMomentum.end(), Array<T,2> ( T(),T() ) );

    // Compute density and momentum for cells that will switch state empty->interface.
    //   It is sufficient to do this is bulk+0.
    //   This loop performs read-only access to the lattice.
    plint i=0;
    typename std::set<Node>::iterator iEle = param.emptyToInterface().begin();
    for ( ; iEle != param.emptyToInterface().end(); ++iEle, ++i )
    {
        Node node = *iEle;
        plint iX=node[0];
        plint iY=node[1];

        // If non-bulk elements are left in the list, disregard to avoid accessing undefined neighbors.
        if ( contained ( iX,iY,domain ) )
        {
            // For initialization of the new cell, compute average density
            //   and momentum on neighbors.
            T averageDensity = T ( 0 );
            Array<T,2> averageMomentum ( T ( 0 ),T ( 0 ) );
            T sumWeights = ( T ) 0;
            for ( plint iPop=1; iPop<D::q; ++iPop )
            {
                plint nextX = iX+D::c[iPop][0];
                plint nextY = iY+D::c[iPop][1];

                // Warning: it is not accounted for the fact that neighbors can have excess mass. It
                //   might be good to account for this in the future.
                if ( isWet ( param.flag ( nextX,nextY ) ) )
                {
                    T weight = D::t[iPop];
                    sumWeights += weight;
                    averageDensity += weight * param.getDensity ( nextX,nextY );
                    averageMomentum += weight * param.getMomentum ( nextX,nextY );
                }
            }
            T invSum = T ( 1 ) /sumWeights;
            averageDensity  *= invSum;
            averageMomentum *= invSum;
            newDensity[i] = averageDensity;
            newMomentum[i] = averageMomentum;
        }
    }

    // Elements that have switched state empty->interface are initialized at equilibrium.
    //   It is sufficient to initialize them in bulk+0.
    //   This loop performs write-only access on the lattice.
    i=0;
    iEle = param.emptyToInterface().begin();
    for ( ; iEle != param.emptyToInterface().end(); ++iEle, ++i )
    {
        Node node = *iEle;

        plint iX=node[0];
        plint iY=node[1];

        // If non-bulk elements are left in the list, disregard to avoid accessing undefined neighbors.
        if ( contained ( iX,iY,domain ) )
        {
            T averageDensity = newDensity[i];
            Array<T,2> averageMomentum = newMomentum[i];

            param.attributeDynamics ( iX,iY, dynamicsTemplate->clone() );

            iniCellAtEquilibrium ( param.cell ( iX,iY ), averageDensity, averageMomentum/averageDensity );
            param.setForce ( iX,iY, force );
            // Change density, but leave mass and volumeFraction at 0, as they are later
            //   recomputed (Warning: this is probably correct, but there remains a small doubt).
            param.setMomentum ( iX,iY, averageMomentum );
            param.setDensity ( iX,iY, averageDensity );
            param.mass ( iX,iY ) = T();
            param.volumeFraction ( iX,iY ) = T();
            param.flag ( iX,iY ) = interface;
        }
    }
}

/* *************** Class FreeSurfaceRemoveFalseInterfaceCells2D ******************************************* */

template< typename T,template<typename U> class Descriptor>
void FreeSurfaceRemoveFalseInterfaceCells2D<T,Descriptor>
::processGenericBlocks ( Box2D domain,std::vector<AtomicBlock2D*> atomicBlocks )
{
    typedef Descriptor<T> D;
    typedef typename InterfaceLists<T,Descriptor>::Node Node;
    using namespace twoPhaseFlag;

    FreeSurfaceProcessorParam2D<T,Descriptor> param ( atomicBlocks );

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
                    }
                }
            }
        }
    }

    for ( pluint i=0; i<interfaceToFluidNodes.size(); ++i )
    {
        Node const& pos = interfaceToFluidNodes[i];
        param.flag ( pos[0],pos[1] ) = fluid;
    }
    for ( pluint i=0; i<interfaceToEmptyNodes.size(); ++i )
    {
        Node const& pos = interfaceToEmptyNodes[i];
        param.flag ( pos[0],pos[1] ) = empty;
    }
}

/* *************** Class FreeSurfaceEqualMassExcessReDistribution2D ******************************************* */

template< typename T,template<typename U> class Descriptor>
void FreeSurfaceEqualMassExcessReDistribution2D<T,Descriptor>
::processGenericBlocks ( Box2D domain,std::vector<AtomicBlock2D*> atomicBlocks )
{
    typedef Descriptor<T> D;
    typedef typename InterfaceLists<T,Descriptor>::Node Node;
    using namespace twoPhaseFlag;
    FreeSurfaceProcessorParam2D<T,Descriptor> param ( atomicBlocks );

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
                        param.mass ( nextX,nextY ) / param.getDensity ( nextX,nextY );
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
}

template< typename T,template<typename U> class Descriptor>
void FreeSurfaceInterfaceFilter2D<T,Descriptor>
::processGenericBlocks ( Box2D domain,std::vector<AtomicBlock2D*> atomicBlocks )
{
    typedef Descriptor<T> D;
    using namespace twoPhaseFlag;
    FreeSurfaceProcessorParam2D<T,Descriptor> param ( atomicBlocks );

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

    std::vector<T> newRho;
    std::vector<Array<T,2> > newJ;
    for ( pluint i=0; i<interfaceNodes.size(); ++i )
    {
        Array<plint,2> node = interfaceNodes[i];
        plint iX=node[0];
        plint iY=node[1];
        T rho = T();
        Array<T,2> j;
        j.resetToZero();
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
                    j += newWeight*param.getMomentum ( px,py );
                    weight += newWeight;
                    ++numNodes;
                }
            }
        }
        PLB_ASSERT ( numNodes>0 );
        rho /= weight;
        j /= weight;
        newRho.push_back ( rho );
        newJ.push_back ( j );
    }
    for ( pluint i=0; i<interfaceNodes.size(); ++i )
    {
        Array<plint,2> node = interfaceNodes[i];
        plint iX=node[0];
        plint iY=node[1];

        param.setDensity ( iX,iY, newRho[i] );
        param.setMomentum ( iX,iY, newJ[i] );
    }
}

/* *************** Class TwoPhaseComputeStatistics2D ******************************************* */

template< typename T,template<typename U> class Descriptor>
void TwoPhaseComputeStatistics2D<T,Descriptor>
::processGenericBlocks ( Box2D domain,std::vector<AtomicBlock2D*> atomicBlocks )
{
    typedef Descriptor<T> D;
    using namespace twoPhaseFlag;

    FreeSurfaceProcessorParam2D<T,Descriptor> param ( atomicBlocks );

    for ( plint iX=domain.x0; iX<=domain.x1; ++iX )
    {
        for ( plint iY=domain.y0; iY<=domain.y1; ++iY )
        {
            if ( isWet ( param.flag ( iX,iY ) ) )
            {
                param.addToTotalMass ( param.mass ( iX,iY ) );
                if ( param.flag ( iX,iY ) ==interface )
                {
                    param.addToInterfaceCells ( 1 );
                }
            }
        }
    }
}

}  // namespace plb

#endif  // FREE_SURFACE_MODEL_2D_HH


