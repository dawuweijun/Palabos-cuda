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

/* Main author: Orestis Malaspinas */

#ifndef VOXELIZER2D_HH
#define VOXELIZER2D_HH

#include "core/globalDefs.h"
#include "offLattice/voxelizer2D.h"
#include "atomicBlock/dataField2D.h"
#include "multiBlock/multiBlockGenerator2D.h"

namespace plb
{

namespace voxelFlag
{
inline int invert ( int arg )
{
    switch ( arg )
    {
    case inside:
        return outside;
    case outside:
        return inside;
    case innerBorder:
        return outerBorder;
    case outerBorder:
        return innerBorder;
    case undetermined:
        return undetermined;
    default:
        PLB_ASSERT ( false );
    }
    return undetermined;
}
inline int bulkFlag ( int arg )
{
    if ( arg==innerBorder || arg==inside )
    {
        return inside;
    }
    else if ( arg==outerBorder || arg==outside )
    {
        return outside;
    }
    else
    {
        return undetermined;
    }
}
inline int borderFlag ( int arg )
{
    if ( arg==inside || arg==innerBorder )
    {
        return innerBorder;
    }
    else if ( arg==outside || arg==outerBorder )
    {
        return outerBorder;
    }
    else
    {
        return undetermined;
    }
}
inline bool insideFlag ( int arg )
{
    return arg==inside || arg==innerBorder;
}
inline bool outsideFlag ( int arg )
{
    return arg==outside || arg==outerBorder;
}

}  // namespace voxelFlag

template<typename T>
std::auto_ptr<MultiScalarField2D<int> > voxelize (
    SegmentPolygonMesh2D<T> const& mesh,
    plint symmetricLayer, plint borderWidth )
{
    Array<T,2> xRange, yRange;
    mesh.computeBoundingBox ( xRange, yRange );
    // Creation of the multi-scalar field. The +1 is because if the resolution is N,
    //   the number of nodes is N+1.
    plint nx = ( plint ) ( xRange[1] - xRange[0] ) + 1 + 2*symmetricLayer;
    plint ny = ( plint ) ( yRange[1] - yRange[0] ) + 1 + 2*symmetricLayer;

    return voxelize ( mesh, Box2D ( 0,nx-1, 0,ny-1 ), borderWidth );
}

template<typename T>
std::auto_ptr<MultiScalarField2D<int> > voxelize (
    SegmentPolygonMesh2D<T> const& mesh,
    Box2D const& domain, plint borderWidth )
{
    // As initial seed, a one-cell layer around the outer boundary is tagged
    //   as ouside cells.
    plint envelopeWidth=1;
    std::auto_ptr<MultiScalarField2D<int> > voxelMatrix
        = generateMultiScalarField<int> ( domain, voxelFlag::outside, envelopeWidth );
    setToConstant ( *voxelMatrix, voxelMatrix->getBoundingBox().enlarge ( -1 ),
                    voxelFlag::undetermined );

    MultiContainerBlock2D hashContainer ( *voxelMatrix );
    std::vector<MultiBlock2D*> container_arg;
    container_arg.push_back ( &hashContainer );
    applyProcessingFunctional (
        new CreateSegmentHash<T> ( mesh ),
        hashContainer.getBoundingBox(), container_arg );

    std::vector<MultiBlock2D*> flag_hash_arg;
    flag_hash_arg.push_back ( voxelMatrix.get() );
    flag_hash_arg.push_back ( &hashContainer );

    voxelMatrix->resetFlags(); // Flags are used internally by VoxelizeMeshFunctional2D.
    while ( !allFlagsTrue ( voxelMatrix.get() ) )
    {
        applyProcessingFunctional (
            new VoxelizeMeshFunctional2D<T> ( mesh ),
            voxelMatrix->getBoundingBox(), flag_hash_arg );
    }

    detectBorderLine ( *voxelMatrix, voxelMatrix->getBoundingBox(), borderWidth );

    return std::auto_ptr<MultiScalarField2D<int> > ( voxelMatrix );
}

template<typename T>
std::auto_ptr<MultiScalarField2D<int> > voxelize (
    SegmentPolygonMesh2D<T> const& mesh,
    Box2D const& domain, plint borderWidth, Box2D seed )
{
    // As initial seed, a one-cell layer around the outer boundary is tagged
    //   as ouside cells.
    plint envelopeWidth=1;

    std::auto_ptr<MultiScalarField2D<int> > voxelMatrix
        = generateMultiScalarField<int> ( domain, voxelFlag::undetermined, envelopeWidth );
    setToConstant ( *voxelMatrix, seed, voxelFlag::outside );

    MultiContainerBlock2D hashContainer ( *voxelMatrix );
    std::vector<MultiBlock2D*> container_arg;
    container_arg.push_back ( &hashContainer );
    applyProcessingFunctional (
        new CreateSegmentHash<T> ( mesh ),
        hashContainer.getBoundingBox(), container_arg );

    std::vector<MultiBlock2D*> flag_hash_arg;
    flag_hash_arg.push_back ( voxelMatrix.get() );
    flag_hash_arg.push_back ( &hashContainer );

    voxelMatrix->resetFlags(); // Flags are used internally by VoxelizeMeshFunctional2D.
    plint maxIteration=100;
    plint i=0;
    while ( !allFlagsTrue ( voxelMatrix.get() ) && i<maxIteration )
    {
        applyProcessingFunctional (
            new VoxelizeMeshFunctional2D<T> ( mesh ),
            voxelMatrix->getBoundingBox(), flag_hash_arg );
        ++i;
    }
    if ( i==maxIteration )
    {
        pcout << "Warning: Voxelization failed." << std::endl;
    }

    detectBorderLine ( *voxelMatrix, voxelMatrix->getBoundingBox(), borderWidth );

    return std::auto_ptr<MultiScalarField2D<int> > ( voxelMatrix );
}


template<typename T>
std::auto_ptr<MultiScalarField2D<int> > revoxelize (
    SegmentPolygonMesh2D<T> const& mesh,
    MultiScalarField2D<int>& oldVoxelMatrix,
    MultiContainerBlock2D& hashContainer, plint borderWidth )
{
    // As initial seed, a one-cell layer around the outer boundary is tagged
    //   as ouside cells.
    Box2D domain ( oldVoxelMatrix.getBoundingBox() );
    std::auto_ptr<MultiScalarField2D<int> > voxelMatrix (
        new MultiScalarField2D<int> ( ( MultiBlock2D& ) oldVoxelMatrix ) );
    setToConstant ( *voxelMatrix, domain, voxelFlag::outside );
    setToConstant ( *voxelMatrix, voxelMatrix->getBoundingBox().enlarge ( -1 ),
                    voxelFlag::undetermined );

    std::vector<MultiBlock2D*> flag_hash_arg;
    flag_hash_arg.push_back ( voxelMatrix.get() );
    flag_hash_arg.push_back ( &hashContainer );

    voxelMatrix->resetFlags(); // Flags are used internally by VoxelizeMeshFunctional2D.
    while ( !allFlagsTrue ( voxelMatrix.get() ) )
    {
        applyProcessingFunctional (
            new VoxelizeMeshFunctional2D<T> ( mesh ),
            voxelMatrix->getBoundingBox(), flag_hash_arg );
    }

    detectBorderLine ( *voxelMatrix, voxelMatrix->getBoundingBox(), borderWidth );

    return std::auto_ptr<MultiScalarField2D<int> > ( voxelMatrix );
}


/* ******** VoxelizeMeshFunctional2D ************************************* */

template<typename T>
VoxelizeMeshFunctional2D<T>::VoxelizeMeshFunctional2D (
    SegmentPolygonMesh2D<T> const& mesh_ )
    : mesh ( mesh_ )
{ }

template<typename T>
bool VoxelizeMeshFunctional2D<T>::distanceToSurface (
    AtomicContainerBlock2D& hashContainer,
    Array<T,2> const& point, T& distance, bool& isBehind ) const
{
    T maxDistance = sqrt ( 3 );
    Array<T,2> xRange ( point[0]-maxDistance, point[0]+maxDistance );
    Array<T,2> yRange ( point[1]-maxDistance, point[1]+maxDistance );
    SegmentHash<T> segmentHash ( hashContainer );
    std::vector<plint> possibleSegments;
    segmentHash.getSegments ( xRange, yRange, possibleSegments );

    T    tmpDistance;
    bool tmpIsBehind;
    bool segmentFound = false;

    for ( pluint iPossible=0; iPossible<possibleSegments.size(); ++iPossible )
    {
        plint iSegment = possibleSegments[iPossible];
        mesh.distanceToSegment (
            point, iSegment, tmpDistance, tmpIsBehind );
        if ( !segmentFound || tmpDistance<distance )
        {
            distance = tmpDistance;
            isBehind = tmpIsBehind;
            segmentFound = true;
        }
    }
    return segmentFound;
}

template<typename T>
bool VoxelizeMeshFunctional2D<T>::checkIfFacetsCrossed (
    AtomicContainerBlock2D& hashContainer,
    Array<T,2> const& point1, Array<T,2> const& point2,
    T& distance, plint& whichSegment )
{
    Array<T,2> xRange (
        std::min ( point1[0], point2[0] ),
        std::max ( point1[0], point2[0] ) );
    Array<T,2> yRange (
        std::min ( point1[1], point2[1] ),
        std::max ( point1[1], point2[1] ) );
    SegmentHash<T> segmentHash ( hashContainer );
    std::vector<plint> possibleSegments;
    segmentHash.getSegments ( xRange, yRange, possibleSegments );

    int flag = 0; // Check for crossings inside the point1-point2 segment.
    Array<T,2> intersection; // Dummy variable.
    Array<T,2> normal;       // Dummy variable.
    T tmpDistance;           // Dummy variable.

    if ( global::counter ( "voxelizer-debug" ).getCount() ==1 )
    {
        std::cout << "{";
    }
    std::vector<T> crossings;
    for ( pluint iPossible=0; iPossible<possibleSegments.size(); ++iPossible )
    {
        plint iSegment = possibleSegments[iPossible];
        if ( mesh.pointOnSegment ( point1, point2, flag, iSegment, intersection, normal, tmpDistance ) ==1 )
        {
            if ( global::counter ( "voxelizer-debug" ).getCount() ==1 )
            {
                std::cout << "(" << iSegment << ";" << tmpDistance << ")";
            }
            crossings.push_back ( tmpDistance );
            if ( crossings.size() ==1 || tmpDistance<distance )
            {
                distance = tmpDistance;
                whichSegment = iSegment;
            }
        }
    }
    if ( global::counter ( "voxelizer-debug" ).getCount() ==1 )
    {
        std::cout << "}";
    }

    if ( crossings.size() ==0 )
    {
        return false;
    }
    else
    {
        bool hasCrossed = true;
        for ( pluint iCrossing=1; iCrossing<crossings.size(); ++iCrossing )
        {
            //const T eps1 = std::numeric_limits<double>::epsilon()*1.e2;
            //if ( !util::fpequal(crossings[iCrossing], crossings[iCrossing-1], eps1) )
            const T eps1 = std::numeric_limits<double>::epsilon() *1.e4;
            if ( fabs ( crossings[iCrossing]-crossings[iCrossing-1] ) >eps1 )
            {
                hasCrossed = !hasCrossed;
            }
        }
        return hasCrossed;
    }
}

template<typename T>
bool VoxelizeMeshFunctional2D<T>::createVoxelizationRange (
    Box2D const& domain, ScalarField2D<int>& voxels,
    Array<plint,2>& xRange, Array<plint,2>& yRange )
{
    // The purpose of the first three loops is to locate the eight
    //   corners of the cube. One voxel per corner would be insufficient
    //   because a potential seed is situated differently, depending on
    //   whether it is on the boundary of the multi-block or somewhere inside.
    for ( plint dx=0; dx<=+1; ++dx )
    {
        plint xMin = domain.x0+dx*domain.getNx()-1;
        plint xMax = domain.x0+dx*domain.getNx();
        for ( plint dy=0; dy<=+1; ++dy )
        {
            plint yMin = domain.y0+dy*domain.getNy()-1;
            plint yMax = domain.y0+dy*domain.getNy();
            // Locate a potential seed in one of the corners.
            for ( plint iX=xMin; iX<=xMax; ++iX )
            {
                for ( plint iY=yMin; iY<=yMax; ++iY )
                {
                    if ( voxels.get ( iX,iY ) != voxelFlag::undetermined )
                    {
                        xRange[0] = domain.x0+dx* ( domain.getNx()-1 );
                        xRange[1] = domain.x0+ ( 1-dx ) * ( domain.getNx()-1 );
                        yRange[0] = domain.y0+dy* ( domain.getNy()-1 );
                        yRange[1] = domain.y0+ ( 1-dy ) * ( domain.getNy()-1 );
                        return true;
                    }
                }
            }

        }
    }
    return false;
}

template<typename T>
void VoxelizeMeshFunctional2D<T>::printOffender (
    ScalarField2D<int> const& voxels,
    AtomicContainerBlock2D& hashContainer,
    Dot2D pos )
{
    std::set<plint> segments;
    Dot2D offset = voxels.getLocation();
    Dot2D pos_ = pos+offset;
    std::cout << "Position (" << pos_.x << "," << pos_.y << ")" << std::endl;
    for ( plint dx=-1; dx<=+1; ++dx )
    {
        for ( plint dy=-1; dy<=+1; ++dy )
        {
            if ( ! ( dx==0 && dy==0 ) )
            {
                Dot2D neigh = pos+offset+Dot2D ( dx,dy );
                int typeOfNeighbor = voxels.get ( pos.x+dx,pos.y+dy );
                if ( typeOfNeighbor!=voxelFlag::undetermined )
                {
                    T distance;
                    plint whichSegment;
                    Array<T,2> p1 ( pos_.x,pos_.y );
                    Array<T,2> p2 ( neigh.x,neigh.y );
                    global::counter ( "voxelizer-debug" ).increment ( 1 );
                    bool crossed = checkIfFacetsCrossed (
                                       hashContainer, p1, p2, distance, whichSegment );
                    global::counter ( "voxelizer-debug" ).reset();
                    std::cout << "Neighbor ("
                              << dx << "," << dy
                              << "); is "
                              << ( voxelFlag::insideFlag ( typeOfNeighbor ) ? "inside" : "outside" );
                    if ( crossed )
                    {
                        segments.insert ( whichSegment );
                        std::cout
                                << " inters. at distance " << distance
                                << " with segment " << whichSegment << std::endl;
                    }
                    else
                    {
                        std::cout << " no inters." << std::endl;
                    }
                }
            }
        }
    }
    std::set<plint>::iterator it = segments.begin();
    for ( ; it!=segments.end(); ++it )
    {
        std::cout << "Segment " << *it << " [" << std::flush;
        Array<T,2> p0 = mesh.getVertex ( *it, 0 );
        Array<T,2> p1 = mesh.getVertex ( *it, 1 );
        Array<T,2> p2 = mesh.getVertex ( *it, 2 );
        std::cout << p0[0] << " " << p1[0] << " " << p2[0] << " " << p0[0] << "], ["
                  << p0[1] << " " << p1[1] << " " << p2[1] << " " << p0[1] << "], ["
                  << p0[2] << " " << p1[2] << " " << p2[2] << " " << p0[2] << "]" << std::endl;
    }
}

template<typename T>
bool VoxelizeMeshFunctional2D<T>::voxelizeFromNeighbor (
    ScalarField2D<int> const& voxels,
    AtomicContainerBlock2D& hashContainer,
    Dot2D pos, Dot2D neighbor, int& voxelType )
{
    int verificationLevel = 0;
    Dot2D offset = voxels.getLocation();
    int typeOfNeighbor = voxels.get ( neighbor.x,neighbor.y );
    if ( typeOfNeighbor==voxelFlag::undetermined )
    {
        return true;
    }
    // If there is no verification and the voxel has already been voxelized,
    //   it is not being re-voxelized here.
    if ( verificationLevel==0 )
    {
        if ( voxelType!=voxelFlag::undetermined )
        {
            return true;
        }
    }
    Dot2D pos_ = pos+offset;
    Dot2D neighbor_ = neighbor+offset;
    Array<T,2> point1 ( ( T ) pos_.x, ( T ) pos_.y );
    Array<T,2> point2 ( ( T ) neighbor_.x, ( T ) neighbor_.y );
    int newVoxelType = voxelFlag::undetermined;
    T distance1, distance2, distance3, distance4;
    bool isBehind1, isBehind2;
    plint whichSegment1, whichSegment2;
    if ( checkIfFacetsCrossed ( hashContainer, point1, point2, distance1, whichSegment1 ) )
    {
        newVoxelType = voxelFlag::invert ( typeOfNeighbor );
        // Additional consistency checks only at the ultimate level of verification.
        if ( verificationLevel==2 )
        {
            PLB_ASSERT ( distance1 < sqrt ( ( T ) 3 ) + ( T ) 0.0001 );
#ifdef PLB_DEBUG
            bool ok = checkIfFacetsCrossed ( hashContainer, point2, point1, distance2, whichSegment2 );
#else
            ( void ) checkIfFacetsCrossed ( hashContainer, point2, point1, distance2, whichSegment2 );
#endif
            PLB_ASSERT ( ok );
            PLB_ASSERT ( distance2 < sqrt ( ( T ) 3 ) + ( T ) 0.0001 );

#ifdef PLB_DEBUG
            bool ok1 = distanceToSurface ( hashContainer, point1, distance3, isBehind1 );
#else
            ( void ) distanceToSurface ( hashContainer, point1, distance3, isBehind1 );
#endif

            PLB_ASSERT ( ok1 );
            PLB_ASSERT ( distance1 < sqrt ( ( T ) 3 ) + ( T ) 0.0001 );
            // Attention: At this moment, the following consistency check fails sometimes,
            //   god knows why. It might be that there is a bug in the method
            //   mesh.distanceToSurface.
            PLB_ASSERT ( ( voxelFlag::insideFlag ( newVoxelType ) && isBehind1 ) ||
                         ( voxelFlag::outsideFlag ( newVoxelType ) && !isBehind1 ) );

#ifdef PLB_DEBUG
            bool ok2 = distanceToSurface ( hashContainer, point2, distance4, isBehind2 );
#else
            ( void ) distanceToSurface ( hashContainer, point2, distance4, isBehind2 );
#endif
            PLB_ASSERT ( ok2 );
            PLB_ASSERT ( distance2 < sqrt ( ( T ) 3 ) + ( T ) 0.0001 );
            PLB_ASSERT ( ( voxelFlag::insideFlag ( typeOfNeighbor ) && isBehind2 ) ||
                         ( voxelFlag::outsideFlag ( typeOfNeighbor ) && !isBehind2 ) );
        }
    }
    else
    {
        newVoxelType = typeOfNeighbor;
    }
    int oldVoxelType = voxelType;
    voxelType = newVoxelType;
    if ( oldVoxelType == voxelFlag::undetermined )
    {
        return true;
    }
    else
    {
        return oldVoxelType == newVoxelType;
    }
}

template<typename T>
void VoxelizeMeshFunctional2D<T>::processGenericBlocks (
    Box2D domain, std::vector<AtomicBlock2D*> blocks )
{
    PLB_PRECONDITION ( blocks.size() ==2 );
    ScalarField2D<int>* voxels =
        dynamic_cast<ScalarField2D<int>*> ( blocks[0] );
    PLB_ASSERT ( voxels );
    AtomicContainerBlock2D* container =
        dynamic_cast<AtomicContainerBlock2D*> ( blocks[1] );
    PLB_ASSERT ( container );

    // Return if this block is already voxelized.
    if ( voxels->getFlag() )
    {
        return;
    }

    Array<plint,2> xRange, yRange;
    if ( !createVoxelizationRange ( domain, *voxels, xRange, yRange ) )
    {
        // If no seed has been found in the envelope, just return and wait
        //   for the next round.
        return;
    }

    // Specify if the loops go in positive or negative direction.
    plint xIncr = xRange[1]>xRange[0] ? 1 : -1;
    plint yIncr = yRange[1]>yRange[0] ? 1 : -1;
    // The ranges are closed on both ends. Here, the range[1] end
    //   is converted to an open one so we can use != checks in the loops.
    xRange[1] += xIncr;
    yRange[1] += yIncr;
    for ( plint iX=xRange[0]; iX!=xRange[1]; iX+=xIncr )
    {
        for ( plint iY=yRange[0]; iY!=yRange[1]; iY+=yIncr )
        {
            Dot2D pos ( iX,iY );
            int voxelType = voxels->get ( iX,iY );
            if ( voxelType==voxelFlag::undetermined )
            {
                for ( plint dx=-1; dx<=+1; ++dx )
                {
                    for ( plint dy=-1; dy<=+1; ++dy )
                    {
                        if ( ! ( dx==0 && dy==0 ) )
                        {
                            Dot2D neighbor ( iX+dx, iY+dy );
                            bool ok = voxelizeFromNeighbor (
                                          *voxels, *container,
                                          pos, neighbor, voxelType );
                            if ( !ok )
                            {
                                printOffender ( *voxels, *container, pos );
                            }
                            PLB_ASSERT ( ok );
                        }
                    }
                }
                voxels->get ( iX,iY ) = voxelType;
            }
        }
    }
    // Indicate that this atomic-block has been voxelized.
    voxels->setFlag ( true );
}

template<typename T>
VoxelizeMeshFunctional2D<T>* VoxelizeMeshFunctional2D<T>::clone() const
{
    return new VoxelizeMeshFunctional2D<T> ( *this );
}

template<typename T>
void VoxelizeMeshFunctional2D<T>::getTypeOfModification ( std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::staticVariables;  // Voxels
    modified[1] = modif::nothing; // Hash Container
}

template<typename T>
BlockDomain::DomainT VoxelizeMeshFunctional2D<T>::appliesTo() const
{
    return BlockDomain::bulk;
}



/* ******** DetectBorderLineFunctional2D ************************************* */

template<typename T>
void detectBorderLine ( MultiScalarField2D<T>& voxelMatrix,
                        Box2D const& domain, plint borderWidth )
{
    applyProcessingFunctional ( new DetectBorderLineFunctional2D<T> ( borderWidth ),
                                domain, voxelMatrix );
}

template<typename T>
DetectBorderLineFunctional2D<T>::DetectBorderLineFunctional2D ( plint borderWidth_ )
    : borderWidth ( borderWidth_ )
{ }

template<typename T>
void DetectBorderLineFunctional2D<T>::process (
    Box2D domain, ScalarField2D<T>& voxels )
{
    for ( plint iX = domain.x0; iX <= domain.x1; ++iX )
    {
        for ( plint iY = domain.y0; iY <= domain.y1; ++iY )
        {
            for ( plint dx=-borderWidth; dx<=borderWidth; ++dx )
                for ( plint dy=-borderWidth; dy<=borderWidth; ++dy )
                    if ( ! ( dx==0 && dy==0 ) )
                    {
                        plint nextX = iX + dx;
                        plint nextY = iY + dy;
                        if ( contained ( Dot2D ( nextX,nextY ),voxels.getBoundingBox() ) )
                        {
                            if ( voxelFlag::outsideFlag ( voxels.get ( iX,iY ) ) &&
                                    voxelFlag::insideFlag ( voxels.get ( nextX,nextY ) ) )
                            {
                                voxels.get ( iX,iY ) = voxelFlag::outerBorder;
                            }
                            if ( voxelFlag::insideFlag ( voxels.get ( iX,iY ) ) &&
                                    voxelFlag::outsideFlag ( voxels.get ( nextX,nextY ) ) )
                            {
                                voxels.get ( iX,iY ) = voxelFlag::innerBorder;
                            }
                        }
                    }
        }
    }
}

template<typename T>
DetectBorderLineFunctional2D<T>* DetectBorderLineFunctional2D<T>::clone() const
{
    return new DetectBorderLineFunctional2D<T> ( *this );
}

template<typename T>
void DetectBorderLineFunctional2D<T>::getTypeOfModification ( std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::staticVariables;
}

template<typename T>
BlockDomain::DomainT DetectBorderLineFunctional2D<T>::appliesTo() const
{
    return BlockDomain::bulk;
}

} // namespace plb

#endif  // VOXELIZER2D_HH
