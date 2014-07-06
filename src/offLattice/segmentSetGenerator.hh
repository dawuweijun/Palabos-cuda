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

/* Main author: Dimitrios Kontaxakis */

#ifndef SEGMENT_SET_GENERATOR_HH
#define SEGMENT_SET_GENERATOR_HH

#include "core/globalDefs.h"
#include "offLattice/segmentSetGenerator.h"
#include "offLattice/segmentPolygonMesh2D.h"
#include "offLattice/segmentBoundary2D.h"
#include "offLattice/segmentPolygonMesh2D.h"

namespace plb
{

template<typename T>
SegmentSet<T>* constructCircle ( Array<T,2> const& center, T radius, plint minNumOfSegments )
{
    std::vector<typename SegmentSet<T>::Segment> segments;
#ifdef PLB_DEBUG
    static const T eps = std::numeric_limits<T>::epsilon();
#endif
    PLB_ASSERT ( radius > ( T ) 0.0 && !util::fpequal ( radius, ( T ) 0.0, eps ) &&
                 minNumOfSegments > 3 );

    // Create a triangularized unit circle

    // Initial 4 vertices

    Array<T,2> va(1.0,0.0);

    Array<T,2> vb( 0.0,1.0);

    Array<T,2> vc(-1.0,0.0);

    Array<T,2> vd(0.0 ,-1 );

    // Initial 4 segments

    typename SegmentSet<T>::Segment tmp;

    tmp[0] = vd;
    tmp[1] = va;
    segments.push_back ( tmp );

    tmp[0] = va;
    tmp[1] = vb;
    segments.push_back ( tmp );

    tmp[0] = vb;
    tmp[1] = vc;
    segments.push_back ( tmp );

    tmp[0] = vc;
    tmp[1] = vd;
    segments.push_back ( tmp );

    // Perform refinement iterations

    plint size;
    while ( ( size = segments.size() ) < minNumOfSegments )
    {
        for ( plint i = 0; i < size; i++ )
        {
            va = segments[i][0];
            vb = segments[i][1];

            vc = ( T ) 0.5 * ( va + vb );

            vc /= norm ( vc );

            segments[i][0] = va;
            segments[i][1] = vc;

            tmp[0] = vc;
            tmp[1] = vb;

            segments.push_back ( tmp );
        }
    }

    // Scale and translate the mesh

    SegmentSet<T> *segmentSet=new SegmentSet<T> ( segments );

    segmentSet->scale ( radius );
    segmentSet->translate ( center );

    return segmentSet;
}

template<typename T>
void addSurface (
    Array<T,2> const& lowerCorner,
    Array<T,2> const& delta1, plint n1, Array<T,2> const& delta2, plint n2,
    std::vector<typename SegmentSet<T>::Segment>& segments )
{
    Array<T,2> pos1 ( lowerCorner );
    for ( plint i1=0; i1<n1; ++i1, pos1+=delta1 )
    {
        Array<T,2> pos2 ( pos1 );
        for ( plint i2=0; i2<n2; ++i2, pos2+=delta2 )
        {
            typename SegmentSet<T>::Segment segment;
            segment[0] = pos2;
            segment[1] = pos2+delta1;
            segments.push_back ( segment );
            segment[0] += delta1+delta2;
            std::swap ( segment[1], segment[2] );
            segments.push_back ( segment );
        }
    }
}

template<typename T>
SegmentSet<T> patchTubes ( SegmentSet<T> const& geometryWithOpenings, plint sortDirection, std::vector<T> patchLengths )
{
    typedef typename SegmentSet<T>::Segment Segment;

    std::vector<Segment> fullGeometry ( geometryWithOpenings.getSegments() );

    DEFscaledMesh2D<T>* defMesh = new DEFscaledMesh2D<T> ( geometryWithOpenings );
    SegmentPolygonMesh2D<T>& mesh = defMesh->getMesh();

    std::vector<Lid> holes = mesh.closeHoles();
    std::sort ( holes.begin(), holes.end(), LidLessThan2D<T> ( sortDirection, mesh ) );

    PLB_ASSERT ( holes.size() == patchLengths.size() );

    for ( pluint iHole=0; iHole<holes.size(); ++iHole )
    {
        Array<T,2> baryCenter = computeGeometricCenter2D ( mesh,holes[iHole] );
        plint numHoleVertices = ( plint ) holes[iHole].boundaryVertices.size();

        Array<T,2> normal = computeNormal ( mesh, holes[iHole] );
        T aveRadius = computeGeometricRadius2D ( mesh,holes[iHole] );

        Array<T,2> nextCenter = baryCenter + normal*aveRadius* ( T ) 10./ ( T ) numHoleVertices;

        plint numInletPoints = numHoleVertices;
        bool oddNumber = numHoleVertices%2==1;
        if ( oddNumber ) numInletPoints--; // Must be even for cylinder construction algorithm.
        std::vector<Array<T,2> > inletPoints;
        plint numPointsOnLength = numInletPoints*patchLengths[iHole]/aveRadius/8;
        if ( numPointsOnLength<3 ) numPointsOnLength = 3;
        SegmentSet<T> piece = constructCircle ( nextCenter, normal, aveRadius, aveRadius, patchLengths[iHole],
                                                numPointsOnLength, numInletPoints/2, inletPoints );
        std::vector<Segment> pieceSegments = piece.getSegments();

        plint newId = 0;
        T minDistance = std::numeric_limits<T>::max();
        plint minDistanceId = -1;
        for ( plint i=0; i<numHoleVertices; ++i )
        {
            plint iVertex = holes[iHole].boundaryVertices[i];
            Array<T,2> p1 = mesh.getVertex ( iVertex );
            T nextDistance = norm ( inletPoints[newId]-p1 );
            if ( nextDistance<minDistance )
            {
                minDistance = nextDistance;
                minDistanceId = i;
            }
        }
        plint oldId = minDistanceId;
        plint newId_p1 = 0;
        for ( plint i=0; i<numInletPoints; ++i )
        {
            plint newId_p1 = ( newId+1 ) % numInletPoints;
            plint oldId_p1 = oldId-1;
            if ( oldId_p1<0 ) oldId_p1 = numHoleVertices-1;

            plint oldVertex1 = holes[iHole].boundaryVertices[oldId];
            plint oldVertex2 = holes[iHole].boundaryVertices[oldId_p1];
            Array<T,2> p1 = mesh.getVertex ( oldVertex1 );
            Array<T,2> p2 = inletPoints[newId];
            Array<T,2> p3 = mesh.getVertex ( oldVertex2 );
            Array<T,2> p4 = inletPoints[newId_p1];

            pieceSegments.push_back ( Segment ( p1,p3,p2 ) );
            pieceSegments.push_back ( Segment ( p2,p3,p4 ) );

            std::swap ( newId, newId_p1 );
            std::swap ( oldId, oldId_p1 );
        }

        if ( oddNumber )
        {
            plint id_a = holes[iHole].boundaryVertices[oldId];
            plint oldId_p2 = oldId-1;
            if ( oldId_p2<0 ) oldId_p2 = numHoleVertices-1;
            plint id_b = holes[iHole].boundaryVertices[oldId_p2];
            plint id_c = newId_p1;
            Array<T,2> a = mesh.getVertex ( id_a );
            Array<T,2> b = mesh.getVertex ( id_b );
            Array<T,2> c = inletPoints[id_c];
            pieceSegments.push_back ( Segment ( a,b,c ) );
        }


        fullGeometry.insert ( fullGeometry.end(), pieceSegments.begin(), pieceSegments.end() );
    }

    return SegmentSet<T> ( fullGeometry, geometryWithOpenings.getPrecision() );
}


template<typename T>
SegmentSet<T> constructRectangle ( T lx, T ly )
{
    std::vector<typename SegmentSet<T>::Segment> segments;

    Array< T, 2  > &pA ( 0.,0. );
    Array< T, 2  >  pB ( lx,0. );
    Array< T, 2  > &pC ( lx,ly );
    Array< T, 2  >  pD ( 0 ,ly );

    segments.push_back ( SegmentSet<T>::Segment ( pA,pB ) );
    segments.push_back ( SegmentSet<T>::Segment ( pB,pC ) );
    segments.push_back ( SegmentSet<T>::Segment ( pC,pD ) );
    segments.push_back ( SegmentSet<T>::Segment ( pD,pA ) );

    return SegmentSet<T> ( segments );
}


} // namespace plb

#endif  // SEGMENT_SET_GENERATOR_HH

