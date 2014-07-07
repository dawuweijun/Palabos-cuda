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

#ifndef SEGMENT_SET_HH
#define SEGMENT_SET_HH

#include "offLattice/segmentSet.h"
#include "core/util.h"
#include <stdlib.h>
#include <algorithm>
#include <limits>
#include <cstdio>
#include <cmath>

namespace plb
{

template<typename T>
SegmentSet<T>::SegmentSet ( Precision precision_ )
    : minSegLength ( std::numeric_limits<T>::max() ),
      maxSegLength ( std::numeric_limits<T>::min() )
{
    PLB_ASSERT ( precision_ == FLT || precision_ == DBL || precision_ == LDBL );
    precision = precision_;

    boundingCuboid.lowerLeftCorner  = Array<T,2> ( ( T ) 0.0, ( T ) 0.0 );
    boundingCuboid.upperRightCorner = Array<T,2> ( ( T ) 0.0, ( T ) 0.0 );
}

template<typename T>
SegmentSet<T>::SegmentSet ( std::vector<Segment> const& segments_, Precision precision_ )
    : segments ( segments_ ),
      minSegLength ( std::numeric_limits<T>::max() ),
      maxSegLength ( std::numeric_limits<T>::min() )
{
    PLB_ASSERT ( precision_ == FLT || precision_ == DBL || precision_ == LDBL );
    precision = precision_;

    computeMinMaxSegments();
    computeBoundingCuboid();
}

template<typename T>
SegmentSet<T>::SegmentSet ( std::string fname, Precision precision_, SurfaceGeometryFileFormat fformat )
    : minSegLength ( std::numeric_limits<T>::max() ),
      maxSegLength ( std::numeric_limits<T>::min() )
{
    PLB_ASSERT ( precision_ == FLT || precision_ == DBL || precision_ == LDBL );
    PLB_ASSERT ( fformat == STL );
    precision = precision_;

    switch ( fformat )
    {
    case STL:
    default:
        readSTL ( fname );
        break;
    }

    // TODO: Check if the next call is actually needed, since the min and max edges are
    //       computed in the readAsciiSTL and readBinarySTL functions as well.
    computeMinMaxSegments();

    computeBoundingCuboid();
}

template<typename T>
std::vector<typename SegmentSet<T>::Segment> const&
SegmentSet<T>::getSegments() const
{
    return segments;
}

template<typename T>
void SegmentSet<T>::setPrecision ( Precision precision_ )
{
    PLB_ASSERT ( precision_ == FLT || precision_ == DBL || precision_ == LDBL );
    precision = precision_;
}

template<typename T>
void SegmentSet<T>::readSTL ( std::string fname )
{
    char buf[256];
    FILE *fp = fopen ( fname.c_str(), "r" );
    PLB_ASSERT ( fp != NULL ); // The input file cannot be read.

    char *sp = fgets ( buf, 256, fp );
    PLB_ASSERT ( sp != NULL ); // The input file cannot be read.
    rewind ( fp );

    if ( strstr ( buf, "solid" ) != NULL )
    {
        readAsciiSTL ( fp );
    }
    else
    {
        readBinarySTL ( fp );
    }

    fclose ( fp );
}

template<typename T>
void SegmentSet<T>::readAsciiSTL ( FILE* fp )
{
    char buf[256];
    char *cp, *sp;

    sp = fgets ( buf, 256, fp );
    PLB_ASSERT ( sp != NULL ); // The input file is badly structured.

    char fmt[32];
    bool failed = false;
    if ( sizeof ( T ) == sizeof ( float ) )
        strcpy ( fmt, "%f%f" );
    else if ( sizeof ( T ) == sizeof ( double ) )
        strcpy ( fmt, "%lf%lf" );
    else if ( sizeof ( T ) == sizeof ( long double ) )
        strcpy ( fmt, "%Lf%Lf" );
    else
        failed = true;

    PLB_ASSERT ( !failed ); // The input file cannot be read.

    cp = strstr ( buf, "solid" );
    PLB_ASSERT ( cp != NULL ); // The input file is badly structured.

    while ( cp != NULL && !failed )
    {
        if ( fgets ( buf, 256, fp ) == NULL )
        {
            failed = true;
        }

        do
        {
            if ( ( cp = strstr ( buf, "segment normal" ) ) == NULL )
            {
                failed = true;
            }
            cp += 14;
            Array<T,2> n;
            if ( sscanf ( cp, fmt, &n[0], &n[1] ) != 2 )
            {
                failed = true;
            }

            if ( fgets ( buf, 256, fp ) == NULL || strstr ( buf, "outer loop" ) == NULL )
            {
                failed = true;
            }

            Segment segment;
            T segLength;
            for ( int i = 0; i < 2; i++ )
            {
                if ( fgets ( buf, 256, fp ) == NULL ||
                        ( cp = strstr ( buf, "vertex" ) ) == NULL )
                {
                    failed = true;
                }
                cp += 6;
                segment[i][0] = T();
                segment[i][1] = T();
                if ( sscanf ( cp, fmt,
                              &segment[i][0],
                              &segment[i][1] ) != 2 )
                {
                    failed = true;
                }
            }

            if ( fgets ( buf, 256, fp ) == NULL || strstr ( buf, "endloop" ) == NULL )
            {
                failed = true;
            }

            if ( fgets ( buf, 256, fp ) == NULL || strstr ( buf, "endsegment" ) == NULL )
            {
                failed = true;
            }

            if ( checkNoAbort ( segment, n ) )
            {
                segments.push_back ( segment );

                computeSegmentLength ( segments.size()-1, segLength );
                minSegLength = std::min ( minSegLength, segLength );
                maxSegLength = std::max ( maxSegLength, segLength );
            }

            if ( fgets ( buf, 256, fp ) == NULL )
            {
                failed = true;
            }
            cp = strstr ( buf, "endsolid" );
        }
        while ( cp == NULL && !failed );

        if ( fgets ( buf, 256, fp ) == NULL )
            break;

        cp = strstr ( buf, "solid" );
    }
    PLB_ASSERT ( !failed ); // The input file is badly structured.
}

template<typename T>
void SegmentSet<T>::readBinarySTL ( FILE* fp )
{
    char buf[256];
    unsigned int nt;
    float array[2];
    unsigned short abc;
    bool failed = false;

    int count = 0;
    while ( fread ( buf, sizeof ( char ), 80, fp ) == 80 &&
            fread ( &nt, sizeof ( unsigned int ), 1, fp ) == 1 && !failed )
    {
        count++;
        T segLength;
        for ( unsigned it = 0; it < nt && !failed; it++ )
        {
            if ( fread ( array, sizeof ( float ), 2, fp ) != 2 )
            {
                failed = true;
            }
            Array<T,2> n;
            n[0] = array[0];
            n[1] = array[1];

            Segment segment;
            for ( int i = 0; i < 2 && !failed; i++ )
            {
                if ( fread ( array, sizeof ( float ), 2, fp ) != 2 )
                {
                    failed = true;
                }
                segment[i][0] = T();
                segment[i][1] = T();
                segment[i][0] = array[0];
                segment[i][1] = array[1];
            }

            if ( fread ( &abc, sizeof ( unsigned short ), 1, fp ) != 1 )
            {
                failed = true;
            }

            if ( checkNoAbort ( segment, n ) )
            {
                segments.push_back ( segment );

                computeSegmentLength ( segments.size()-1, segLength );
                minSegLength = std::min ( minSegLength, segLength );
                maxSegLength = std::max ( maxSegLength, segLength );
            }
        }
    }

    if ( count == 0 )
        failed = true;

    PLB_ASSERT ( !failed ); // The input file is badly structured.
}

/// Make some optional checks and fix segment orientation.
template<typename T>
void SegmentSet<T>::check ( Segment& segment, Array<T,2> const& n )
{
    T eps = getEpsilon<T> ( precision );

    Array<T,2> v01 = segment[1] - segment[0];

    T segLength = sqrt ( VectorTemplateImpl<T,2>::normSqr ( v01 ) );

    if ( util::fpequal ( segLength, ( T ) 0.0, eps ) )
    {
        PLB_ASSERT ( false ); // The length of the segment too short.
    }

    if ( v01[0]*n[1]-v01[1]*n[0] < ( T ) 0.0 )
    {
        std::swap ( segment[0],segment[1] );
    }
}

/// Make some optional checks and fix segment orientation.
template<typename T>
bool SegmentSet<T>::checkNoAbort ( Segment& segment, Array<T,2> const& n )
{
    T eps = getEpsilon<T> ( precision );

    Array<T,2> v01 = segment[1] - segment[0];

    T segLength = sqrt ( VectorTemplateImpl<T,2>::normSqr ( v01 ) );

    if ( util::fpequal ( segLength, ( T ) 0.0, eps ) )
    {
        return false; // The length of the segment too short.
    }

    if ( v01[0]*n[1]-v01[1]*n[0] < ( T ) 0.0 )
    {
        std::swap ( segment[0],segment[1] );
    }
    return true;
}

template<typename T>
void SegmentSet<T>::translate ( Array<T,2> const& vector )
{
    T eps = std::numeric_limits<T>::epsilon();


    if ( util::fpequal ( ( T ) sqrt ( VectorTemplateImpl<T,2>::normSqr ( vector ) ), ( T ) 0.0, eps ) )
        return;

    plint size = segments.size();
    if ( size == 0 )
        return;

    for ( plint i = 0; i < size; i++ )
    {
        segments[i][0] += vector;
        segments[i][1] += vector;
    }

    boundingCuboid.lowerLeftCorner  += vector;
    boundingCuboid.upperRightCorner += vector;
}

template<typename T>
void SegmentSet<T>::scale ( T alpha )
{
    T eps = std::numeric_limits<T>::epsilon();
    if ( util::fpequal ( alpha, ( T ) 1.0, eps ) )
        return;

    plint size = segments.size();
    if ( size == 0 )
        return;

    for ( plint i = 0; i < size; i++ )
    {
        for ( int j = 0; j < 2; j++ )
        {
            segments[i][j] *= alpha;
        }
    }
    minSegLength *= alpha;
    maxSegLength *= alpha;

    boundingCuboid.lowerLeftCorner  *= alpha;
    boundingCuboid.upperRightCorner *= alpha;
}

template<typename T>
void SegmentSet<T>::rotate ( T theta )
{
    T eps = std::numeric_limits<T>::epsilon();

    plint size = segments.size();
    if ( size == 0 )
        return;

    T a[2][2];
    a[0][0] =  cos ( theta );
    a[0][1] = -sin ( theta );
    a[1][0] =  sin ( theta );
    a[1][1] =  cos ( theta );


    for ( plint iSegment = 0; iSegment < size; iSegment++ )
    {
        for ( int iVertex = 0; iVertex < 2; iVertex++ )
        {
            Array<T,2> x = segments[iSegment][iVertex];
            for ( int i = 0; i < 2; i++ )
            {
                segments[iSegment][iVertex][i] = ( T ) 0.0;
                for ( int j = 0; j < 2; j++ )
                {
                    segments[iSegment][iVertex][i] += a[i][j]*x[j];
                }
            }
        }
    }

    computeBoundingCuboid();
}

template<typename T>
void SegmentSet<T>::merge ( std::vector<SegmentSet<T>*> meshes )
{
    PLB_ASSERT ( meshes.size() != 0 );

    segments.assign ( meshes[0]->getSegments().begin(), meshes[0]->getSegments().end() );
    minSegLength = meshes[0]->getMinEdgeLength();
    maxSegLength = meshes[0]->getMaxEdgeLength();
    boundingCuboid = meshes[0]->getBoundingCuboid();
    for ( pluint i = 1; i < meshes.size(); i++ )
    {
        segments.insert ( segments.end(), meshes[i]->getSegments().begin(), meshes[i]->getSegments().end() );
        minSegLength = std::min ( minSegLength, meshes[i]->getMinEdgeLength() );
        maxSegLength = std::max ( maxSegLength, meshes[i]->getMaxEdgeLength() );

        Cuboid2D<T> bcuboid = meshes[i]->getBoundingCuboid();
        for ( plint j = 0; j < 2; j++ )
        {
            boundingCuboid.lowerLeftCorner[j]  = std::min ( boundingCuboid.lowerLeftCorner[j],
                                                 bcuboid.lowerLeftCorner[j] );
            boundingCuboid.upperRightCorner[j] = std::max ( boundingCuboid.upperRightCorner[j],
                                                 bcuboid.upperRightCorner[j] );
        }
    }
}

template<typename T>
void SegmentSet<T>::append ( SegmentSet<T> const& mesh )
{
    segments.insert ( segments.end(), mesh.getSegments().begin(), mesh.getSegments().end() );
    minSegLength = std::min ( minSegLength, mesh.getMinEdgeLength() );
    maxSegLength = std::max ( maxSegLength, mesh.getMaxEdgeLength() );

    Cuboid2D<T> bcuboid = mesh.getBoundingCuboid();
    for ( plint j = 0; j < 2; j++ )
    {
        boundingCuboid.lowerLeftCorner[j]  = std::min ( boundingCuboid.lowerLeftCorner[j],
                                             bcuboid.lowerLeftCorner[j] );
        boundingCuboid.upperRightCorner[j] = std::max ( boundingCuboid.upperRightCorner[j],
                                             bcuboid.upperRightCorner[j] );
    }
}

template<typename T>
void SegmentSet<T>::refine()
{
    std::vector<Segment> newSegments;

    for ( pluint i = 0; i < segments.size(); ++i )
    {
        Segment const& segment = segments[i];

        Array<T,2> v00 = segment[0];
        Array<T,2> v11 = segment[1];

        Array<T,2> v01 = 0.5 * ( v00 + v11 );

        Segment newSegment;

        newSegment[0] = v00;
        newSegment[1] = v01;
        newSegments.push_back ( newSegment );

        newSegment[0] = v01;
        newSegment[1] = v11;
        newSegments.push_back ( newSegment );
    }

    segments.clear();
    segments = newSegments;

    computeMinMaxSegments();
    computeBoundingCuboid();
}

template<typename T>
void SegmentSet<T>::reverseOrientation()
{
    plint size = segments.size();
    for ( plint i = 0; i < size; i++ )
        std::swap ( segments[i][0],segments[i][1] );
}

template<typename T>
void SegmentSet<T>::writeAsciiSTL ( std::string fname ) const
{
    if ( global::mpi().isMainProcessor() )
    {
        FILE *fp = fopen ( fname.c_str(), "w" );
        PLB_ASSERT ( fp != NULL );

        plint size = segments.size();
        fprintf ( fp, "solid plb\n" );
        for ( plint i = 0; i < size; i++ )
        {
            Array<T,2> v0 = segments[i][0];
            Array<T,2> v1 = segments[i][1];

            Array<T,2> e01 = v1 - v0;
            //顺着线段的正方向看，左侧为正，右侧为负
            Array<T,2> n ( -e01[1],e01[0] );
            n /= sqrt ( VectorTemplateImpl<T,2>::normSqr ( n ) );
            fprintf ( fp, "  segment normal % e % e\n", ( double ) n[0], ( double ) n[1] );
            fprintf ( fp, "    outer loop\n" );
            fprintf ( fp, "      vertex % e % e\n", ( double ) v0[0], ( double ) v0[1] );
            fprintf ( fp, "      vertex % e % e\n", ( double ) v1[0], ( double ) v1[1] );
            fprintf ( fp, "    endloop\n" );
            fprintf ( fp, "  endsegment\n" );
        }
        fprintf ( fp, "endsolid plb\n" );
        fclose ( fp );
    }
}

template<typename T>
void SegmentSet<T>::writeBinarySTL ( std::string fname ) const
{
    if ( global::mpi().isMainProcessor() )
    {
        FILE *fp = fopen ( fname.c_str(), "wb" );
        PLB_ASSERT ( fp != NULL );

        unsigned int nt = segments.size();
        unsigned short abc = 0;
        char buf[80];

        for ( int i = 0; i < 80; i++ )
            buf[i] = '\0';

        fwrite ( buf, sizeof ( char ), 80, fp );
        fwrite ( &nt, sizeof ( unsigned int ), 1, fp );
        for ( unsigned int i = 0; i < nt; i++ )
        {
            Array<T,2> v0 = segments[i][0];
            Array<T,2> v1 = segments[i][1];

            Array<T,2> e01 = v1 - v0;
            //顺着线段的正方向看，左侧为正，右侧为负
            Array<T,2> nrml ( -e01[1],e01[0] );
            nrml /= sqrt ( VectorTemplateImpl<T,2>::normSqr ( nrml ) );

            float n[2];
            n[0] = nrml[0];
            n[1] = nrml[1];
            fwrite ( ( void * ) n, sizeof ( float ), 2, fp );
            float v[2];
            v[0] = v0[0];
            v[1] = v0[1];
            fwrite ( ( void * ) v, sizeof ( float ), 2, fp );
            v[0] = v1[0];
            v[1] = v1[1];
            fwrite ( ( void * ) v, sizeof ( float ), 2, fp );
            fwrite ( &abc, sizeof ( unsigned short ), 1, fp );
        }

        fclose ( fp );
    }
}

template<typename T>
void SegmentSet<T>::writeSVG ( std::string fname ) const
{
    if ( global::mpi().isMainProcessor() )
    {
        plb_ofstream out ( fname.c_str() );

        plint size = segments.size();

        //write header
        out << "<?xml version=\"1.0\" standalone=\"no\"?>\n";
        out << "<!DOCTYPE svg PUBLIC"
            << " \"-//W3C//DTD SVG 1.1//EN\" \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n";
        out << "<svg width=\"100%\" height=\"100%\" version=\"1.1\" xmlns=\"http://www.w3.org/2000/svg\">\n";
        //line set
        for ( plint i = 0; i < size; i++ )
        {
            out<<"<line x1=\""<<segments[i][0][0]
		    <<"\" y1=\""<<segments[i][0][1]
		    <<"\" x2=\""<<segments[i][1][0]
		    <<"\" y2=\""<<segments[i][1][1]
		    <<"\" style=\"stroke:rgb(0,0,0);stroke-width:0.5\"/>\n";
        }
        //write end
        out<<"</svg>";
    }
}

template<typename T>
int SegmentSet<T>::cutSegmentWithLine ( Line2D<T> const& line, Segment const& segment,
                                        SegmentSet<T>& newSegmentSet ) const
{
    T epsilon = getEpsilon<T> ( precision );

    int vertexTags[2];

    // Tag the segment vertices.
    for ( int iVertex = 0; iVertex < 2; iVertex++ )
    {
        Array<T,2> tmp = segment[iVertex] - line.point;
        T norm_tmp = norm ( tmp );
        if ( norm_tmp > epsilon )
        {
            tmp /= norm_tmp;
        }
        else
        {
            tmp[0] = tmp[1] = ( T ) 0.0;
        }
        T dotp = dot ( tmp, line.normal );
        if ( fabs ( dotp ) <= epsilon )
        {
            vertexTags[iVertex] = 0;
        }
        else if ( dotp > ( T ) 0.0 && fabs ( dotp ) > epsilon )
        {
            vertexTags[iVertex] = -1;
        }
        else if ( dotp < ( T ) 0.0 && fabs ( dotp ) > epsilon )
        {
            vertexTags[iVertex] = 1;
        }
        else
        {
            return -1;
        }
    }

    // All three vertices belong to one side of the cut line.
    if ( vertexTags[0] == 1 && vertexTags[1] == 1 )
    {
        newSegmentSet.segments.push_back ( segment );
        return 1;
    }
    else if ( vertexTags[0] == -1 && vertexTags[1] == -1 )
    {
        return 0;
    }

    // One vertex belongs to one side of the cut line and the other one vertice
    //   belong to the other side.
    if ( vertexTags[0] == 1 && vertexTags[1] == -1 )
    {
        Array<T,2> intersection_ij ( ( T ) 0.0, ( T ) 0.0 );
        int rv = 0;
        rv = lineIntersectionWithLine<T> ( line, segment[0], segment[1], precision, intersection_ij );
        if ( rv != 1 )
        {
            return -1;
        }
        Segment newSegment ( segment[0], intersection_ij );
        newSegmentSet.segments.push_back ( newSegment );
        return 1;
    }
    else if ( vertexTags[0] == -1 && vertexTags[1] == 1 )
    {
        Array<T,2> intersection_ij ( ( T ) 0.0, ( T ) 0.0 );
        int rv = 0;
        rv = lineIntersectionWithLine<T> ( line, segment[0], segment[1], precision, intersection_ij );
        if ( rv != 1 )
        {
            return -1;
        }
        Segment newSegment_0 ( intersection_ij, segment[1] );
        newSegmentSet.segments.push_back ( newSegment_0 );
        return 1;
    }

    // Only one vertex belongs to the cut line.
    if ( vertexTags[0]+vertexTags[1] == 1 )
    {
        newSegmentSet.segments.push_back ( segment );
        return 1;
    }
    else
    {
        return 0;
    }


// All 2 vertices belong to the cut line.
    if ( vertexTags[0] == 0 && vertexTags[1] == 0 )
    {
        newSegmentSet.segments.push_back ( segment );
        return 1;
    }

    return -1;
}

template<typename T>
int SegmentSet<T>::cutWithLine ( Line2D<T> const& line, SegmentSet<T>& newSegmentSet ) const
{
    T epsilon = getEpsilon<T> ( precision );

    T norm_normal = norm ( line.normal );
    PLB_ASSERT ( norm_normal > epsilon ); // The cut line normal vector cannot have zero magnitude.
    Line2D<T> newLine;
    newLine.point = line.point;
    newLine.normal = line.normal / norm_normal;

    newSegmentSet.segments.resize ( 0 );

    newSegmentSet.precision = precision;

    for ( pluint iSegment = 0; iSegment < segments.size(); iSegment++ )
    {
        if ( cutSegmentWithLine ( newLine, segments[iSegment], newSegmentSet ) == -1 )
        {
            return -1;
        }
    }

    if ( newSegmentSet.segments.size() != 0 )
    {
        newSegmentSet.computeMinMaxEdges();
        newSegmentSet.computeBoundingCuboid();
    }

    if ( newSegmentSet.segments.size() == 0 || newSegmentSet.segments.size() == segments.size() )
    {
        return 0;
    }

    return 1;
}

template<typename T>
int SegmentSet<T>::cutWithLine (
    Line2D<T> const& line, Cuboid2D<T> const& cuboid, SegmentSet<T>& newSegmentSet ) const
{
    T epsilon = getEpsilon<T> ( precision );

    T norm_normal = norm ( line.normal );
    PLB_ASSERT ( norm_normal > epsilon ); // The cut line normal vector cannot have zero magnitude.
    Line2D<T> newLine;
    newLine.point = line.point;
    newLine.normal = line.normal / norm_normal;

    T norm_diagonal = norm ( cuboid.upperRightCorner - cuboid.lowerLeftCorner );
    PLB_ASSERT ( norm_diagonal > epsilon ); // The diagonal of the cuboid cannot have zero length.

    newSegmentSet.segments.resize ( 0 );

    newSegmentSet.precision = precision;

    for ( pluint iSegment = 0; iSegment < segments.size(); iSegment++ )
    {
        Segment const& segment = segments[iSegment];

        Array<T,2> vertices[2];
        vertices[0] = segment[0];
        vertices[1] = segment[1];

        // Check if the segment is fully contained in the cuboid.
        int isNotFullyContained = 0;
        for ( int iVertex = 0; iVertex < 2; iVertex++ )
        {
            Array<T,2> diff_l;
            diff_l = vertices[iVertex] - cuboid.lowerLeftCorner;

            Array<T,2> diff_u;
            diff_u = vertices[iVertex] - cuboid.upperRightCorner;

            if ( ( diff_l[0] < ( T ) 0.0 && fabs ( diff_l[0] ) > epsilon ) ||
                    ( diff_l[1] < ( T ) 0.0 && fabs ( diff_l[1] ) > epsilon ) ||
                    ( diff_u[0] > ( T ) 0.0 && fabs ( diff_u[0] ) > epsilon ) ||
                    ( diff_u[1] > ( T ) 0.0 && fabs ( diff_u[1] ) > epsilon ) )
            {
                isNotFullyContained = 1;
                break;
            }
        }

        if ( isNotFullyContained )
        {
            newSegmentSet.segments.push_back ( segment );
            continue;
        }

        if ( cutSegmentWithLine ( newLine, segment, newSegmentSet ) == -1 )
            return -1;
    }

    if ( newSegmentSet.segments.size() != 0 )
    {
        newSegmentSet.computeMinMaxSegments();
        newSegmentSet.computeBoundingCuboid();
    }

    if ( newSegmentSet.segments.size() == 0 || newSegmentSet.segments.size() == segments.size() )
    {
        return 0;
    }

    return 1;
}

template<typename T>
void SegmentSet<T>::computeMinMaxSegments()
{
    T segLength;
    for ( pluint i=0; i<segments.size(); ++i )
    {
        computeSegmentLength ( i, segLength );
        minSegLength = std::min ( minSegLength, segLength );
        maxSegLength = std::max ( maxSegLength, segLength );
    }
}

template<typename T>
void SegmentSet<T>::computeSegmentLength ( plb::pluint iSegment, T& segLength ) const
{
    PLB_ASSERT ( iSegment<segments.size() );
    Segment const& segment = segments[iSegment];
    segLength = norm ( segment[1]-segment[0] );
}

template<typename T>
void SegmentSet<T>::computeBoundingCuboid()
{
    T xMin, yMin;
    T xMax, yMax;

    xMin = yMin  =  std::numeric_limits<T>::max();
    xMax = yMax  = -std::numeric_limits<T>::max();
    for ( pluint i=0; i<segments.size(); ++i )
    {
        Segment const& segment = segments[i];

        xMin = std::min ( xMin, std::min ( segment[0][0], segment[1][0] ) );
        yMin = std::min ( yMin, std::min ( segment[0][1], segment[1][1] ) );

        xMax = std::max ( xMax, std::max ( segment[0][0], segment[1][0] ) );
        yMax = std::max ( yMax, std::max ( segment[0][1], segment[1][1] ) );
    }
    boundingCuboid.lowerLeftCorner  = Array<T,2> ( xMin, yMin );
    boundingCuboid.upperRightCorner = Array<T,2> ( xMax, yMax );
}

} // namespace plb

#endif  // SEGMENT_SET_HH




