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
#include <algorithm>
#include <limits>
#include <cstdio>
#include <cmath>

namespace plb {

template<typename T>
SegmentSet<T>::SegmentSet(Precision precision_)
    : minEdgeLength(std::numeric_limits<T>::max()),
      maxEdgeLength(std::numeric_limits<T>::min())
{
    PLB_ASSERT(precision_ == FLT || precision_ == DBL || precision_ == LDBL);
    precision = precision_;

    boundingCuboid.lowerLeftCorner  = Array<T,2>((T) 0.0, (T) 0.0, (T) 0.0);
    boundingCuboid.upperRightCorner = Array<T,2>((T) 0.0, (T) 0.0, (T) 0.0);
}

template<typename T>
SegmentSet<T>::SegmentSet(std::vector<Segment> const& segments_, Precision precision_)
    : segments(segments_),
      minEdgeLength(std::numeric_limits<T>::max()),
      maxEdgeLength(std::numeric_limits<T>::min())
{
    PLB_ASSERT(precision_ == FLT || precision_ == DBL || precision_ == LDBL);
    precision = precision_;

    computeMinMaxEdges();
    computeBoundingCuboid();
}

template<typename T>
SegmentSet<T>::SegmentSet(std::string fname, Precision precision_, SurfaceGeometryFileFormat fformat)
    : minEdgeLength(std::numeric_limits<T>::max()),
      maxEdgeLength(std::numeric_limits<T>::min())
{
    PLB_ASSERT(precision_ == FLT || precision_ == DBL || precision_ == LDBL);
    PLB_ASSERT(fformat == STL);
    precision = precision_;

    switch (fformat) {
    case STL: default:
        readSTL(fname);
        break;
    }

    // TODO: Check if the next call is actually needed, since the min and max edges are
    //       computed in the readAsciiSTL and readBinarySTL functions as well.
    computeMinMaxEdges();

    computeBoundingCuboid();
}

template<typename T>
std::vector<typename SegmentSet<T>::Segment> const&
    SegmentSet<T>::getSegments() const
{
    return segments;
}

template<typename T>
void SegmentSet<T>::setPrecision(Precision precision_)
{
    PLB_ASSERT(precision_ == FLT || precision_ == DBL || precision_ == LDBL);
    precision = precision_;
}

template<typename T>
void SegmentSet<T>::readSTL(std::string fname)
{
    char buf[256];
    FILE *fp = fopen(fname.c_str(), "r");
    PLB_ASSERT(fp != NULL); // The input file cannot be read.

    char *sp = fgets(buf, 256, fp);
    PLB_ASSERT(sp != NULL); // The input file cannot be read.
    rewind(fp);

    if (strstr(buf, "solid") != NULL) {
        readAsciiSTL(fp);
    }
    else {
        readBinarySTL(fp);
    }

    fclose(fp);
}

template<typename T>
void SegmentSet<T>::readAsciiSTL(FILE* fp) {
    char buf[256];
    char *cp, *sp;

    sp = fgets(buf, 256, fp);
    PLB_ASSERT(sp != NULL); // The input file is badly structured.

    char fmt[32];
    bool failed = false;
    if (sizeof(T) == sizeof(float))
        strcpy(fmt, "%f%f%f");
    else if (sizeof(T) == sizeof(double))
        strcpy(fmt, "%lf%lf%lf");
    else if (sizeof(T) == sizeof(long double))
        strcpy(fmt, "%Lf%Lf%Lf");
    else
        failed = true;

    PLB_ASSERT(!failed); // The input file cannot be read.

    cp = strstr(buf, "solid");
    PLB_ASSERT(cp != NULL); // The input file is badly structured.

    while (cp != NULL && !failed) {
        if (fgets(buf, 256, fp) == NULL) {
            failed = true;
        }

        do {
            if ((cp = strstr(buf, "facet normal")) == NULL) {
                failed = true;
            }
            cp += 12;
            Array<T,2> n;
            if (sscanf(cp, fmt, &n[0], &n[1], &n[2]) != 3) {
                failed = true;
            }

            if (fgets(buf, 256, fp) == NULL || strstr(buf, "outer loop") == NULL)
            {
                failed = true;
            }

            Segment segment;
            T nextMin, nextMax;
            for (int i = 0; i < 3; i++) {
                if ( fgets(buf, 256, fp) == NULL ||
                     (cp = strstr(buf, "vertex")) == NULL )
                {
                    failed = true;
                }
                cp += 6;
                segment[i][0] = T();
                segment[i][1] = T();
                segment[i][2] = T();
                if (sscanf( cp, fmt,
                            &segment[i][0],
                            &segment[i][1],
                            &segment[i][2] ) != 3)
                {
                    failed = true;
                }
            }

            if (fgets(buf, 256, fp) == NULL || strstr(buf, "endloop") == NULL) {
                failed = true;
            }

            if (fgets(buf, 256, fp) == NULL || strstr(buf, "endfacet") == NULL) {
                failed = true;
            }

            if (checkNoAbort(segment, n)) {
                segments.push_back(segment);

                computeMinMaxEdge(segments.size()-1, nextMin, nextMax);
                minEdgeLength = std::min(minEdgeLength, nextMin);
                maxEdgeLength = std::max(maxEdgeLength, nextMax);
            }

            if (fgets(buf, 256, fp) == NULL) {
                failed = true;
            }
            cp = strstr(buf, "endsolid");
        }
        while (cp == NULL && !failed);

        if (fgets(buf, 256, fp) == NULL)
            break;

        cp = strstr(buf, "solid");
    }
    PLB_ASSERT(!failed); // The input file is badly structured.
}

template<typename T>
void SegmentSet<T>::readBinarySTL(FILE* fp)
{
    char buf[256];
    unsigned int nt;
    float array[3];
    unsigned short abc;
    bool failed = false;

    int count = 0;
    while (fread(buf, sizeof(char), 80, fp) == 80 &&
           fread(&nt, sizeof(unsigned int), 1, fp) == 1 && !failed) {
        count++;
        T nextMin, nextMax;
        for (unsigned it = 0; it < nt && !failed; it++) {
            if (fread(array, sizeof(float), 3, fp) != 3) {
                failed = true;
            }
            Array<T,2> n;
            n[0] = array[0];
            n[1] = array[1];
            n[2] = array[2];

            Segment segment;
            for (int i = 0; i < 3 && !failed; i++) {
                if (fread(array, sizeof(float), 3, fp) != 3) {
                    failed = true;
                }
                segment[i][0] = T();
                segment[i][1] = T();
                segment[i][2] = T();
                segment[i][0] = array[0];
                segment[i][1] = array[1];
                segment[i][2] = array[2];
            }

            if (fread(&abc, sizeof(unsigned short), 1, fp) != 1) {
                failed = true;
            }

            if (checkNoAbort(segment, n)) {
                segments.push_back(segment);

                computeMinMaxEdge(segments.size()-1, nextMin, nextMax);
                minEdgeLength = std::min(minEdgeLength, nextMin);
                maxEdgeLength = std::max(maxEdgeLength, nextMax);
            }
        }
    }

    if (count == 0) 
        failed = true;

    PLB_ASSERT(!failed); // The input file is badly structured.
}

/// Make some optional checks and fix segment orientation.
template<typename T>
void SegmentSet<T>::check(Segment& segment, Array<T,2> const& n)
{
    T eps = getEpsilon<T>(precision);

    Array<T,2> v01 = segment[1] - segment[0];
    Array<T,2> v02 = segment[2] - segment[0];
    Array<T,2> v12 = segment[2] - segment[1];

    T norm01 = sqrt(VectorTemplateImpl<T,2>::normSqr(v01));
    T norm02 = sqrt(VectorTemplateImpl<T,2>::normSqr(v02));
    T norm12 = sqrt(VectorTemplateImpl<T,2>::normSqr(v12));

    if (util::fpequal(norm01, (T) 0.0, eps) || util::fpequal(norm02, (T) 0.0, eps) ||
        util::fpequal(norm12, (T) 0.0, eps))
    {
        PLB_ASSERT(false); // One of the segments is degenerate.
    }

    Array<T,2> cn;
    crossProduct(v01, v02, cn);

    T norm_cn = sqrt(VectorTemplateImpl<T,2>::normSqr(cn));

    if (util::fpequal(norm_cn, (T) 0.0, eps))
    {
        PLB_ASSERT(false); // One of the segments has zero area.
    }

    T dot = VectorTemplateImpl<T,2>::scalarProduct(cn,n);

    if (dot < (T) 0.0) {
        std::swap(segment[1],segment[2]);
    }
}

/// Make some optional checks and fix segment orientation.
template<typename T>
bool SegmentSet<T>::checkNoAbort(Segment& segment, Array<T,2> const& n)
{
    T eps = getEpsilon<T>(precision);

    Array<T,2> v01 = segment[1] - segment[0];
    Array<T,2> v02 = segment[2] - segment[0];
    Array<T,2> v12 = segment[2] - segment[1];

    T norm01 = sqrt(VectorTemplateImpl<T,2>::normSqr(v01));
    T norm02 = sqrt(VectorTemplateImpl<T,2>::normSqr(v02));
    T norm12 = sqrt(VectorTemplateImpl<T,2>::normSqr(v12));

    if (util::fpequal(norm01, (T) 0.0, eps) || util::fpequal(norm02, (T) 0.0, eps) ||
        util::fpequal(norm12, (T) 0.0, eps))
    {
        return false; // One of the segments is degenerate.
    }

    Array<T,2> cn;
    crossProduct(v01, v02, cn);

    T norm_cn = sqrt(VectorTemplateImpl<T,2>::normSqr(cn));

    if (util::fpequal(norm_cn, (T) 0.0, eps))
    {
        return false; // One of the segments has zero area.
    }

    T dot = VectorTemplateImpl<T,2>::scalarProduct(cn,n);

    if (dot < (T) 0.0) {
        std::swap(segment[1],segment[2]);
    }
    return true;
}

template<typename T>
void SegmentSet<T>::translate(Array<T,2> const& vector)
{
    T eps = std::numeric_limits<T>::epsilon();
    if (util::fpequal((T) sqrt(VectorTemplateImpl<T,2>::normSqr(vector)), (T) 0.0, eps))
        return;

    plint size = segments.size();
    if (size == 0)
        return;

    for (plint i = 0; i < size; i++) {
        for (int j = 0; j < 3; j++) {
            segments[i][j] += vector;
        }
    }

    boundingCuboid.lowerLeftCorner  += vector;
    boundingCuboid.upperRightCorner += vector;
}

template<typename T>
void SegmentSet<T>::scale(T alpha)
{
    T eps = std::numeric_limits<T>::epsilon();
    if (util::fpequal(alpha, (T) 1.0, eps))
        return;

    plint size = segments.size();
    if (size == 0)
        return;

    for (plint i = 0; i < size; i++) {
        for (int j = 0; j < 3; j++) {
            segments[i][j] *= alpha;
        }
    }
    minEdgeLength *= alpha;
    maxEdgeLength *= alpha;

    boundingCuboid.lowerLeftCorner  *= alpha;
    boundingCuboid.upperRightCorner *= alpha;
}

template<typename T>
void SegmentSet<T>::rotateAtOrigin(Array<T,2> const& normedAxis, T theta) {
    plint size = (plint)segments.size();
    for (plint iSegment = 0; iSegment < size; iSegment++) {
        for (int iVertex = 0; iVertex < 3; iVertex++) {
            segments[iSegment][iVertex] =
                plb::rotateAtOrigin (
                    segments[iSegment][iVertex], normedAxis, theta );
        }
    }

    computeBoundingCuboid();
}

template<typename T>
void SegmentSet<T>::rotate(T phi, T theta, T psi)
{   
    T eps = std::numeric_limits<T>::epsilon();
    T pi = acos(-1.0);

    PLB_ASSERT((theta > (T) 0.0 || util::fpequal(theta, (T) 0.0, eps)) &&
               (theta < pi  || util::fpequal(theta, pi, eps)));

    plint size = segments.size();
    if (size == 0)
        return;

    T a[3][3];
    a[0][0] =  (T) 1.0;
    a[0][1] =  (T) 0.0;
    a[0][2] =  (T) 0.0;
    a[1][0] =  (T) 0.0;
    a[1][1] =  cos(theta);
    a[1][2] = -sin(theta);
    a[2][0] =  (T) 0.0;
    a[2][1] =  sin(theta);
    a[2][2] =  cos(theta);

    T b[3][3];
    b[0][0] =  cos(phi);
    b[0][1] = -sin(phi);
    b[0][2] =  (T) 0.0;
    b[1][0] =  sin(phi);
    b[1][1] =  cos(phi);
    b[1][2] =  (T) 0.0;
    b[2][0] =  (T) 0.0;
    b[2][1] =  (T) 0.0;
    b[2][2] =  (T) 1.0;

    T c[3][3];
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            c[i][j] = (T) 0.0;
            for (int k = 0; k < 3; k++) {
                c[i][j] += a[i][k]*b[k][j];
            }
        }
    }

    b[0][0] =  cos(psi);
    b[0][1] = -sin(psi);
    b[0][2] =  (T) 0.0;
    b[1][0] =  sin(psi);
    b[1][1] =  cos(psi);
    b[1][2] =  (T) 0.0;
    b[2][0] =  (T) 0.0;
    b[2][1] =  (T) 0.0;
    b[2][2] =  (T) 1.0;

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            a[i][j] = (T) 0.0;
            for (int k = 0; k < 3; k++) {
                a[i][j] += b[i][k]*c[k][j];
            }
        }
    }

    for (plint iSegment = 0; iSegment < size; iSegment++) {
        for (int iVertex = 0; iVertex < 3; iVertex++) {
            Array<T,2> x = segments[iSegment][iVertex];
            for (int i = 0; i < 3; i++) {
                segments[iSegment][iVertex][i] = (T) 0.0;
                for (int j = 0; j < 3; j++) {
                    segments[iSegment][iVertex][i] += a[i][j]*x[j];
                }
            }
        }
    }

    computeBoundingCuboid();
}

template<typename T>
void SegmentSet<T>::merge(std::vector<SegmentSet<T>*> meshes)
{
    PLB_ASSERT(meshes.size() != 0);

    segments.assign(meshes[0]->getSegments().begin(), meshes[0]->getSegments().end());
    minEdgeLength = meshes[0]->getMinEdgeLength();
    maxEdgeLength = meshes[0]->getMaxEdgeLength();
    boundingCuboid = meshes[0]->getBoundingCuboid();
    for (pluint i = 1; i < meshes.size(); i++) {
        segments.insert(segments.end(), meshes[i]->getSegments().begin(), meshes[i]->getSegments().end());
        minEdgeLength = std::min(minEdgeLength, meshes[i]->getMinEdgeLength());
        maxEdgeLength = std::max(maxEdgeLength, meshes[i]->getMaxEdgeLength());

        Cuboid2D<T> bcuboid = meshes[i]->getBoundingCuboid();
        for (plint j = 0; j < 3; j++) {
            boundingCuboid.lowerLeftCorner[j]  = std::min(boundingCuboid.lowerLeftCorner[j],
                                                          bcuboid.lowerLeftCorner[j]);
            boundingCuboid.upperRightCorner[j] = std::max(boundingCuboid.upperRightCorner[j],
                                                          bcuboid.upperRightCorner[j]);
        }
    }
}

template<typename T>
void SegmentSet<T>::append(SegmentSet<T> const& mesh)
{
    segments.insert(segments.end(), mesh.getSegments().begin(), mesh.getSegments().end());
    minEdgeLength = std::min(minEdgeLength, mesh.getMinEdgeLength());
    maxEdgeLength = std::max(maxEdgeLength, mesh.getMaxEdgeLength());

    Cuboid2D<T> bcuboid = mesh.getBoundingCuboid();
    for (plint j = 0; j < 2; j++) {
        boundingCuboid.lowerLeftCorner[j]  = std::min(boundingCuboid.lowerLeftCorner[j],
                                                      bcuboid.lowerLeftCorner[j]);
        boundingCuboid.upperRightCorner[j] = std::max(boundingCuboid.upperRightCorner[j],
                                                      bcuboid.upperRightCorner[j]);
    }
}

template<typename T>
void SegmentSet<T>::refine()
{
    std::vector<Segment> newSegments;

    for (pluint i = 0; i < segments.size(); ++i) {
        Segment const& segment = segments[i];

        Array<T,2> v00 = segment[0];
        Array<T,2> v11 = segment[1];
        Array<T,2> v22 = segment[2];

        Array<T,2> v01 = 0.5 * (v00 + v11);
        Array<T,2> v12 = 0.5 * (v11 + v22);
        Array<T,2> v20 = 0.5 * (v22 + v00);

        Segment newSegment;

        newSegment[0] = v01;
        newSegment[1] = v12;
        newSegment[2] = v20;
        newSegments.push_back(newSegment);

        newSegment[0] = v00;
        newSegment[1] = v01;
        newSegment[2] = v20;
        newSegments.push_back(newSegment);

        newSegment[0] = v01;
        newSegment[1] = v11;
        newSegment[2] = v12;
        newSegments.push_back(newSegment);

        newSegment[0] = v20;
        newSegment[1] = v12;
        newSegment[2] = v22;
        newSegments.push_back(newSegment);
    }

    segments.clear();
    segments = newSegments;

    computeMinMaxEdges();
    computeBoundingCuboid();
}

template<typename T>
void SegmentSet<T>::reverseOrientation()
{
    plint size = segments.size();
    for (plint i = 0; i < size; i++) 
        std::swap(segments[i][1],segments[i][2]);
}

template<typename T>
void SegmentSet<T>::writeAsciiSTL(std::string fname) const
{
    if (global::mpi().isMainProcessor()) {
        FILE *fp = fopen(fname.c_str(), "w");
        PLB_ASSERT(fp != NULL);

        plint size = segments.size();
        fprintf(fp, "solid plb\n");
        for (plint i = 0; i < size; i++) {
            Array<T,2> v0 = segments[i][0];
            Array<T,2> v1 = segments[i][1];
            Array<T,2> v2 = segments[i][2];

            Array<T,2> e01 = v1 - v0;
            Array<T,2> e02 = v2 - v0;

            Array<T,2> n;
            crossProduct(e01, e02, n);
            n /= sqrt(VectorTemplateImpl<T,2>::normSqr(n));
            fprintf(fp, "  facet normal % e % e % e\n", (double) n[0], (double) n[1], (double) n[2]);
            fprintf(fp, "    outer loop\n");
            fprintf(fp, "      vertex % e % e % e\n", (double) v0[0], (double) v0[1], (double) v0[2]);
            fprintf(fp, "      vertex % e % e % e\n", (double) v1[0], (double) v1[1], (double) v1[2]);
            fprintf(fp, "      vertex % e % e % e\n", (double) v2[0], (double) v2[1], (double) v2[2]);
            fprintf(fp, "    endloop\n");
            fprintf(fp, "  endfacet\n");
        }
        fprintf(fp, "endsolid plb\n");
        fclose(fp);
    }
}

template<typename T>
void SegmentSet<T>::writeBinarySTL(std::string fname) const
{
    if (global::mpi().isMainProcessor()) {
        FILE *fp = fopen(fname.c_str(), "wb");
        PLB_ASSERT(fp != NULL);

        unsigned int nt = segments.size();
        unsigned short abc = 0;
        char buf[80];

        for (int i = 0; i < 80; i++)
            buf[i] = '\0';

        fwrite(buf, sizeof(char), 80, fp);
        fwrite(&nt, sizeof(unsigned int), 1, fp);
        for (unsigned int i = 0; i < nt; i++) {
            Array<T,2> v0 = segments[i][0];
            Array<T,2> v1 = segments[i][1];
            Array<T,2> v2 = segments[i][2];

            Array<T,2> e01 = v1 - v0;
            Array<T,2> e02 = v2 - v0;

            Array<T,2> nrml;
            crossProduct(e01, e02, nrml);
            nrml /= sqrt(VectorTemplateImpl<T,2>::normSqr(nrml));

            float n[3];
            n[0] = nrml[0];
            n[1] = nrml[1];
            n[2] = nrml[2];
            fwrite((void *) n, sizeof(float), 3, fp);
            float v[3];
            v[0] = v0[0];
            v[1] = v0[1];
            v[2] = v0[2];
            fwrite((void *) v, sizeof(float), 3, fp);
            v[0] = v1[0];
            v[1] = v1[1];
            v[2] = v1[2];
            fwrite((void *) v, sizeof(float), 3, fp);
            v[0] = v2[0];
            v[1] = v2[1];
            v[2] = v2[2];
            fwrite((void *) v, sizeof(float), 3, fp);
            fwrite(&abc, sizeof(unsigned short), 1, fp);
        }

        fclose(fp);
    }
}

template<typename T>
int SegmentSet<T>::cutSegmentWithLine(Line<T> const& line, Segment const& segment,
        SegmentSet<T>& newSegmentSet) const
{
    T epsilon = getEpsilon<T>(precision);

    int vertexTags[3];

    // Tag the segment vertices.
    for (int iVertex = 0; iVertex < 3; iVertex++) {
        Array<T,2> tmp = segment[iVertex] - line.point;
        T norm_tmp = norm(tmp);
        if (norm_tmp > epsilon) {
            tmp /= norm_tmp;
        } else {
            tmp[0] = tmp[1] = tmp[2] = (T) 0.0;
        }
        T dotp = dot(tmp, line.normal);
        if (fabs(dotp) <= epsilon) {
            vertexTags[iVertex] = 0;
        } else if (dotp > (T) 0.0 && fabs(dotp) > epsilon) {
            vertexTags[iVertex] = -1;
        } else if (dotp < (T) 0.0 && fabs(dotp) > epsilon) {
            vertexTags[iVertex] = 1;
        } else {
            return -1;
        }
    }

    // All three vertices belong to one side of the cut line.
    if (vertexTags[0] == 1 && vertexTags[1] == 1 && vertexTags[2] == 1) {
        newSegmentSet.segments.push_back(segment);
        return 1;
    } else if (vertexTags[0] == -1 && vertexTags[1] == -1 && vertexTags[2] == -1) {
        return 0;
    }

    // One vertex belongs to one side of the cut line and the other two vertices
    //   belong to the other side.
    for (int i = 0; i < 3; i++) {
        int j = (i + 1) != 3 ? (i + 1) : 0;
        int k = (j + 1) != 3 ? (j + 1) : 0;

        if (vertexTags[i] == 1 && vertexTags[j] == -1 && vertexTags[k] == -1) {
            Array<T,2> intersection_ij((T) 0.0, (T) 0.0, (T) 0.0), intersection_ik((T) 0.0, (T) 0.0, (T) 0.0);
            int rv = 0;
            rv = lineIntersectionWithLine<T>(line, segment[i], segment[j], precision, intersection_ij);
            if (rv != 1) {
                return -1;
            }
            rv = lineIntersectionWithLine<T>(line, segment[i], segment[k], precision, intersection_ik);
            if (rv != 1) {
                return -1;
            }
            Segment newSegment(segment[i], intersection_ij, intersection_ik);
            newSegmentSet.segments.push_back(newSegment);
            return 1;
        } else if (vertexTags[i] == -1 && vertexTags[j] == 1 && vertexTags[k] == 1) {
            Array<T,2> intersection_ij((T) 0.0, (T) 0.0, (T) 0.0), intersection_ik((T) 0.0, (T) 0.0, (T) 0.0);
            int rv = 0;
            rv = lineIntersectionWithLine<T>(line, segment[i], segment[j], precision, intersection_ij);
            if (rv != 1) {
                return -1;
            }
            rv = lineIntersectionWithLine<T>(line, segment[i], segment[k], precision, intersection_ik);
            if (rv != 1) {
                return -1;
            }
            Segment newSegment_0(segment[k], intersection_ij, segment[j]);
            Segment newSegment_1(segment[k], intersection_ik, intersection_ij);
            newSegmentSet.segments.push_back(newSegment_0);
            newSegmentSet.segments.push_back(newSegment_1);
            return 1;
        }
    }

    // Only one vertex belongs to the cut line.
    for (int i = 0; i < 3; i++) {
        int j = (i + 1) != 3 ? (i + 1) : 0;
        int k = (j + 1) != 3 ? (j + 1) : 0;

        if (vertexTags[i] == 0) {
            if (vertexTags[j] == 1 && vertexTags[k] == 1) {
                newSegmentSet.segments.push_back(segment);
                return 1;
            } else if (vertexTags[j] == -1 && vertexTags[k] == -1) {
                return 0;
            } else if (vertexTags[j] == 1 && vertexTags[k] == -1) {
                Array<T,2> intersection((T) 0.0, (T) 0.0, (T) 0.0);
                int rv = 0;
                rv = lineIntersectionWithLine<T>(line, segment[j], segment[k], precision, intersection);
                if (rv != 1) {
                    return -1;
                }
                Segment newSegment(segment[i], segment[j], intersection);
                newSegmentSet.segments.push_back(newSegment);
                return 1;
            } else if (vertexTags[j] == -1 && vertexTags[k] == 1) {
                Array<T,2> intersection((T) 0.0, (T) 0.0, (T) 0.0);
                int rv = 0;
                rv = lineIntersectionWithLine<T>(line, segment[j], segment[k], precision, intersection);
                if (rv != 1) {
                    return -1;
                }
                Segment newSegment(segment[i], intersection, segment[k]);
                newSegmentSet.segments.push_back(newSegment);
                return 1;
            }
        }
    }

    // Only two of the three vertices belong to the cut line.
    for (int i = 0; i < 3; i++) {
        int j = (i + 1) != 3 ? (i + 1) : 0;
        int k = (j + 1) != 3 ? (j + 1) : 0;

        if (vertexTags[i] == 0 && vertexTags[j] == 0) {
            if (vertexTags[k] == 1) {
                newSegmentSet.segments.push_back(segment);
                return 1;
            } else if (vertexTags[k] == -1) {
                return 0;
            }
        }
    }

    // All 3 vertices belong to the cut line.
    if (vertexTags[0] == 0 && vertexTags[1] == 0 && vertexTags[2] == 0) {
        newSegmentSet.segments.push_back(segment);
        return 1;
    }

    return -1;
}

template<typename T>
int SegmentSet<T>::cutWithLine(Line<T> const& line, SegmentSet<T>& newSegmentSet) const
{
    T epsilon = getEpsilon<T>(precision);

    T norm_normal = norm(line.normal);
    PLB_ASSERT(norm_normal > epsilon); // The cut line normal vector cannot have zero magnitude.
    Line<T> newLine;
    newLine.point = line.point;
    newLine.normal = line.normal / norm_normal;

    newSegmentSet.segments.resize(0);

    newSegmentSet.precision = precision;

    for (pluint iSegment = 0; iSegment < segments.size(); iSegment++) {
        if (cutSegmentWithLine(newLine, segments[iSegment], newSegmentSet) == -1) {
            return -1;
        }
    }

    if (newSegmentSet.segments.size() != 0) {
        newSegmentSet.computeMinMaxEdges();
        newSegmentSet.computeBoundingCuboid();
    }

    if (newSegmentSet.segments.size() == 0 || newSegmentSet.segments.size() == segments.size()) {
        return 0;
    }

    return 1;
}

template<typename T>
int SegmentSet<T>::cutWithLine (
        Line<T> const& line, Cuboid2D<T> const& cuboid, SegmentSet<T>& newSegmentSet ) const
{
    T epsilon = getEpsilon<T>(precision);

    T norm_normal = norm(line.normal);
    PLB_ASSERT(norm_normal > epsilon); // The cut line normal vector cannot have zero magnitude.
    Line<T> newLine;
    newLine.point = line.point;
    newLine.normal = line.normal / norm_normal;

    T norm_diagonal = norm(cuboid.upperRightCorner - cuboid.lowerLeftCorner);
    PLB_ASSERT(norm_diagonal > epsilon); // The diagonal of the cuboid cannot have zero length.

    newSegmentSet.segments.resize(0);

    newSegmentSet.precision = precision;

    for (pluint iSegment = 0; iSegment < segments.size(); iSegment++) {
        Segment const& segment = segments[iSegment];

        Array<T,2> vertices[3];
        vertices[0] = segment[0];
        vertices[1] = segment[1];
        vertices[2] = segment[2];

        // Check if the segment is fully contained in the cuboid.
        int isNotFullyContained = 0;
        for (int iVertex = 0; iVertex < 3; iVertex++) {
            Array<T,2> diff_l;
            diff_l = vertices[iVertex] - cuboid.lowerLeftCorner;

            Array<T,2> diff_u;
            diff_u = vertices[iVertex] - cuboid.upperRightCorner;

            if ((diff_l[0] < (T) 0.0 && fabs(diff_l[0]) > epsilon) ||
                (diff_l[1] < (T) 0.0 && fabs(diff_l[1]) > epsilon) ||
                (diff_l[2] < (T) 0.0 && fabs(diff_l[2]) > epsilon) ||
                (diff_u[0] > (T) 0.0 && fabs(diff_u[0]) > epsilon) ||
                (diff_u[1] > (T) 0.0 && fabs(diff_u[1]) > epsilon) ||
                (diff_u[2] > (T) 0.0 && fabs(diff_u[2]) > epsilon)) {
                isNotFullyContained = 1;
                break;
            }
        }

        if (isNotFullyContained) {
            newSegmentSet.segments.push_back(segment);
            continue;
        }

        if (cutSegmentWithLine(newLine, segment, newSegmentSet) == -1)
            return -1;
    }

    if (newSegmentSet.segments.size() != 0) {
        newSegmentSet.computeMinMaxEdges();
        newSegmentSet.computeBoundingCuboid();
    }

    if (newSegmentSet.segments.size() == 0 || newSegmentSet.segments.size() == segments.size()) {
        return 0;
    }

    return 1;
}

template<typename T>
void SegmentSet<T>::computeMinMaxEdges() {
    T nextMin, nextMax;
    for (pluint i=0; i<segments.size(); ++i) {
        computeMinMaxEdge(i, nextMin, nextMax);
        minEdgeLength = std::min(minEdgeLength, nextMin);
        maxEdgeLength = std::max(maxEdgeLength, nextMax);
    }
}

template<typename T>
void SegmentSet<T>::computeMinMaxEdge(pluint iSegment, T& minEdge, T& maxEdge) const {
    PLB_ASSERT( iSegment<segments.size() );
    Segment const& segment = segments[iSegment];
    T edge1 = norm(segment[1]-segment[0]);
    T edge2 = norm(segment[2]-segment[1]);
    T edge3 = norm(segment[0]-segment[2]);
    minEdge = std::min(edge1, std::min(edge2, edge3));
    maxEdge = std::max(edge1, std::max(edge2, edge3));
}

template<typename T>
void SegmentSet<T>::computeBoundingCuboid() {
    T xMin, yMin;
    T xMax, yMax;

    xMin = yMin  =  std::numeric_limits<T>::max();
    xMax = yMax  = -std::numeric_limits<T>::max();
    for (pluint i=0; i<segments.size(); ++i) {
        Segment const& segment = segments[i];

        xMin = std::min(xMin, std::min(segment[0][0], segment[1][0]));
        yMin = std::min(yMin, std::min(segment[0][1], segment[1][1]));

        xMax = std::max(xMax, std::max(segment[0][0], segment[1][0]));
        yMax = std::max(yMax, std::max(segment[0][1], segment[1][1]));
    }
    boundingCuboid.lowerLeftCorner  = Array<T,2>(xMin, yMin);
    boundingCuboid.upperRightCorner = Array<T,2>(xMax, yMax);
}

} // namespace plb

#endif  // SEGMENT_SET_HH
