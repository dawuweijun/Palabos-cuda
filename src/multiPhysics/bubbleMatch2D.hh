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

#ifndef BUBBLE_MATCH_2D_HH
#define BUBBLE_MATCH_2D_HH

#include "multiPhysics/bubbleMatch2D.h"
#include "offLattice/makeSparse2D.h"
#include "parallelism/mpiManager.h"
#include "atomicBlock/reductiveDataProcessingFunctional2D.h"
#include "multiBlock/reductiveMultiDataProcessorWrapper2D.h"
#include "atomicBlock/atomicContainerBlock2D.h"
#include "atomicBlock/dataProcessingFunctional2D.h"
#include "multiPhysics/freeSurfaceUtil2D.h"


namespace plb
{


template<typename T>
void BubbleMatch2D::execute ( MultiScalarField2D<int>& flag, MultiScalarField2D<T>& volumeFraction )
{
    bubbleBucketFill ( flag );
    pluint numBubbles = countAndTagBubbles();
    bubbleVolume.clear();
    bubbleCenter.clear();
    bubbleAnalysis ( flag, volumeFraction, numBubbles );
}

template<typename T>
void BubbleMatch2D::bubbleAnalysis ( MultiScalarField2D<int>& flag,
                                     MultiScalarField2D<T>& volumeFraction, pluint numBubbles )
{
    std::vector<MultiBlock2D*> args;
    args.push_back ( tagMatrix );
    args.push_back ( &flag );
    args.push_back ( bubbleAnalysisContainer );
    args.push_back ( &volumeFraction );
    applyProcessingFunctional ( new AnalyzeBubbles2D<T> ( numBubbles, matchEmpty ),
                                bubbleAnalysisContainer->getBoundingBox(), args );
    computeBubbleData ( numBubbles );
}

/* *************** Class AnalyzeBubbles2D ******************************** */

template<typename T>
AnalyzeBubbles2D<T>::AnalyzeBubbles2D ( pluint numBubbles_, bool matchEmpty_ )
    : numBubbles ( numBubbles_ ),
      matchEmpty ( matchEmpty_ )
{ }

template<typename T>
AnalyzeBubbles2D<T>* AnalyzeBubbles2D<T>::clone() const
{
    return new AnalyzeBubbles2D<T> ( *this );
}

template<typename T>
void AnalyzeBubbles2D<T>::processGenericBlocks ( Box2D domain,std::vector<AtomicBlock2D*> atomicBlocks )
{
    PLB_ASSERT ( atomicBlocks.size() ==4 );
    ScalarField2D<plint>* pTagMatrix = dynamic_cast<ScalarField2D<plint>*> ( atomicBlocks[0] );
    PLB_ASSERT ( pTagMatrix );
    ScalarField2D<plint>& tagMatrix = *pTagMatrix;

    ScalarField2D<int>* pFlagMatrix = dynamic_cast<ScalarField2D<int>*> ( atomicBlocks[1] );
    PLB_ASSERT ( pFlagMatrix );
    ScalarField2D<int>& flagMatrix = *pFlagMatrix;

    AtomicContainerBlock2D* pDataBlock = dynamic_cast<AtomicContainerBlock2D*> ( atomicBlocks[2] );
    PLB_ASSERT ( pDataBlock );
    AtomicContainerBlock2D& dataBlock = *pDataBlock;
    BubbleAnalysisData2D* pData = dynamic_cast<BubbleAnalysisData2D*> ( dataBlock.getData() );
    PLB_ASSERT ( pData );
    BubbleAnalysisData2D& data = *pData;

    ScalarField2D<T>* pVolumeFraction = dynamic_cast<ScalarField2D<T>*> ( atomicBlocks[3] );
    PLB_ASSERT ( pVolumeFraction );
    ScalarField2D<T>& volumeFraction = *pVolumeFraction;

    Dot2D flagOffset = computeRelativeDisplacement ( tagMatrix, flagMatrix );
    Dot2D vfOffset = computeRelativeDisplacement ( tagMatrix, volumeFraction );
    Dot2D absOfs = tagMatrix.getLocation();

    std::vector<double> bubbleVolume ( numBubbles );
    std::fill ( bubbleVolume.begin(), bubbleVolume.end(), 0. );
    std::vector<Array<double,3> > bubbleCenter ( numBubbles );
    std::fill ( bubbleCenter.begin(), bubbleCenter.end(), Array<double,3> ( 0.,0.,0. ) );

    for ( plint iX=domain.x0; iX<=domain.x1; ++iX )
    {
        for ( plint iY=domain.y0; iY<=domain.y1; ++iY )
        {
            plint tag = tagMatrix.get ( iX,iY );
            if ( tag>=0 )
            {
                if ( ( matchEmpty && flagMatrix.get ( iX+flagOffset.x,iY+flagOffset.y ) ==twoPhaseFlag::empty ) ||
                        ( !matchEmpty && flagMatrix.get ( iX+flagOffset.x,iY+flagOffset.y ) ==twoPhaseFlag::fluid ) )
                {
                    PLB_ASSERT ( tag < ( plint ) bubbleVolume.size() );
                    bubbleVolume[tag] += 1.0;
                    bubbleCenter[tag] += Array<double,3> ( ( double ) iX+absOfs.x, ( double ) iY+absOfs.y );
                }
                else if ( flagMatrix.get ( iX+flagOffset.x,iY+flagOffset.y ) ==twoPhaseFlag::interface )
                {
                    PLB_ASSERT ( tag < ( plint ) bubbleVolume.size() );
                    double vf = ( double ) volumeFraction.get ( iX+vfOffset.x,iY+vfOffset.y );
                    if ( matchEmpty )
                    {
                        vf = 1.0 - vf;
                    }
                    bubbleVolume[tag] += vf;
                    bubbleCenter[tag] += vf*Array<double,3> ( ( double ) iX+absOfs.x, ( double ) iY+absOfs.y );
                }
                else
                {
                    PLB_ASSERT ( false );
                }
            }
        }
    }

    data.bubbleVolume = bubbleVolume;
    data.bubbleCenter = bubbleCenter;
}

}  // namespace plb

#endif  // BUBBLE_MATCH_2D_HH

