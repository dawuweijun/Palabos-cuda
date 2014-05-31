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

#ifndef BUBBLE_MATCH_2D_H
#define BUBBLE_MATCH_2D_H

#include "core/globalDefs.h"
#include "offLattice/makeSparse2D.h"
#include "parallelism/mpiManager.h"
#include "atomicBlock/reductiveDataProcessingFunctional2D.h"
#include "atomicBlock/atomicContainerBlock2D.h"
#include "atomicBlock/dataProcessingFunctional2D.h"


namespace plb
{

class BubbleMPIdata2D
{
public:
    BubbleMPIdata2D ( MultiBlock2D& block );
    std::vector<plint> const& getLocalIds() const;
private:
    void computeLocalIds ( MultiBlock2D& block );
private:
    std::vector<plint> localIds;
};

class BubbleMatch2D
{
public:
    BubbleMatch2D ( MultiBlock2D& templ, bool matchEmpty_=true );
    ~BubbleMatch2D();
    template<typename T>
    void execute ( MultiScalarField2D<int>& flag, MultiScalarField2D<T>& volumeFraction );
    MultiScalarField2D<plint>* getTagMatrix()
    {
        return tagMatrix;
    }
    void setTagMatrix ( MultiScalarField2D<plint>* newTagMatrix )
    {
        tagMatrix = newTagMatrix;
    }
    std::vector<double> const&  getBubbleVolume()
    {
        return bubbleVolume;
    }
    std::vector<Array<double,2> > const&  getBubbleCenter()
    {
        return bubbleCenter;
    }
    pluint numBubbles() const
    {
        return bubbleVolume.size();
    }
private:
    // Re-assign a continuously numbered ID to the detected bubbles.
    pluint countAndTagBubbles();
    // Computes the volumes and centers of all new bubbles.
    template<typename T>
    void bubbleAnalysis ( MultiScalarField2D<int>& flag, MultiScalarField2D<T>& volumeFraction, pluint numBubbles );
    // Implements all required MPI operations needed to compute the bubble volume and centers,
    // after calling AnalyzeBubbles2D.
    void computeBubbleData ( pluint numBubbles );
    // Implements all required MPI operations needed to compute the global IDs of the current
    // bubbbles, after calling CollectBubbleTags2D.
    plint globalBubbleIds();
    // Prepare the bubble container for the next time iteration.
    void resetBubbleContainer();
    // Parallel bucket-fill algorithm to assign a unique ID to every contiguous region.
    void bubbleBucketFill ( MultiScalarField2D<int>& flag );
private:
    BubbleMatch2D ( BubbleMatch2D const& rhs ) : mpiData ( rhs.mpiData )
    {
        PLB_ASSERT ( false );
    }
    BubbleMatch2D& operator= ( BubbleMatch2D const& rhs )
    {
        PLB_ASSERT ( false );
        return *this;
    }
private:
    MultiContainerBlock2D *bubbleContainer, *bubbleAnalysisContainer, *bubbleRemapContainer;
    BubbleMPIdata2D mpiData;
    MultiScalarField2D<plint> *tagMatrix;
    std::vector<double> bubbleVolume;
    std::vector<Array<double,2> > bubbleCenter;
    bool matchEmpty;
    static const plint maxNumBubbles = 100000;
};



/**
 * Data for the bubble counter, associated to one block.
 * It holds information that changes during time:
 *  nextCellId: next available ID, if a new cell-type is found.
 *  retagging: a map that relates equivalent IDs, when two domains
 *             are found to be contiguous.
 *  maxNumBubbles: an upper bound for the number of bubbles per block,
 *                 so every block can create a globally unique bubble ID.
 **/
class BubbleCounterData2D : public ContainerBlockData
{
public:
    BubbleCounterData2D ( plint maxNumBubbles_ );
    virtual BubbleCounterData2D* clone() const;
    // The assumption here is that cell 0 is a cell that needs to
    // be tagged ("a bubble cell"). Depending on the tag of neighboring
    // cells, either one of the neighbor's tag or a new tag is assigned.
    // There is a conflict if cell 0 has been previously converted in
    // a way that is incompatible with the neighbor's tags.
    bool convertCell (
        plint& tag0,
        plint tag1,  plint tag2,  plint tag3,
        plint tag4,  plint tag5,  plint tag6,
        plint tag7,  plint tag8 );
    // The assumption here is that cell 0 is a cell that needs to
    // be tagged ("a bubble cell"). Depending on the tag of neighboring
    // cell 1, either cell 1's tag or a new tag is assigned.
    // There is a conflict if cell 0 has been previously converted in
    // a way that is incompatible with cell 1's tag.
    bool processNeighbor ( plint& tag0, plint tag1 );
    plint getNextTag();
    void reset();
    // Get the "real" value (after remapping) of a given tag.
    plint convertTag ( plint tag ) const;
    // Important: it is assumed that oldTag and newTag are not
    // themselves mapped. This means that for arbitrary tags tag1 and
    // tag2 you should not call registerConflict(tag1,tag2), but
    // registerConflict(convertTag(tag1), convertTag(tag2).

    void registerConflict ( plint oldTag, plint newTag );
private:
    plint nextCellId;
    std::map<plint,plint> retagging;
private:
    plint maxNumBubbles;
};

class BubbleRemapData2D : public ContainerBlockData
{
public:
    BubbleRemapData2D ( plint maxNumBubbles_=0 )
        : maxNumBubbles ( maxNumBubbles_ )
    { }
    virtual BubbleRemapData2D* clone() const;
    std::vector<plint>& getUniqueTags()
    {
        return uniqueTags;
    }
    std::vector<plint> const& getUniqueTags() const
    {
        return uniqueTags;
    }
    std::map<plint,plint>& getTagRemap()
    {
        return tagRemap;
    }
    bool isMyTag ( plint tag );
private:
    plint maxNumBubbles;
    std::vector<plint> uniqueTags;
    std::map<plint,plint> tagRemap;
};

struct BubbleAnalysisData2D : public ContainerBlockData
{
    virtual BubbleAnalysisData2D* clone() const
    {
        return new BubbleAnalysisData2D ( *this );
    }
    std::vector<double> bubbleVolume;
    std::vector<Array<double,2> > bubbleCenter;
};

class CountBubbleIteration2D : public PlainReductiveBoxProcessingFunctional2D
{
public:
    CountBubbleIteration2D ( bool matchEmpty_ );
    virtual void processGenericBlocks ( Box2D domain, std::vector<AtomicBlock2D*> atomicBlocks );
    virtual CountBubbleIteration2D* clone() const;
    virtual void getTypeOfModification ( std::vector<modif::ModifT>& modified ) const
    {
        modified[0] = modif::staticVariables; // tags.
        modified[1] = modif::nothing;         // flags.
        modified[2] = modif::nothing;         // data.
    }
    plint getNumConflicts() const;
private:
    plint numConflictsId;
    bool matchEmpty;
};

template<typename T>
class AnalyzeBubbles2D : public BoxProcessingFunctional2D
{
public:
    AnalyzeBubbles2D ( pluint numBubbles_, bool matchEmpty_ );
    virtual void processGenericBlocks ( Box2D domain, std::vector<AtomicBlock2D*> atomicBlocks );
    virtual AnalyzeBubbles2D<T>* clone() const;
    virtual void getTypeOfModification ( std::vector<modif::ModifT>& modified ) const
    {
        modified[0] = modif::staticVariables; // tags.
        modified[1] = modif::nothing;         // flags.
        modified[2] = modif::nothing;         // data.
        modified[3] = modif::nothing;         // volume fraction.
    }
private:
    pluint numBubbles;
    bool matchEmpty;
};


// Converts the information about the overall available bubble tags available
// to all processors.
class CollectBubbleTags2D : public BoxProcessingFunctional2D
{
public:
    virtual void processGenericBlocks ( Box2D domain, std::vector<AtomicBlock2D*> atomicBlocks );
    virtual CollectBubbleTags2D* clone() const;
    virtual void getTypeOfModification ( std::vector<modif::ModifT>& modified ) const
    {
        modified[0] = modif::nothing;  // tags.
        modified[1] = modif::nothing;  // data.
    }
};

// Assign a new tag to all bubble cells (they must have been uniquely tagged previously).
// The only field in the BubbleCounterData2D which is used here is tagRemap.
class ApplyTagRemap2D : public BoxProcessingFunctional2D
{
public:
    virtual void processGenericBlocks ( Box2D domain, std::vector<AtomicBlock2D*> atomicBlocks );
    virtual ApplyTagRemap2D* clone() const;
    virtual void getTypeOfModification ( std::vector<modif::ModifT>& modified ) const
    {
        modified[0] = modif::staticVariables; // tags.
        modified[1] = modif::nothing;         // data.
    }
};

}  // namespace plb

#endif  // BUBBLE_MATCH_2D_H

