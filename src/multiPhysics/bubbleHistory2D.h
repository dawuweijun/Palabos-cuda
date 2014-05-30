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

#ifndef BUBBLE_HISTORY_2D_H
#define BUBBLE_HISTORY_2D_H

#include "core/globalDefs.h"
#include "multiPhysics/bubbleMatch2D.h"
#include "offLattice/makeSparse2D.h"
#include "parallelism/mpiManager.h"
#include "atomicBlock/reductiveDataProcessingFunctional2D.h"
#include "atomicBlock/atomicContainerBlock2D.h"
#include "atomicBlock/dataProcessingFunctional2D.h"
#include "multiPhysics/freeSurfaceUtil2D.h"
#include "multiPhysics/freeSurfaceModel2D.h"
#include <sstream>
#include <string>

namespace plb {

struct BubbleTransition2D;
class BubbleInfo2D;
struct FullBubbleRecord;

template<typename T>
class BubbleHistory2D {
public:
    BubbleHistory2D(MultiBlock2D& templ);
    ~BubbleHistory2D();
    void transition(BubbleMatch2D& bubbleMatch, plint iterationStep, T newBubbleVolumeCorrection);
    // Based on the bubble information (original volume, current volume), update the pressure field
    // for all bubbles.
    void updateBubblePressure(MultiScalarField2D<T>& outsideDensity, T rhoEmpty);
    void freeze();
    void timeHistoryLog(std::string fName);
    void fullBubbleLog(std::string fName);
    std::map<plint,BubbleInfo2D> const& getBubbles() const { return bubbles; }
    std::map<plint,BubbleInfo2D>& getBubbles() { return bubbles; }
private:
    // Implement the time evolution of the bubbles:
    // - Reassign an ID to the new bubbles, compatible with the old ID.
    // - Create a data structure for the new bubbles.
    // - Print a message if bubbles are created/deleted.
    void matchAndRemapBubbles( BubbleMatch2D& bubbleMatch, MultiScalarField2D<plint>& tagMatrix1,
                               MultiScalarField2D<plint>& tagMatrix2, T newBubbleVolumeCorrection, plint iterationStep);
    // Create the map "newToAllOldMap" which maps the temporary new IDs (continuously
    // numbered) to the "old" IDs (either existing ones or newly attributed).
    void correlateBubbleIds (
            MultiScalarField2D<plint>& tagMatrix1, MultiScalarField2D<plint>& tagMatrix2,
            std::vector<std::vector<plint> >& newToAllOldMap, pluint numBubbles );
    // Updates the (non-parallel) data structure which holds the overview of the currently
    // available bubbles.
    void updateBubbleInformation (
            std::vector<BubbleTransition2D>& bubbleTransitions, std::vector<double> const& bubbleVolume,
            std::vector<Array<double,3> > const& bubbleCenter, T newBubbleVolumeCorrection, plint iterationStep );
    void computeNewBubbles (
            std::set<plint>& oldIDs, std::set<plint>& newIDs,
            std::vector<double> const& bubbleVolume, std::vector<Array<double,3> > const& bubbleCenter,
            T newBubbleVolumeCorrection, std::map<plint,BubbleInfo2D>& newBubbles, std::map<plint,plint>& newToFinal );
    void updateBubbleLog (
            BubbleTransition2D& bubbleTransition, std::vector<double> const& bubbleVolume, plint iterationStep,
            std::map<plint,BubbleInfo2D>& newBubbles, std::map<plint,plint>& newToFinal );
    static void computeBubbleTransitions (
            std::vector<std::vector<plint> > const& newToAllOldMap, std::map<plint,std::vector<plint> > const& oldToAllNewMap,
            std::vector<BubbleTransition2D>& bubbleTransitions );
private:
    BubbleHistory2D(BubbleHistory2D<T> const& rhs) { PLB_ASSERT( false ); }
    BubbleHistory2D<T>& operator=(BubbleHistory2D<T> const& rhs) { PLB_ASSERT( false ); return *this; }
private:
    MultiContainerBlock2D *bubbleAnalysisContainer, *bubbleCorrelationContainer, *bubbleRemapContainer;
    BubbleMPIdata2D mpiData;
    MultiScalarField2D<plint> *oldTagMatrix;
    std::map<plint,BubbleInfo2D> bubbles;
    plint nextBubbleID;
    //    <Iteration>   <IDs of created bubbles>,<IDs of vanished bubbles>
    std::map<plint, std::pair<std::vector<plint>,std::vector<plint> > > timeHistory;
    std::vector<FullBubbleRecord> fullBubbleRecord;
};

struct BubbleCorrelationData2D : public ContainerBlockData {
    virtual BubbleCorrelationData2D* clone() const {
        return new BubbleCorrelationData2D(*this);
    }
    // Maps new bubble tags to old bubble tags at successive iterations.
    // This is needed because a new tag may map to several old tags, due
    // to a bubble merging process.
    std::vector<plint> newToOldMap0, newToOldMap1;
};

class BubbleInfo2D {
public:
    BubbleInfo2D()
        : referenceVolume(0.),
          currentVolume(0.),
          frozen(false)
    { }
    BubbleInfo2D(double volume, Array<double,3> const& center_)
        : referenceVolume(volume),
          currentVolume(volume),
          center(center_),
          frozen(false)
    { }
    void freeze() { frozen=true; }
    void setVolume(double newVolume) {
        if (!frozen) {
            currentVolume = newVolume;
        }
    }
    void setReferenceVolume(double newReferenceVolume) {
        if (!frozen) {
            referenceVolume = newReferenceVolume;
        }
    }
    double getVolumeRatio() const {
        static const double epsilon = std::numeric_limits<double>::epsilon()*1.e4;
        if (fabs(currentVolume)>epsilon) {
            return referenceVolume/currentVolume;
        }
        else {
            return 1.0;
        }
    }
    double getReferenceVolume() const { return referenceVolume; }
    double getVolume() const { return currentVolume; }
    Array<double,3> const& getCenter() const { return center; }
    bool isFrozen() const { return frozen; }
private:
    double referenceVolume;
    double currentVolume;
    Array<double,3> center;
    bool frozen;
};

struct BubbleTransition2D {
    bool empty() const {
        return oldIDs.empty() && newIDs.empty();
    }
    std::string description() const {
        std::stringstream sstream;
        if (oldIDs.empty() && newIDs.empty()) { }
        else if (oldIDs.empty()) {
            PLB_ASSERT(newIDs.size()==1);
            sstream << "Bubble " << *newIDs.begin() << " created";
        }
        else if (newIDs.empty()) {
            PLB_ASSERT(oldIDs.size()==1);
            sstream << "Bubble " << *oldIDs.begin() << " vanished";
        }
        else if (oldIDs.size()==1) {
            if (newIDs.size()==1) {
                sstream << "Sraight transition from ID "
                        << *oldIDs.begin() << " to ID "
                        << *newIDs.begin();
            }
            else {
                sstream << "Splitting from ID "
                        << *oldIDs.begin() << " to IDs";
                std::set<plint,plint>::const_iterator it = newIDs.begin();
                for (; it!=newIDs.end(); ++it) {
                    sstream << " " << *it;
                }
            }
        }
        else if (newIDs.size()==1) {
            sstream << "Merging IDs";
            std::set<plint,plint>::const_iterator it = oldIDs.begin();
            for(; it!=oldIDs.end(); ++it) {
                sstream << " " << *it;
            }
            sstream << " into ID " << *newIDs.begin();
        }
        else {
            sstream << "Transition from IDs";
            std::set<plint,plint>::const_iterator it1 = oldIDs.begin();
            for(; it1!=oldIDs.end(); ++it1) {
                sstream << " " << *it1;
            }
            sstream << " into IDs ";
            std::set<plint,plint>::const_iterator it2 = newIDs.begin();
            for (; it2!=newIDs.end(); ++it2) {
                sstream << " " << *it2;
            }
        }
        return sstream.str();
    }
    std::set<plint> oldIDs, newIDs;
};

struct FullBubbleRecord {
    FullBubbleRecord(double initialVolume_, plint beginIteration_)
        : initialVolume(initialVolume_),
          finalVolume(0.),
          beginIteration(beginIteration_),
          endIteration(beginIteration_),
          frozen(false)
    { }

    std::string description(plint ID) const {
        std::stringstream sstream;
        sstream << "Bubble " << ID
                << ". At it. " << beginIteration
                << ", Vol. " << initialVolume
                << ": "
                << beginTransition.description() << ".";
        if (endIteration>beginIteration) {
            sstream << " Removed at it. " << endIteration
                << ", Vol. " << finalVolume
                << ": "
                << endTransition.description() << ".";
        }
        sstream << " Frozen bubble: " << ( frozen ? "Yes." : "No." );
        return sstream.str();
    }

    double initialVolume, finalVolume;
    plint beginIteration, endIteration;
    BubbleTransition2D beginTransition, endTransition;
    bool frozen;
};


template<typename T>
class CorrelateBubbleIds2D : public BoxProcessingFunctional2D
{
public:
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D*> atomicBlocks);
    virtual CorrelateBubbleIds2D<T>* clone() const;
    virtual void getTypeOfModification (std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::nothing;         // tags before.
        modified[1] = modif::nothing;         // tags after.
        modified[2] = modif::nothing;         // data.
    }
};


template<typename T>
class UpdateBubblePressure2D : public BoxProcessingFunctional2D_SS<plint,T>
{
public:
    UpdateBubblePressure2D(std::map<plint,BubbleInfo2D> const& bubbles_, T rho0_);
    virtual void process(Box2D domain, ScalarField2D<plint>& tags,
                                       ScalarField2D<T>& density);
    virtual UpdateBubblePressure2D* clone() const;
    virtual void getTypeOfModification (std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::nothing;          // tags.
        modified[1] = modif::staticVariables;  // density.
    }
private:
    std::map<plint,BubbleInfo2D> bubbles;
    T rho0;
};

}  // namespace plb

#endif  // BUBBLE_HISTORY_2D_H

