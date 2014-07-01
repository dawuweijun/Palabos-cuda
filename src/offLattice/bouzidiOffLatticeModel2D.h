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

#ifndef BOUZIDI_OFF_LATTICE_MODEL_2D_H
#define BOUZIDI_OFF_LATTICE_MODEL_2D_H

#include "core/globalDefs.h"
#include "offLattice/offLatticeModel2D.h"

namespace plb {

template<typename T, template<typename U> class Descriptor>
class BouzidiOffLatticeModel2D : public OffLatticeModel2D<T,Array<T,2> >
{
public:
    BouzidiOffLatticeModel2D(BoundaryShape2D<T,Array<T,2> >* shape_, int flowType_);
    virtual BouzidiOffLatticeModel2D<T,Descriptor>* clone() const;
    virtual plint getNumNeighbors() const;
    virtual void prepareCell (
            Dot2D const& cellLocation, AtomicContainerBlock2D& container );
    virtual void boundaryCompletion (
            AtomicBlock2D& lattice, AtomicContainerBlock2D& container,
            std::vector<AtomicBlock2D const*> const& args );
    void cellCompletion (
            BlockLattice2D<T,Descriptor>& lattice,
            Dot2D const& boundaryNode,
            std::vector<int> const& solidDirections, std::vector<plint> const& boundaryIds,
            std::vector<bool> const& hasFluidNeighbor, Dot2D const& absoluteOffset,
            Array<T,2>& localForce, std::vector<AtomicBlock2D const*> const& args );
    virtual ContainerBlockData* generateOffLatticeInfo() const;
    virtual Array<T,2> getLocalForce(AtomicContainerBlock2D& container) const;
    void selectComputeStat(bool flag) { computeStat = flag; }
    bool computesStat() const { return computeStat; }
private:
    bool computeStat;
    std::vector<T> invAB;
private:
    /// Store the location of wall nodes, as well as the pattern of missing vs. known
    ///   populations.
    class BouzidiOffLatticeInfo2D : public ContainerBlockData {
    public:
        std::vector<Dot2D> const&               getBoundaryNodes() const
        { return boundaryNodes; }
        std::vector<Dot2D>&                     getBoundaryNodes()
        { return boundaryNodes; }
        std::vector<std::vector<int> > const&   getSolidDirections() const
        { return solidDirections; }
        std::vector<std::vector<int> >&         getSolidDirections()
        { return solidDirections; }
        std::vector<std::vector<plint> > const& getBoundaryIds() const
        { return boundaryIds; }
        std::vector<std::vector<plint> >&       getBoundaryIds()
        { return boundaryIds; }
        std::vector<std::vector<bool> > const&  getHasFluidNeighbor() const
        { return hasFluidNeighbor; }
        std::vector<std::vector<bool> >&        getHasFluidNeighbor()
        { return hasFluidNeighbor; }
        Array<T,2> const&                       getLocalForce() const
        { return localForce; }
        Array<T,2>&                             getLocalForce()
        { return localForce; }
        virtual BouzidiOffLatticeInfo2D* clone() const {
            return new BouzidiOffLatticeInfo2D(*this);
        }
    private:
        std::vector<Dot2D>               boundaryNodes;
        std::vector<std::vector<int> >   solidDirections;
        std::vector<std::vector<plint> > boundaryIds;
        std::vector<std::vector<bool> >  hasFluidNeighbor;
        Array<T,2>                       localForce;
    };
};

}  // namespace plb

#endif  // BOUZIDI_OFF_LATTICE_MODEL_2D_H

