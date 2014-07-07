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

#ifndef GUO_ADV_DIFF_OFF_LATTICE_MODEL_2D_H
#define GUO_ADV_DIFF_OFF_LATTICE_MODEL_2D_H

#include "core/globalDefs.h"
#include "core/array.h"
#include "offLattice/offLatticeModel2D.h"
#include "offLattice/boundaryShapes2D.h"

namespace plb {

template<typename T, template<typename U> class Descriptor>
class GuoAdvDiffOffLatticeModel2D : public OffLatticeModel2D<T,Array<T,2> >
{
public:
    GuoAdvDiffOffLatticeModel2D(BoundaryShape2D<T,Array<T,2> >* shape_, int flowType_);
    GuoAdvDiffOffLatticeModel2D(GuoAdvDiffOffLatticeModel2D<T,Descriptor> const& rhs);
    GuoAdvDiffOffLatticeModel2D<T,Descriptor>& operator=(GuoAdvDiffOffLatticeModel2D<T,Descriptor> const& rhs);
    virtual GuoAdvDiffOffLatticeModel2D<T,Descriptor>* clone() const;
    virtual plint getNumNeighbors() const;
    virtual void prepareCell (
            Dot2D const& cellLocation, AtomicContainerBlock2D& container );
    virtual void boundaryCompletion (
            AtomicBlock2D& lattice, AtomicContainerBlock2D& container,
            std::vector<AtomicBlock2D const*> const& args );
    virtual ContainerBlockData* generateOffLatticeInfo() const;
    virtual Array<T,2> getLocalForce(AtomicContainerBlock2D& container) const { return Array<T,2>(T(),T(),T()); }
    void selectSecondOrder(bool flag) { secondOrderFlag = flag; }
    bool usesSecondOrder() const { return secondOrderFlag; }
private:
    void cellCompletion (
            BlockLattice2D<T,Descriptor>& lattice,
            Dot2D const& guoNode,
            std::vector<std::pair<int,int> > const& dryNodeFluidDirections,
            std::vector<plint> const& dryNodeIds, Dot2D const& absoluteOffset );
    void computeRhoBarJNeq (
              BlockLattice2D<T,Descriptor> const& lattice, Dot2D const& guoNode,
              Dot2D const& fluidDirection, int depth, Array<T,2> const& wallNode, T delta,
              Array<T,2> wallData, OffBoundary::Type bdType, Array<T,2> const& wallNormal,
              T& rhoBar, Array<T,Descriptor<T>::d>& jNeq ) const;
private:
    bool secondOrderFlag;
private:
    /// Store the location of wall nodes, as well as the pattern of missing vs. known
    ///   populations.
    class GuoAdvDiffOffLatticeInfo2D : public ContainerBlockData {
    public:
        std::vector<Dot2D> const&                               getDryNodes() const
        { return dryNodes; }
        std::vector<Dot2D>&                                     getDryNodes()
        { return dryNodes; }
        std::vector<std::vector<std::pair<int,int> > > const&   getDryNodeFluidDirections() const
        { return dryNodeFluidDirections; }
        std::vector<std::vector<std::pair<int,int> > >&         getDryNodeFluidDirections()
        { return dryNodeFluidDirections; }
        std::vector<std::vector<plint> > const&                 getDryNodeIds() const
        { return dryNodeIds; }
        std::vector<std::vector<plint> >&                       getDryNodeIds()
        { return dryNodeIds; }
        virtual GuoAdvDiffOffLatticeInfo2D* clone() const {
            return new GuoAdvDiffOffLatticeInfo2D(*this);
        }
    private:
        std::vector<Dot2D> dryNodes;
        std::vector<std::vector<std::pair<int,int> > >   dryNodeFluidDirections;
        std::vector<std::vector<plint> >                 dryNodeIds;
    };
};

}  // namespace plb

#endif  // GUO_ADV_DIFF_OFF_LATTICE_MODEL_2D_H
