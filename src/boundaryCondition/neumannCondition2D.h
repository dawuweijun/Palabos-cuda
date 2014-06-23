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

/** \file
 * Neumann boundary conditions -- header file.
 */
#ifndef NEUMANN_CONDITION_2D_H
#define NEUMANN_CONDITION_2D_H

#include "core/globalDefs.h"
#include "core/blockLatticeBase2D.h"
#include "atomicBlock/blockLattice2D.h"
#include "atomicBlock/dataProcessorWrapper2D.h"

namespace plb {

template<typename T, template<typename U> class Descriptor, int direction, int orientation>
class CopyUnknownPopulationsFunctional2D : public BoxProcessingFunctional2D_L<T,Descriptor>
{
public:
    virtual void process(Box2D domain, BlockLattice2D<T,Descriptor>& lattice);
    virtual CopyUnknownPopulationsFunctional2D<T,Descriptor,direction,orientation>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, template<typename U> class Descriptor, int normalX, int normalY>
class CopyAllPopulationsFunctional2D : public BoxProcessingFunctional2D_L<T,Descriptor>
{
public:
    virtual void process(Box2D domain, BlockLattice2D<T,Descriptor>& lattice);
    virtual CopyAllPopulationsFunctional2D<T,Descriptor,normalX,normalY>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, template<typename U> class Descriptor, int normalX, int normalY>
class CopyVelocityFunctional2D : public BoxProcessingFunctional2D_L<T,Descriptor>
{
public:
    virtual void process(Box2D domain, BlockLattice2D<T,Descriptor>& lattice);
    virtual CopyVelocityFunctional2D<T,Descriptor,normalX,normalY>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, template<typename U> class Descriptor, int normalX, int normalY>
class CopyTangentialVelocityFunctional2D : public BoxProcessingFunctional2D_L<T,Descriptor>
{
public:
    virtual void process(Box2D domain, BlockLattice2D<T,Descriptor>& lattice);
    virtual CopyTangentialVelocityFunctional2D<T,Descriptor,normalX,normalY>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, template<typename U> class Descriptor, int normalX, int normalY>
class CopyNormalVelocityFunctional2D : public BoxProcessingFunctional2D_L<T,Descriptor>
{
public:
    virtual void process(Box2D domain, BlockLattice2D<T,Descriptor>& lattice);
    virtual CopyNormalVelocityFunctional2D<T,Descriptor,normalX,normalY>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, template<typename U> class Descriptor, int normalX, int normalY>
class CopyDensityFunctional2D : public BoxProcessingFunctional2D_L<T,Descriptor>
{
public:
    virtual void process(Box2D domain, BlockLattice2D<T,Descriptor>& lattice);
    virtual CopyDensityFunctional2D<T,Descriptor,normalX,normalY>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

/// Outflow boundary condition that works with external rhoBar and J from the previous time step.
///     Caution is necessary with the usage of this data processor. The rhoBar and J provided,
///     must contain the values of the fields at the previous time step. This means that if
///     the "external-rhoBar-J-collide-and-stream" is used at a processor level 0, then
///     this data processor must be integrated at processor level 1, and the data processor
///     to recompute the new rhoBar and J should be integrated at processor level 2.
template<typename T, template<typename U> class Descriptor>
class VirtualOutlet2D : public BoxProcessingFunctional2D
{
public:
    /* Type 0: Close to FluidPressureOutlet2D (imposes a strict pressure).
     * Type 1: Laplacian filter / extrapolation on the pressure.
     **/
    VirtualOutlet2D(T outsideDensity_, Box2D globalDomain_, int type_=1);
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D*> blocks);
    virtual VirtualOutlet2D<T,Descriptor>* clone() const
    {
        return new VirtualOutlet2D<T,Descriptor>(*this);
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const
    {
        modified[0] = modif::staticVariables; // Block lattice.
        modified[1] = modif::nothing;         // RhoBar.
        modified[2] = modif::nothing;         // J.
    }
private:
    T outsideDensity;   // Boundary condition for the density (usually 1.0).
    Box2D globalDomain; // The globalDomain must be at most as big as the whole simulation
                        // domain for non-periodic problems, and bigger than the whole simulation
                        // domain plus the envelope (per periodic direction) for periodic problems.
    int type;           // If type = 0 then this is very close to FluidPressureOutlet2D.
                        // If type = 1 some times gives the best results.
};

}  // namespace plb

#endif  // NEUMANN_CONDITION_2D_H
