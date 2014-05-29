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

#ifndef FREE_SURFACE_BOUNDARY_CONDITION_2D_H
#define FREE_SURFACE_BOUNDARY_CONDITION_2D_H

namespace plb {

template<typename T, template<typename U> class Descriptor>
class FreeSurfaceFadingArea2D : public BoxProcessingFunctional2D_L<T,Descriptor>
{
public:
    FreeSurfaceFadingArea2D(T factor_);
    virtual void process(Box2D domain, BlockLattice2D<T,Descriptor>& lattice);
    virtual FreeSurfaceFadingArea2D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::staticVariables;  // Fluid
    }
private:
    T factor;
};

template< typename T,template<typename U> class Descriptor>
class PouringLiquid2D : public BoxProcessingFunctional2D {
public:
    PouringLiquid2D(Dynamics<T,Descriptor>* dynamicsTemplate_, Array<T,3> injectionVelocity_)
        : dynamicsTemplate(dynamicsTemplate_), injectionVelocity(injectionVelocity_)
    { }
    PouringLiquid2D(PouringLiquid2D<T,Descriptor> const& rhs)
        : dynamicsTemplate(rhs.dynamicsTemplate->clone()),
          injectionVelocity(rhs.injectionVelocity)
    { }
    PouringLiquid2D<T,Descriptor>* operator=(PouringLiquid2D<T,Descriptor> const& rhs)
    { 
        PouringLiquid2D<T,Descriptor>(rhs).swap(*this);
        return *this;
    }
    void swap(PouringLiquid2D<T,Descriptor>& rhs) {
        std::swap(dynamicsTemplate, rhs.dynamicsTemplate);
        std::swap(injectionVelocity, rhs.injectionVelocity);
    }
    virtual ~PouringLiquid2D() {
        delete dynamicsTemplate;
    }
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D*> atomicBlocks);
    virtual PouringLiquid2D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::dataStructure;    // Fluid
        modified[1] = modif::staticVariables;  // rhoBar.
        modified[2] = modif::staticVariables;  // j.
        modified[3] = modif::staticVariables;  // Mass
        modified[4] = modif::staticVariables;  // Mass-fraction
        modified[5] = modif::staticVariables;  // Flag-status
        modified[6] = modif::nothing;          // Normal
        modified[7] = modif::nothing;          // Interface lists
        modified[8] = modif::nothing;          // Curvature.
        modified[9] = modif::nothing;          // Outside density.
    }
private:
    Dynamics<T,Descriptor>* dynamicsTemplate;
    Array<T,3> injectionVelocity;
};

template< typename T,template<typename U> class Descriptor>
class RemoveMass2D : public BoxProcessingFunctional2D {
public:
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D*> atomicBlocks);
    virtual RemoveMass2D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::dataStructure;    // Fluid
        modified[1] = modif::staticVariables;  // rhoBar.
        modified[2] = modif::staticVariables;  // j.
        modified[3] = modif::staticVariables;  // Mass
        modified[4] = modif::staticVariables;  // Mass-fraction
        modified[5] = modif::staticVariables;  // Flag-status
        modified[6] = modif::nothing;          // Normal
        modified[7] = modif::nothing;          // Interface lists
        modified[8] = modif::nothing;          // Curvature.
        modified[9] = modif::nothing;          // Outside density.
    }
};

template< typename T,template<typename U> class Descriptor>
class ShortenBounceBack2D : public BoxProcessingFunctional2D {
public:
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D*> atomicBlocks);
    virtual ShortenBounceBack2D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::dataStructure;    // Fluid
        modified[1] = modif::staticVariables;  // rhoBar.
        modified[2] = modif::staticVariables;  // j.
        modified[3] = modif::staticVariables;  // Mass
        modified[4] = modif::staticVariables;  // Mass-fraction
        modified[5] = modif::staticVariables;  // Flag-status
        modified[6] = modif::nothing;          // Normal
        modified[7] = modif::nothing;          // Interface lists
        modified[8] = modif::nothing;          // Curvature.
        modified[9] = modif::nothing;          // Outside density.
    }
};

}  // namespace plb

#endif  // FREE_SURFACE_BOUNDARY_CONDITION_2D_H

