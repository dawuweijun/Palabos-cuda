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

#ifndef IMMERSED_WALLS_2D_HH
#define IMMERSED_WALLS_2D_HH

#include "core/globalDefs.h"
#include "core/array.h"
#include "atomicBlock/dataProcessingFunctional2D.h"
#include "atomicBlock/dataField2D.h"

#include "immersedWalls2D.h"

namespace plb {

/* ******** ReduceAxialTorqueImmersed2D ************************************ */

template<typename T>
ReduceAxialTorqueImmersed2D<T>::ReduceAxialTorqueImmersed2D (
        Array<T,2> const& center_, Array<T,2> const& unitaryAxis_, int reductionFlag_ )
    : center(center_),
      unitaryAxis(unitaryAxis_),
      sum_torque_ids (
            Array<plint,2> (
                this->getStatistics().subscribeSum(),
                this->getStatistics().subscribeSum(),
                this->getStatistics().subscribeSum() ) ),
      reductionFlag(reductionFlag_)
{ }

template<typename T>
void ReduceAxialTorqueImmersed2D<T>::processGenericBlocks (
        Box2D domain, std::vector<AtomicBlock2D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==1 );
    AtomicContainerBlock2D* container = dynamic_cast<AtomicContainerBlock2D*>(blocks[0]);
    PLB_ASSERT( container );

    ImmersedWallData2D<T>* wallData = 
        dynamic_cast<ImmersedWallData2D<T>*>( container->getData() );
    PLB_ASSERT(wallData);
    std::vector< Array<T,2> > const& vertices = wallData->vertices;
    std::vector< Array<T,2> > const& g = wallData->g;
    std::vector<int> const& flags = wallData->flags;
    Array<T,2> offset = wallData->offset;
    PLB_ASSERT( vertices.size()==g.size() );
    PLB_ASSERT( vertices.size()==flags.size() );

    for (pluint i=0; i<vertices.size(); ++i) {
        Array<T,2> vertex = vertices[i];
        if ( flags[i]==reductionFlag &&
             closedOpenContained(vertex, domain) )
        {
            Array<T,2> physVertex = vertex+offset;
            Array<T,2> r(physVertex-center);
            r -= dot(r,unitaryAxis)*unitaryAxis;
            Array<T,2> torque(crossProduct(r,g[i]));
            this->getStatistics().gatherSum(sum_torque_ids[0], torque[0]);
            this->getStatistics().gatherSum(sum_torque_ids[1], torque[1]);
        }
    }
}

template<typename T>
ReduceAxialTorqueImmersed2D<T>* ReduceAxialTorqueImmersed2D<T>::clone() const {
    return new ReduceAxialTorqueImmersed2D<T>(*this);
}

template<typename T>
void ReduceAxialTorqueImmersed2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing; // Container Block.
}

template<typename T>
BlockDomain::DomainT ReduceAxialTorqueImmersed2D<T>::appliesTo() const {
    return BlockDomain::bulk;
}

template<typename T>
Array<T,2> ReduceAxialTorqueImmersed2D<T>::getSumTorque() const {
    return Array<T,2> (
            this->getStatistics().getSum(sum_torque_ids[0]),
            this->getStatistics().getSum(sum_torque_ids[1]) );
}

/* ******** ReduceImmersedForce2D ************************************ */

template<typename T>
ReduceImmersedForce2D<T>::ReduceImmersedForce2D(int reductionFlag_)
    : sum_g_ids(
            Array<plint,2> (
                this->getStatistics().subscribeSum(),
                this->getStatistics().subscribeSum(),
                this->getStatistics().subscribeSum() ) ),
      reductionFlag(reductionFlag_)
{ }

template<typename T>
void ReduceImmersedForce2D<T>::processGenericBlocks (
        Box2D domain, std::vector<AtomicBlock2D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==1 );
    AtomicContainerBlock2D* container = dynamic_cast<AtomicContainerBlock2D*>(blocks[0]);
    PLB_ASSERT( container );

    ImmersedWallData2D<T>* wallData = 
        dynamic_cast<ImmersedWallData2D<T>*>( container->getData() );
    PLB_ASSERT(wallData);
    std::vector< Array<T,2> > const& vertices = wallData->vertices;
    std::vector< Array<T,2> > const& g = wallData->g;
    std::vector<int> const& flags = wallData->flags;
    PLB_ASSERT( vertices.size()==g.size() );
    PLB_ASSERT( vertices.size()==flags.size() );

    for (pluint i=0; i<vertices.size(); ++i) {
        Array<T,2> vertex = vertices[i];
        if ( flags[i]==reductionFlag &&
             closedOpenContained(vertex, domain) )
        {
            this->getStatistics().gatherSum(sum_g_ids[0], g[i][0]);
            this->getStatistics().gatherSum(sum_g_ids[1], g[i][1]);
        }
    }
}

template<typename T>
ReduceImmersedForce2D<T>* ReduceImmersedForce2D<T>::clone() const {
    return new ReduceImmersedForce2D<T>(*this);
}

template<typename T>
void ReduceImmersedForce2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing; // Container Block.
}

template<typename T>
BlockDomain::DomainT ReduceImmersedForce2D<T>::appliesTo() const {
    return BlockDomain::bulk;
}

template<typename T>
Array<T,2> ReduceImmersedForce2D<T>::getSumG() const {
    return Array<T,2> (
            this->getStatistics().getSum(sum_g_ids[0]),
            this->getStatistics().getSum(sum_g_ids[1]) );
}

/* ******** ReduceImmersedArea2D ************************************ */

template<typename T>
ReduceImmersedArea2D<T>::ReduceImmersedArea2D(int reductionFlag_)
    : sum_area_id(this->getStatistics().subscribeSum()),
      reductionFlag(reductionFlag_)
{ }

template<typename T>
void ReduceImmersedArea2D<T>::processGenericBlocks (
        Box2D domain, std::vector<AtomicBlock2D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==1 );
    AtomicContainerBlock2D* container = dynamic_cast<AtomicContainerBlock2D*>(blocks[0]);
    PLB_ASSERT( container );

    ImmersedWallData2D<T>* wallData = 
        dynamic_cast<ImmersedWallData2D<T>*>( container->getData() );
    PLB_ASSERT(wallData);
    std::vector< Array<T,2> > const& vertices = wallData->vertices;
    std::vector<T> const& areas = wallData->areas;
    std::vector<int> const& flags = wallData->flags;
    PLB_ASSERT( vertices.size()==areas.size() );
    PLB_ASSERT( vertices.size()==flags.size() );

    for (pluint i=0; i<vertices.size(); ++i) {
        Array<T,2> vertex = vertices[i];
        if ( flags[i]==reductionFlag &&
             closedOpenContained(vertex, domain) )
        {
            this->getStatistics().gatherSum(sum_area_id, areas[i]);
        }
    }
}

template<typename T>
ReduceImmersedArea2D<T>* ReduceImmersedArea2D<T>::clone() const {
    return new ReduceImmersedArea2D<T>(*this);
}

template<typename T>
void ReduceImmersedArea2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing; // Container Block.
}

template<typename T>
BlockDomain::DomainT ReduceImmersedArea2D<T>::appliesTo() const {
    return BlockDomain::bulk;
}

template<typename T>
T ReduceImmersedArea2D<T>::getSumArea() const {
    return this->getStatistics().getSum(sum_area_id);
}

/* ******** InamuroIteration2D ************************************ */

template<typename T, class VelFunction>
InamuroIteration2D<T,VelFunction>::InamuroIteration2D(VelFunction velFunction_, T tau_)
    : velFunction(velFunction_),
      tau(tau_)
{ }

template<typename T, class VelFunction>
void InamuroIteration2D<T,VelFunction>::processGenericBlocks (
        Box2D domain, std::vector<AtomicBlock2D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==3 );
    ScalarField2D<T>* rhoBar = dynamic_cast<ScalarField2D<T>*>(blocks[0]);
    TensorField2D<T,2>* j = dynamic_cast<TensorField2D<T,2>*>(blocks[1]);
    AtomicContainerBlock2D* container = dynamic_cast<AtomicContainerBlock2D*>(blocks[2]);
    PLB_ASSERT( rhoBar );
    PLB_ASSERT( j );
    PLB_ASSERT( container );

    Dot2D ofsJ = computeRelativeDisplacement(*rhoBar, *j);

    ImmersedWallData2D<T>* wallData = 
        dynamic_cast<ImmersedWallData2D<T>*>( container->getData() );
    PLB_ASSERT(wallData);
    Array<T,2> absOffset = wallData->offset;

    std::vector< Array<T,2> > const& vertices = wallData->vertices;
    std::vector<T> const& areas = wallData->areas;
    PLB_ASSERT( vertices.size()==areas.size() );
    std::vector<Array<T,2> > deltaG(vertices.size());
    std::vector<Array<T,2> >& g = wallData->g;
    PLB_ASSERT( vertices.size()==g.size() );

    // In this iteration, the force is computed for every vertex.
    for (pluint i=0; i<vertices.size(); ++i) {
        Array<T,2> const& vertex = vertices[i];
        Array<plint,2> intPos ( (plint)vertex[0], (plint)vertex[1] );
        Array<T,2> averageJ; averageJ.resetToZero();
        T averageRhoBar = T();
        // Use the weighting function to compute the average momentum
        // and the average density on the surface vertex.
        // x   x . x   x
        for (plint dx=-1; dx<=+2; ++dx) {
            for (plint dy=-1; dy<=+2; ++dy) {
                    Array<plint,2> pos(intPos+Array<plint,2>(dx,dy));
                    T nextRhoBar = rhoBar->get(pos[0], pos[1]);
                    Array<T,2> nextJ = j->get(pos[0]+ofsJ.x, pos[1]+ofsJ.y);
                    Array<T,2> r(pos[0]-vertex[0],pos[1]-vertex[1]);
                    T W = inamuroDeltaFunction<T>().W(r);
                    averageJ += W*nextJ;
                    averageRhoBar += W*nextRhoBar;
            }
        }
        //averageJ += 0.5*g[i];
        Array<T,2> wallVelocity = velFunction(vertex+absOffset);
        deltaG[i] = areas[i]*((averageRhoBar+1.)*wallVelocity-averageJ);
        //g[i] += deltaG[i];
        g[i] += deltaG[i]/(1.0+averageRhoBar);
    }
    
    // In this iteration, the force is applied from every vertex to the grid nodes.
    for (pluint i=0; i<vertices.size(); ++i) {
        Array<T,2> const& vertex = vertices[i];
        Array<plint,2> intPos ( (plint)vertex[0], (plint)vertex[1] );
        for (plint dx=-1; dx<=+2; ++dx) {
            for (plint dy=-1; dy<=+2; ++dy) {
                    Array<plint,2> pos(intPos+Array<plint,2>(dx,dy));
                    Array<T,2> nextJ = j->get(pos[0]+ofsJ.x, pos[1]+ofsJ.y);
                    Array<T,2> r(pos[0]-vertex[0],pos[1]-vertex[1]);
                    T W = inamuroDeltaFunction<T>().W(r);
                    nextJ += tau*W*deltaG[i];
                    j->get(pos[0]+ofsJ.x, pos[1]+ofsJ.y) = nextJ;
            }
        }
    }
}

template<typename T, class VelFunction>
InamuroIteration2D<T,VelFunction>* InamuroIteration2D<T,VelFunction>::clone() const {
    return new InamuroIteration2D<T,VelFunction>(*this);
}

template<typename T, class VelFunction>
void InamuroIteration2D<T,VelFunction>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;          // RhoBar
    modified[1] = modif::staticVariables;  // J
    modified[2] = modif::nothing;          // Container Block with triangle data.
}

template<typename T, class VelFunction>
BlockDomain::DomainT InamuroIteration2D<T,VelFunction>::appliesTo() const {
    return BlockDomain::bulk;
}

/* ******** IndexedInamuroIteration2D ************************************ */

template<typename T, class VelFunction>
IndexedInamuroIteration2D<T,VelFunction>::IndexedInamuroIteration2D(VelFunction velFunction_, T tau_)
    : velFunction(velFunction_),
      tau(tau_)
{ }

template<typename T, class VelFunction>
void IndexedInamuroIteration2D<T,VelFunction>::processGenericBlocks (
        Box2D domain, std::vector<AtomicBlock2D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==3 );
    ScalarField2D<T>* rhoBar = dynamic_cast<ScalarField2D<T>*>(blocks[0]);
    TensorField2D<T,2>* j = dynamic_cast<TensorField2D<T,2>*>(blocks[1]);
    AtomicContainerBlock2D* container = dynamic_cast<AtomicContainerBlock2D*>(blocks[2]);
    PLB_ASSERT( rhoBar );
    PLB_ASSERT( j );
    PLB_ASSERT( container );

    Dot2D ofsJ = computeRelativeDisplacement(*rhoBar, *j);

    ImmersedWallData2D<T>* wallData = 
        dynamic_cast<ImmersedWallData2D<T>*>( container->getData() );
    PLB_ASSERT(wallData);

    std::vector< Array<T,2> > const& vertices = wallData->vertices;
    std::vector<T> const& areas = wallData->areas;
    PLB_ASSERT( vertices.size()==areas.size() );
    std::vector<Array<T,2> > deltaG(vertices.size());
    std::vector<Array<T,2> >& g = wallData->g;
    PLB_ASSERT( vertices.size()==g.size() );
    std::vector<pluint> const& globalVertexIds = wallData->globalVertexIds;
    PLB_ASSERT( vertices.size()==globalVertexIds.size() );

    for (pluint i=0; i<vertices.size(); ++i) {
        Array<T,2> const& vertex = vertices[i];
        Array<plint,2> intPos ( (plint)vertex[0], (plint)vertex[1] );
        Array<T,2> averageJ; averageJ.resetToZero();
        T averageRhoBar = T();
        // x   x . x   x
        for (plint dx=-1; dx<=+2; ++dx) {
            for (plint dy=-1; dy<=+2; ++dy) {
                    Array<plint,2> pos(intPos+Array<plint,2>(dx,dy));
                    T nextRhoBar = rhoBar->get(pos[0], pos[1]);
                    Array<T,2> nextJ = j->get(pos[0]+ofsJ.x, pos[1]+ofsJ.y);
                    Array<T,2> r(pos[0]-vertex[0],pos[1]-vertex[1]);
                    T W = inamuroDeltaFunction<T>().W(r);
                    averageJ += W*nextJ;
                    averageRhoBar += W*nextRhoBar;
            }
        }
        //averageJ += 0.5*g[i];
        Array<T,2> wallVelocity = velFunction(globalVertexIds[i]);
        deltaG[i] = areas[i]*((averageRhoBar+1.)*wallVelocity-averageJ);
        //g[i] += deltaG[i];
        g[i] += deltaG[i]/(1.0+averageRhoBar);
    }
    
    for (pluint i=0; i<vertices.size(); ++i) {
        Array<T,2> const& vertex = vertices[i];
        Array<plint,2> intPos ( (plint)vertex[0], (plint)vertex[1] );
        for (plint dx=-1; dx<=+2; ++dx) {
            for (plint dy=-1; dy<=+2; ++dy) {
                    Array<plint,2> pos(intPos+Array<plint,2>(dx,dy));
                    Array<T,2> nextJ = j->get(pos[0]+ofsJ.x, pos[1]+ofsJ.y);
                    Array<T,2> r(pos[0]-vertex[0],pos[1]-vertex[1]);
                    T W = inamuroDeltaFunction<T>().W(r);
                    nextJ += tau*W*deltaG[i];
                    j->get(pos[0]+ofsJ.x, pos[1]+ofsJ.y) = nextJ;
            }
        }
    }
}

template<typename T, class VelFunction>
IndexedInamuroIteration2D<T,VelFunction>* IndexedInamuroIteration2D<T,VelFunction>::clone() const {
    return new IndexedInamuroIteration2D<T,VelFunction>(*this);
}

template<typename T, class VelFunction>
void IndexedInamuroIteration2D<T,VelFunction>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;          // RhoBar
    modified[1] = modif::staticVariables;  // J
    modified[2] = modif::nothing;          // Container Block with triangle data.
}

template<typename T, class VelFunction>
BlockDomain::DomainT IndexedInamuroIteration2D<T,VelFunction>::appliesTo() const {
    return BlockDomain::bulk;
}

/* ******** ConstVelInamuroIteration2D ************************************ */

template<typename T>
ConstVelInamuroIteration2D<T>::ConstVelInamuroIteration2D(Array<T,2> const& wallVelocity_, T tau_)
    : wallVelocity(wallVelocity_),
      tau(tau_)
{ }

template<typename T>
void ConstVelInamuroIteration2D<T>::processGenericBlocks (
        Box2D domain, std::vector<AtomicBlock2D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==3 );
    ScalarField2D<T>* rhoBar = dynamic_cast<ScalarField2D<T>*>(blocks[0]);
    TensorField2D<T,2>* j = dynamic_cast<TensorField2D<T,2>*>(blocks[1]);
    AtomicContainerBlock2D* container = dynamic_cast<AtomicContainerBlock2D*>(blocks[2]);
    PLB_ASSERT( rhoBar );
    PLB_ASSERT( j );
    PLB_ASSERT( container );

    Dot2D ofsJ = computeRelativeDisplacement(*rhoBar, *j);

    ImmersedWallData2D<T>* wallData = 
        dynamic_cast<ImmersedWallData2D<T>*>( container->getData() );
    PLB_ASSERT(wallData);
    std::vector< Array<T,2> > const& vertices = wallData->vertices;
    std::vector<T> const& areas = wallData->areas;
    PLB_ASSERT( vertices.size()==areas.size() );
    std::vector<Array<T,2> > deltaG(vertices.size());
    std::vector<Array<T,2> >& g = wallData->g;
    PLB_ASSERT( vertices.size()==g.size() );

    for (pluint i=0; i<vertices.size(); ++i) {
        Array<T,2> const& vertex = vertices[i];
        Array<plint,2> intPos ( (plint)vertex[0], (plint)vertex[1] );
        Array<T,2> averageJ; averageJ.resetToZero();
        T averageRhoBar = T();
        // x   x . x   x
        for (plint dx=-1; dx<=+2; ++dx) {
            for (plint dy=-1; dy<=+2; ++dy) {
                    Array<plint,2> pos(intPos+Array<plint,2>(dx,dy));
                    T nextRhoBar = rhoBar->get(pos[0], pos[1]);
                    Array<T,2> nextJ = j->get(pos[0]+ofsJ.x, pos[1]+ofsJ.y);
                    Array<T,2> r(pos[0]-vertex[0],pos[1]-vertex[1]);
                    T W = inamuroDeltaFunction<T>().W(r);
                    averageJ += W*nextJ;
                    averageRhoBar += W*nextRhoBar;
            }
        }
        //averageJ += 0.5*g[i];
        deltaG[i] = areas[i]*((averageRhoBar+1.)*wallVelocity-averageJ);
        //g[i] += deltaG[i];
        g[i] += deltaG[i]/(1.0+averageRhoBar);
    }
    
    for (pluint i=0; i<vertices.size(); ++i) {
        Array<T,2> const& vertex = vertices[i];
        Array<plint,2> intPos ( (plint)vertex[0], (plint)vertex[1] );
        for (plint dx=-1; dx<=+2; ++dx) {
            for (plint dy=-1; dy<=+2; ++dy) {
                    Array<plint,2> pos(intPos+Array<plint,2>(dx,dy));
                    Array<T,2> nextJ = j->get(pos[0]+ofsJ.x, pos[1]+ofsJ.y);
                    Array<T,2> r(pos[0]-vertex[0],pos[1]-vertex[1]);
                    T W = inamuroDeltaFunction<T>().W(r);
                    nextJ += tau*W*deltaG[i];
                    j->get(pos[0]+ofsJ.x, pos[1]+ofsJ.y) = nextJ;
            }
        }
    }
}

template<typename T>
ConstVelInamuroIteration2D<T>* ConstVelInamuroIteration2D<T>::clone() const {
    return new ConstVelInamuroIteration2D<T>(*this);
}

template<typename T>
void ConstVelInamuroIteration2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;          // RhoBar
    modified[1] = modif::staticVariables;  // J
    modified[2] = modif::nothing;          // Container Block with triangle data.
}

template<typename T>
BlockDomain::DomainT ConstVelInamuroIteration2D<T>::appliesTo() const {
    return BlockDomain::bulk;
}



/* ******** ComputeImmersedBoundaryForce2D ************************************ */

template<typename T>
void ComputeImmersedBoundaryForce2D<T>::processGenericBlocks (
        Box2D domain, std::vector<AtomicBlock2D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==2 );
    TensorField2D<T,2>* force = dynamic_cast<TensorField2D<T,2>*>(blocks[0]);
    AtomicContainerBlock2D* container = dynamic_cast<AtomicContainerBlock2D*>(blocks[1]);
    PLB_ASSERT(force);
    PLB_ASSERT(container);

    ImmersedWallData2D<T>* wallData = 
        dynamic_cast<ImmersedWallData2D<T>*>( container->getData() );
    PLB_ASSERT(wallData);
    std::vector< Array<T,2> > const& vertices = wallData->vertices;
    std::vector<Array<T,2> >& g = wallData->g;
    PLB_ASSERT( vertices.size()==g.size() );

    for (plint iX = 0; iX < force->getNx(); iX++) {
        for (plint iY = 0; iY < force->getNy(); iY++) {
            for (plint iZ = 0; iZ < force->getNz(); iZ++) {
                force->get(iX, iY, iZ).resetToZero();
            }
        }
    }

    for (pluint i=0; i<vertices.size(); ++i) {
        Array<T,2> const& vertex = vertices[i];
        Array<plint,2> intPos ( (plint)vertex[0], (plint)vertex[1] );
        for (plint dx=-1; dx<=+2; ++dx) {
            for (plint dy=-1; dy<=+2; ++dy) {
                    Array<plint,2> pos(intPos+Array<plint,2>(dx,dy));
                    Array<T,2> r(pos[0]-vertex[0],pos[1]-vertex[1]);
                    T W = inamuroDeltaFunction<T>().W(r);
                    force->get(pos[0], pos[1]) += W*g[i];
            }
        }
    }
}

template<typename T>
ComputeImmersedBoundaryForce2D<T>* ComputeImmersedBoundaryForce2D<T>::clone() const {
    return new ComputeImmersedBoundaryForce2D<T>(*this);
}

template<typename T>
void ComputeImmersedBoundaryForce2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;  // Force
    modified[1] = modif::nothing;          // Container Block with triangle data.
}

template<typename T>
BlockDomain::DomainT ComputeImmersedBoundaryForce2D<T>::appliesTo() const {
    return BlockDomain::bulk;
}


/* ******** InstantiateImmersedWallData2D ************************************ */

template<typename T>
InstantiateImmersedWallData2D<T>::InstantiateImmersedWallData2D (
            std::vector< Array<T,2> > const& vertices_,
            std::vector<T> const& areas_,
            std::vector< Array<T,2> > const& normals_)
    : vertices(vertices_),
      areas(areas_),
      normals(normals_)
{
    PLB_ASSERT(vertices.size() == areas.size());
    PLB_ASSERT(normals.size()==0 || normals.size() == areas.size());
}

template<typename T>
void InstantiateImmersedWallData2D<T>::processGenericBlocks (
        Box2D domain, std::vector<AtomicBlock2D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==1 );
    AtomicContainerBlock2D* container = dynamic_cast<AtomicContainerBlock2D*>(blocks[0]);
    PLB_ASSERT( container );
    bool useNormals = normals.size()>0;
    Dot2D location = container->getLocation();
    Array<T,2> offset(location.x,location.y);
    ImmersedWallData2D<T>* wallData = new ImmersedWallData2D<T>;
    Box2D extendedEnvelope(domain.enlarge(2));

    for (pluint i=0; i<vertices.size(); ++i) {
        Array<T,2> vertex = vertices[i]-offset;
        if (contained(vertex, extendedEnvelope)) {
            wallData->vertices.push_back(vertex);
            wallData->areas.push_back(areas[i]);
            if (useNormals) {
                wallData->normals.push_back(normals[i]);
            }
            wallData->g.push_back(Array<T,2>(0.,0.,0.));
            wallData->globalVertexIds.push_back(i);
        }
    }
    wallData->flags = std::vector<int>(wallData->vertices.size(), 0);
    wallData->offset = offset;
    container->setData(wallData);
}

template<typename T>
InstantiateImmersedWallData2D<T>* InstantiateImmersedWallData2D<T>::clone() const {
    return new InstantiateImmersedWallData2D<T>(*this);
}

template<typename T>
void InstantiateImmersedWallData2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;  // Container Block with triangle data.
}

template<typename T>
BlockDomain::DomainT InstantiateImmersedWallData2D<T>::appliesTo() const {
    return BlockDomain::bulk;
}


/* ******** InstantiateImmersedWallDataWithTagging2D ************************************ */

template<typename T>
InstantiateImmersedWallDataWithTagging2D<T>::InstantiateImmersedWallDataWithTagging2D (
            std::vector< Array<T,2> > const& vertices_,
            std::vector<T> const& areas_, int fluidFlag_ )
    : vertices(vertices_),
      areas(areas_),
      fluidFlag(fluidFlag_)
{
    PLB_ASSERT(vertices.size() == areas.size());
}

template<typename T>
void InstantiateImmersedWallDataWithTagging2D<T>::processGenericBlocks (
        Box2D domain, std::vector<AtomicBlock2D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==2 );
    AtomicContainerBlock2D* container = dynamic_cast<AtomicContainerBlock2D*>(blocks[0]);
    PLB_ASSERT( container );
    Dot2D location = container->getLocation();
    Array<T,2> offset(location.x,location.y);

    ScalarField2D<int>* flagMatrix = dynamic_cast<ScalarField2D<int>*>(blocks[1]);
    PLB_ASSERT(flagMatrix);
    Dot2D ofsFlag = computeRelativeDisplacement(*container, *flagMatrix);
    Array<plint,2> flagDispl(ofsFlag.x,ofsFlag.y);

    ImmersedWallData2D<T>* wallData = new ImmersedWallData2D<T>;
    Box2D extendedEnvelope(domain.enlarge(2));

    for (pluint i=0; i<vertices.size(); ++i) {
        Array<T,2> vertex = vertices[i]-offset;
        if (contained(vertex, extendedEnvelope)) {
            wallData->vertices.push_back(vertex);
            wallData->areas.push_back(areas[i]);
            wallData->g.push_back(Array<T,2>(0.,0.,0.));
            wallData->globalVertexIds.push_back(i);
            Array<plint,2> intPos ( (plint)vertex[0], (plint)vertex[1] );
            bool hasFluidNeighbor=false;
            for (plint dx=-1; dx<=+2; ++dx) {
                for (plint dy=-1; dy<=+2; ++dy) {
                        Array<plint,2> pos(intPos+Array<plint,2>(dx,dy)+flagDispl);
                        if(flagMatrix->get(pos[0],pos[1])==fluidFlag) {
                            hasFluidNeighbor=true;
                    }
                }
            }
            if (hasFluidNeighbor) {
                wallData->flags.push_back(0);
            }
            else {
                wallData->flags.push_back(1);
            }
        }
    }
    wallData->offset = offset;
    container->setData(wallData);
}

template<typename T>
InstantiateImmersedWallDataWithTagging2D<T>* InstantiateImmersedWallDataWithTagging2D<T>::clone() const {
    return new InstantiateImmersedWallDataWithTagging2D<T>(*this);
}

template<typename T>
void InstantiateImmersedWallDataWithTagging2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;  // Container Block with triangle data.
    modified[1] = modif::nothing;  // Flag matrix.
}

template<typename T>
BlockDomain::DomainT InstantiateImmersedWallDataWithTagging2D<T>::appliesTo() const {
    return BlockDomain::bulk;
}

/* ******** InstantiateImmersedWallDataWithIndexedTagging2D ************************************ */

template<typename T>
InstantiateImmersedWallDataWithIndexedTagging2D<T>::InstantiateImmersedWallDataWithIndexedTagging2D (
            std::vector< Array<T,2> > const& vertices_,
            std::vector<T> const& areas_, std::vector<int> const& flags_ )
    : vertices(vertices_),
      areas(areas_),
      flags(flags_)
{
    PLB_ASSERT(vertices.size() == areas.size());
    PLB_ASSERT(vertices.size() == flags.size());
}

template<typename T>
void InstantiateImmersedWallDataWithIndexedTagging2D<T>::processGenericBlocks (
        Box2D domain, std::vector<AtomicBlock2D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==1 );
    AtomicContainerBlock2D* container = dynamic_cast<AtomicContainerBlock2D*>(blocks[0]);
    PLB_ASSERT( container );
    Dot2D location = container->getLocation();
    Array<T,2> offset(location.x,location.y);

    ImmersedWallData2D<T>* wallData = new ImmersedWallData2D<T>;
    Box2D extendedEnvelope(domain.enlarge(2));

    for (pluint i=0; i<vertices.size(); ++i) {
        Array<T,2> vertex = vertices[i]-offset;
        if (contained(vertex, extendedEnvelope)) {
            wallData->vertices.push_back(vertex);
            wallData->areas.push_back(areas[i]);
            wallData->g.push_back(Array<T,2>(0.,0.,0.));
            wallData->flags.push_back(flags[i]);
            wallData->globalVertexIds.push_back(i);
        }
    }
    wallData->offset = offset;
    container->setData(wallData);
}

template<typename T>
InstantiateImmersedWallDataWithIndexedTagging2D<T>* InstantiateImmersedWallDataWithIndexedTagging2D<T>::clone() const {
    return new InstantiateImmersedWallDataWithIndexedTagging2D<T>(*this);
}

template<typename T>
void InstantiateImmersedWallDataWithIndexedTagging2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;  // Container Block with triangle data.
}

template<typename T>
BlockDomain::DomainT InstantiateImmersedWallDataWithIndexedTagging2D<T>::appliesTo() const {
    return BlockDomain::bulk;
}

/* ******** TwoPhaseInamuroParam2D ************************************ */

template<typename T>
T TwoPhaseInamuroParam2D<T>::area(plint i) const {
    PLB_ASSERT((pluint)i < numVertices);
    return wallData->areas[i];
}

template<typename T>
Array<T,2>& TwoPhaseInamuroParam2D<T>::g(plint i) {
    PLB_ASSERT((pluint)i < numVertices);
    return wallData->g[i];
}

template<typename T>
Array<T,2> TwoPhaseInamuroParam2D<T>::vertex(plint i) const {
    PLB_ASSERT((pluint)i < numVertices);
    return wallData->vertices[i];
}

template<typename T>
Array<T,2> TwoPhaseInamuroParam2D<T>::absoluteVertex(plint i) const {
    PLB_ASSERT((pluint)i < numVertices);
    return wallData->vertices[i] + absOffset;
}

template<typename T>
Array<plint,2> TwoPhaseInamuroParam2D<T>::intVertex(plint i) const {
    PLB_ASSERT((pluint)i < numVertices);
    Array<T,2> vertex = wallData->vertices[i];
    return Array<plint,2>((plint)vertex[0], (plint)vertex[1]);
}

template<typename T>
T TwoPhaseInamuroParam2D<T>::rhoBar(plint iX, plint iY) const {
    int flag = getFlag(iX,iY);
    if (flag==twoPhaseFlag::empty) {
        return rhoBar2_->get(iX+ofsRhoBar2.x,iY+ofsRhoBar2.y);
    }
    else {
        return rhoBar_->get(iX,iY);
    }
}

template<typename T>
Array<T,2> TwoPhaseInamuroParam2D<T>::j(plint iX, plint iY) const {
    int flag = getFlag(iX,iY);
    if (flag==twoPhaseFlag::empty) {
        return j2_->get(iX+ofsJ2.x,iY+ofsJ2.y);
    }
    else {
        return j_->get(iX+ofsJ.x,iY+ofsJ.y);
    }
}

template<typename T>
void TwoPhaseInamuroParam2D<T>::addToJ(plint iX, plint iY, Array<T,2> deltaJ)
{
    int flag = getFlag(iX,iY);
    if (flag==twoPhaseFlag::interface) {
        j_->get(iX+ofsJ.x,iY+ofsJ.y) += deltaJ;
        j2_->get(iX+ofsJ2.x,iY+ofsJ2.y) += deltaJ;
    }
    else if (flag==twoPhaseFlag::empty) {
        j2_->get(iX+ofsJ2.x,iY+ofsJ2.y) += deltaJ;
    }
    else {
        j_->get(iX+ofsJ.x,iY+ofsJ.y) += deltaJ;
    }
}

template<typename T>
int TwoPhaseInamuroParam2D<T>::getFlag(plint iX, plint iY) const {
    return flag_->get(iX+ofsFlag.x,iY+ofsFlag.y);
}

template<typename T>
pluint TwoPhaseInamuroParam2D<T>::getGlobalVertexId(plint i) const {
    PLB_ASSERT((pluint)i < numVertices);
    return wallData->globalVertexIds[i];
}

template<typename T>
T TwoPhaseInamuroParam2D<T>::getTau(plint iX, plint iY) const {
    T vf = volumeFraction_->get(iX+ofsVF.x,iY+ofsVF.y);
    return tau*vf + tau2*(1.-vf);
}

template<typename T>
TwoPhaseInamuroParam2D<T>::TwoPhaseInamuroParam2D(std::vector<AtomicBlock2D*>& blocks, T tau_, T tau2_)
    : tau(tau_), tau2(tau2_)
{
    PLB_PRECONDITION( blocks.size()==7 );
    rhoBar_ = dynamic_cast<ScalarField2D<T>*>(blocks[0]);
    rhoBar2_ = dynamic_cast<ScalarField2D<T>*>(blocks[1]);
    j_ = dynamic_cast<TensorField2D<T,2>*>(blocks[2]);
    j2_ = dynamic_cast<TensorField2D<T,2>*>(blocks[3]);
    flag_ = dynamic_cast<ScalarField2D<int>*>(blocks[4]);
    volumeFraction_ = dynamic_cast<ScalarField2D<T>*>(blocks[5]);
    container = dynamic_cast<AtomicContainerBlock2D*>(blocks[6]);

    PLB_ASSERT( rhoBar_ );
    PLB_ASSERT( rhoBar2_ );
    PLB_ASSERT( j_ );
    PLB_ASSERT( j2_ );
    PLB_ASSERT( flag_ );
    PLB_ASSERT( volumeFraction_ );
    PLB_ASSERT( container );

    ofsRhoBar2 = computeRelativeDisplacement(*rhoBar_, *rhoBar2_);
    ofsJ = computeRelativeDisplacement(*rhoBar_, *j_);
    ofsJ2 = computeRelativeDisplacement(*rhoBar_, *j2_);
    ofsFlag = computeRelativeDisplacement(*rhoBar_, *flag_);
    ofsVF = computeRelativeDisplacement(*rhoBar_, *volumeFraction_);

    wallData = dynamic_cast<ImmersedWallData2D<T>*>( container->getData() );
    PLB_ASSERT(wallData);
    absOffset = wallData->offset;

    numVertices = wallData->vertices.size();
    PLB_ASSERT( numVertices == wallData->areas.size() );
    PLB_ASSERT( numVertices == wallData->g.size() );
}

/* ******** TwoPhaseInamuroIteration2D ************************************ */

template<typename T, class VelFunction>
TwoPhaseInamuroIteration2D<T,VelFunction>::TwoPhaseInamuroIteration2D(VelFunction velFunction_, T tau_, T tau2_)
    : velFunction(velFunction_),
      tau(tau_),
      tau2(tau2_)
{ }

template<typename T, class VelFunction>
void TwoPhaseInamuroIteration2D<T,VelFunction>::processGenericBlocks (
        Box2D domain, std::vector<AtomicBlock2D*> blocks )
{
    TwoPhaseInamuroParam2D<T> param(blocks, tau, tau2);
    std::vector<Array<T,2> > deltaG(param.getNumVertices());

    for (pluint i=0; i<param.getNumVertices(); ++i) {
        Array<T,2> vertex(param.vertex(i));
        Array<plint,2> intPos(param.intVertex(i));

        Array<T,2> averageJ; averageJ.resetToZero();
        T averageRhoBar = T();
        // x   x . x   x
        for (plint dx=-1; dx<=+2; ++dx) {
            for (plint dy=-1; dy<=+2; ++dy) {
                    Array<plint,2> pos(intPos+Array<plint,2>(dx,dy));
                    T nextRhoBar = param.rhoBar(pos[0], pos[1]);
                    Array<T,2> nextJ = param.j(pos[0], pos[1]);
                    Array<T,2> r(pos[0]-vertex[0],pos[1]-vertex[1]);
                    T W = inamuroDeltaFunction<T>().W(r);
                    averageJ += W*nextJ;
                    averageRhoBar += W*nextRhoBar;
            }
        }
        //averageJ += 0.5*param.g(i);
        Array<T,2> wallVelocity = velFunction(param.absoluteVertex(i));
        deltaG[i] = param.area(i)*((averageRhoBar+1.)*wallVelocity-averageJ);
        //param.g(i) += deltaG[i];
        param.g(i) += deltaG[i]/(1.0+averageRhoBar);
    }
    
    for (pluint i=0; i<param.getNumVertices(); ++i) {
        Array<T,2> vertex(param.vertex(i));
        Array<plint,2> intPos(param.intVertex(i));

        for (plint dx=-1; dx<=+2; ++dx) {
            for (plint dy=-1; dy<=+2; ++dy) {
                    Array<plint,2> pos(intPos+Array<plint,2>(dx,dy));
                    Array<T,2> r(pos[0]-vertex[0],pos[1]-vertex[1]);
                    T W = inamuroDeltaFunction<T>().W(r);
                    param.addToJ(pos[0],pos[1], param.getTau(pos[0],pos[1])*W*deltaG[i]);
            }
        }
    }
}

template<typename T, class VelFunction>
TwoPhaseInamuroIteration2D<T,VelFunction>* TwoPhaseInamuroIteration2D<T,VelFunction>::clone() const {
    return new TwoPhaseInamuroIteration2D<T,VelFunction>(*this);
}

template<typename T, class VelFunction>
void TwoPhaseInamuroIteration2D<T,VelFunction>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;          // RhoBar
    modified[1] = modif::nothing;          // RhoBar2
    modified[2] = modif::staticVariables;  // j
    modified[3] = modif::staticVariables;  // j2
    modified[4] = modif::nothing;          // flag
    modified[5] = modif::nothing;          // volume fraction
    modified[6] = modif::nothing;          // Container Block with triangle data.
}

template<typename T, class VelFunction>
BlockDomain::DomainT TwoPhaseInamuroIteration2D<T,VelFunction>::appliesTo() const {
    return BlockDomain::bulk;
}

/* ******** TwoPhaseIndexedInamuroIteration2D ************************************ */

template<typename T, class VelFunction>
TwoPhaseIndexedInamuroIteration2D<T,VelFunction>::TwoPhaseIndexedInamuroIteration2D(VelFunction velFunction_, T tau_, T tau2_)
    : velFunction(velFunction_),
      tau(tau_),
      tau2(tau2_)
{ }

template<typename T, class VelFunction>
void TwoPhaseIndexedInamuroIteration2D<T,VelFunction>::processGenericBlocks (
        Box2D domain, std::vector<AtomicBlock2D*> blocks )
{
    TwoPhaseInamuroParam2D<T> param(blocks, tau, tau2);
    std::vector<Array<T,2> > deltaG(param.getNumVertices());

    for (pluint i=0; i<param.getNumVertices(); ++i) {
        Array<T,2> vertex(param.vertex(i));
        Array<plint,2> intPos(param.intVertex(i));

        Array<T,2> averageJ; averageJ.resetToZero();
        T averageRhoBar = T();
        // x   x . x   x
        for (plint dx=-1; dx<=+2; ++dx) {
            for (plint dy=-1; dy<=+2; ++dy) {
                    Array<plint,2> pos(intPos+Array<plint,2>(dx,dy));
                    T nextRhoBar = param.rhoBar(pos[0], pos[1]);
                    Array<T,2> nextJ = param.j(pos[0], pos[1]);
                    Array<T,2> r(pos[0]-vertex[0],pos[1]-vertex[1]);
                    T W = inamuroDeltaFunction<T>().W(r);
                    averageJ += W*nextJ;
                    averageRhoBar += W*nextRhoBar;
            }
        }
        //averageJ += 0.5*param.g(i);
        Array<T,2> wallVelocity = velFunction(param.getGlobalVertexId(i));
        deltaG[i] = param.area(i)*((averageRhoBar+1.)*wallVelocity-averageJ);
        //param.g(i) += deltaG[i];
        param.g(i) += deltaG[i]/(1.0+averageRhoBar);
    }
    
    for (pluint i=0; i<param.getNumVertices(); ++i) {
        Array<T,2> vertex(param.vertex(i));
        Array<plint,2> intPos(param.intVertex(i));

        for (plint dx=-1; dx<=+2; ++dx) {
            for (plint dy=-1; dy<=+2; ++dy) {
                    Array<plint,2> pos(intPos+Array<plint,2>(dx,dy));
                    Array<T,2> r(pos[0]-vertex[0],pos[1]-vertex[1]);
                    T W = inamuroDeltaFunction<T>().W(r);
                    param.addToJ(pos[0],pos[1], param.getTau(pos[0],pos[1])*W*deltaG[i]);
            }
        }
    }
}

template<typename T, class VelFunction>
TwoPhaseIndexedInamuroIteration2D<T,VelFunction>* TwoPhaseIndexedInamuroIteration2D<T,VelFunction>::clone() const {
    return new TwoPhaseIndexedInamuroIteration2D<T,VelFunction>(*this);
}

template<typename T, class VelFunction>
void TwoPhaseIndexedInamuroIteration2D<T,VelFunction>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;          // RhoBar
    modified[1] = modif::nothing;          // RhoBar2
    modified[2] = modif::staticVariables;  // j
    modified[3] = modif::staticVariables;  // j2
    modified[4] = modif::nothing;          // flag
    modified[5] = modif::nothing;          // volume fraction
    modified[6] = modif::nothing;          // Container Block with triangle data.
}

template<typename T, class VelFunction>
BlockDomain::DomainT TwoPhaseIndexedInamuroIteration2D<T,VelFunction>::appliesTo() const {
    return BlockDomain::bulk;
}

/* ******** TwoPhaseConstVelInamuroIteration2D ************************************ */

template<typename T>
TwoPhaseConstVelInamuroIteration2D<T>::TwoPhaseConstVelInamuroIteration2D(Array<T,2> const& wallVelocity_, T tau_, T tau2_)
    : wallVelocity(wallVelocity_),
      tau(tau_),
      tau2(tau2_)
{ }

template<typename T>
void TwoPhaseConstVelInamuroIteration2D<T>::processGenericBlocks (
        Box2D domain, std::vector<AtomicBlock2D*> blocks )
{
    TwoPhaseInamuroParam2D<T> param(blocks, tau, tau2);
    std::vector<Array<T,2> > deltaG(param.getNumVertices());

    for (pluint i=0; i<param.getNumVertices(); ++i) {
        Array<T,2> vertex(param.vertex(i));
        Array<plint,2> intPos(param.intVertex(i));

        Array<T,2> averageJ; averageJ.resetToZero();
        T averageRhoBar = T();
        // x   x . x   x
        for (plint dx=-1; dx<=+2; ++dx) {
            for (plint dy=-1; dy<=+2; ++dy) {
                    Array<plint,2> pos(intPos+Array<plint,2>(dx,dy));
                    T nextRhoBar = param.rhoBar(pos[0], pos[1]);
                    Array<T,2> nextJ = param.j(pos[0], pos[1]);
                    Array<T,2> r(pos[0]-vertex[0],pos[1]-vertex[1]);
                    T W = inamuroDeltaFunction<T>().W(r);
                    averageJ += W*nextJ;
                    averageRhoBar += W*nextRhoBar;
            }
        }
        //averageJ += 0.5*param.g(i);
        deltaG[i] = param.area(i)*((averageRhoBar+1.)*wallVelocity-averageJ);
        //param.g(i) += deltaG[i];
        param.g(i) += deltaG[i]/(1.0+averageRhoBar);
    }
    
    for (pluint i=0; i<param.getNumVertices(); ++i) {
        Array<T,2> vertex(param.vertex(i));
        Array<plint,2> intPos(param.intVertex(i));

        for (plint dx=-1; dx<=+2; ++dx) {
            for (plint dy=-1; dy<=+2; ++dy) {
                    Array<plint,2> pos(intPos+Array<plint,2>(dx,dy));
                    Array<T,2> r(pos[0]-vertex[0],pos[1]-vertex[1]);
                    T W = inamuroDeltaFunction<T>().W(r);
                    param.addToJ(pos[0],pos[1], param.getTau(pos[0],pos[1])*W*deltaG[i]);
            }
        }
    }
}

template<typename T>
TwoPhaseConstVelInamuroIteration2D<T>* TwoPhaseConstVelInamuroIteration2D<T>::clone() const {
    return new TwoPhaseConstVelInamuroIteration2D<T>(*this);
}

template<typename T>
void TwoPhaseConstVelInamuroIteration2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;          // RhoBar
    modified[1] = modif::nothing;          // RhoBar2
    modified[2] = modif::staticVariables;  // j
    modified[3] = modif::staticVariables;  // j2
    modified[4] = modif::nothing;          // flag
    modified[5] = modif::nothing;          // Volume fraction
    modified[6] = modif::nothing;          // Container Block with triangle data.
}

template<typename T>
BlockDomain::DomainT TwoPhaseConstVelInamuroIteration2D<T>::appliesTo() const {
    return BlockDomain::bulk;
}

}  // namespace plb

#endif  // IMMERSED_WALLS_2D_HH

