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
 * The dynamics of a 2D block lattice -- header file.
 */
#ifndef BLOCK_LATTICE_2D_H
#define BLOCK_LATTICE_2D_H

#include "core/globalDefs.h"
#include "core/plbDebug.h"
#include "core/cell.h"
#include "atomicBlock/dataField2D.h"
#include "core/blockLatticeBase2D.h"
#include "atomicBlock/atomicBlock2D.h"
#include "core/blockIdentifiers.h"
#include <vector>
#include <map>

namespace plb {

template<typename T, template<typename U> class Descriptor> struct Dynamics;
template<typename T, template<typename U> class Descriptor> class BlockLattice2D;

/**
 * ATTENTION:
 * 1,能否将“std::vector<char>& buffer”移植到KERNEL空间执行？
 *    	答案是不太方便，由于Palabos在设计的时候使用多处使用std::vector<char>& buffer,
 * 如果使用 cudaManagedMemory 容易造成工作量较大，因此应该在BlockLatticeDataTransfer
 * 中使用cudaManagedMemory缓存需要交换的数据，避免在GPU和CPU之间复制BlockLattice中的
 * 所有数据。
 * 2,send-receive 和 attribute 策略有何区别？
 *
 * TODO:
 * 	1,需要BlockLatticeDataTransfer能够访问 BlockLattice2D 的私有数据成员
 * \ref BlockLattice2D::rawData 或者 \ref BlockLattice2D::grid;
 * 	2,BlockDataTransfer2D 需要成为 \ref BlockLattice2D 的友元类。
 * 	3,需要给 \ref BlockLatticeDataTransfer2D 添加 \ref PLB_KERNEL 函数;
 */
template<typename T, template<typename U> class Descriptor>
class BlockLatticeDataTransfer2D : public BlockDataTransfer2D {
public:
    BlockLatticeDataTransfer2D(BlockLattice2D<T,Descriptor>& lattice_);
    virtual plint staticCellSize() const;
    /// Send data from the lattice into a byte-stream.
    virtual void send(Box2D domain, std::vector<char>& buffer, modif::ModifT kind) const;
    /// Receive data from a byte-stream into the lattice.
    virtual void receive(Box2D domain, std::vector<char> const& buffer, modif::ModifT kind);
    virtual void receive(Box2D domain, std::vector<char> const& buffer, modif::ModifT kind, Dot2D offset) {
        receive(domain, buffer, kind);
    }
    /// Receive data from a byte-stream into the block, and re-map IDs for dynamics if exist.
    virtual void receive( Box2D domain, std::vector<char> const& buffer,
                          modif::ModifT kind, std::map<int,std::string> const& foreignIds );
    /// Attribute data between two lattices.
    virtual void attribute(Box2D toDomain, plint deltaX, plint deltaY,
                           AtomicBlock2D const& from, modif::ModifT kind);
    virtual void attribute(Box2D toDomain, plint deltaX, plint deltaY,
                           AtomicBlock2D const& from, modif::ModifT kind, Dot2D absoluteOffset)
    {
        attribute(toDomain, deltaX, deltaY, from, kind);
    }
private:
    void send_static(Box2D domain, std::vector<char>& buffer) const;
    void send_dynamic(Box2D domain, std::vector<char>& buffer) const;
    void send_all(Box2D domain, std::vector<char>& buffer) const;

    void receive_static(Box2D domain, std::vector<char> const& buffer);
    void receive_dynamic(Box2D domain, std::vector<char> const& buffer);
    void receive_all(Box2D domain, std::vector<char> const& buffer);
    void receive_regenerate(Box2D domain, std::vector<char> const& buffer,
                            std::map<int,int> const& idIndirect = (std::map<int,int>()) );

    void attribute_static (
        Box2D toDomain, plint deltaX, plint deltaY,
        BlockLattice2D<T,Descriptor> const& from );
    void attribute_dynamic (
        Box2D toDomain, plint deltaX, plint deltaY,
        BlockLattice2D<T,Descriptor> const& from );
    void attribute_all (
        Box2D toDomain, plint deltaX, plint deltaY,
        BlockLattice2D<T,Descriptor> const& from );
    void attribute_regenerate (
        Box2D toDomain, plint deltaX, plint deltaY,
        BlockLattice2D<T,Descriptor> const& from );
private:
    BlockLattice2D<T,Descriptor>& lattice;
};


/** A block lattice contains a regular array of Cell objects and
 * some useful methods to execute the LB dynamics on the lattice.
 *
 * This class is not intended to be derived from.
 * TODO: 1,查找该类的成员变量的在进行运算时的角色，确定哪一个需要复制到GPU上;
 * 	2,修改BlockLattice2D的部分函数，使其能够在GPU上对rawData进行操作;
 *
 */
template<typename T, template<typename U> class Descriptor>
class BlockLattice2D : public BlockLatticeBase2D<T,Descriptor>,
                       public AtomicBlock2D
{
public:
    /// Construction of an nx_ by ny_ lattice
    BlockLattice2D(plint nx_, plint ny_, Dynamics<T,Descriptor>* backgroundDynamics);
    /// Destruction of the lattice
    ~BlockLattice2D();
    /// Copy construction
    BlockLattice2D(BlockLattice2D<T,Descriptor> const& rhs);
    /// Copy assignment
    BlockLattice2D& operator=(BlockLattice2D<T,Descriptor> const& rhs);
    /// Swap the content of two BlockLattices
    void swap(BlockLattice2D& rhs);
public:
    /// Read/write access to lattice cells
    virtual Cell<T,Descriptor>& get(plint iX, plint iY) {
        PLB_PRECONDITION(iX<this->getNx());
        PLB_PRECONDITION(iY<this->getNy());
        return grid[iX][iY];
    }
    /// Read only access to lattice cells
    virtual Cell<T,Descriptor> const& get(plint iX, plint iY) const {
        PLB_PRECONDITION(iX<this->getNx());
        PLB_PRECONDITION(iY<this->getNy());
        return grid[iX][iY];
    }
    /// Specify wheter statistics measurements are done on given rect. domain
    virtual void specifyStatisticsStatus(Box2D domain, bool status);
    /// Apply collision step to a rectangular domain
    virtual void collide(Box2D domain);
    /// Apply collision step to the whole domain
    PLB_HOST_ONLY virtual void collide();
    /// Apply streaming step to a rectangular domain
    PLB_HOST_ONLY virtual void stream(Box2D domain);
    /// Apply streaming step to the whole domain
    PLB_HOST_ONLY virtual void stream();
    /// Apply first collision, then streaming step to a rectangular domain
    PLB_HOST_ONLY virtual void collideAndStream(Box2D domain);
    /// Apply first collision, then streaming step to the whole domain
    PLB_HOST_ONLY virtual void collideAndStream();
    /// Increment time counter
    virtual void incrementTime();
    /// Get access to data transfer between blocks
    virtual BlockLatticeDataTransfer2D<T,Descriptor>& getDataTransfer();
    /// Get access to data transfer between blocks (const version)
    virtual BlockLatticeDataTransfer2D<T,Descriptor> const& getDataTransfer() const;
public:
    /// Attribute dynamics to a cell.
    void attributeDynamics(plint iX, plint iY, Dynamics<T,Descriptor>* dynamics);
    /// Get a reference to the background dynamics
    Dynamics<T,Descriptor>& getBackgroundDynamics();
    /// Get a const reference to the background dynamics
    Dynamics<T,Descriptor> const& getBackgroundDynamics() const;
    /// Apply streaming step to bulk (non-boundary) cells
    void bulkStream(Box2D domain);
    /// Apply streaming step to boundary cells
    void boundaryStream(Box2D bound, Box2D domain);
    /// Apply collision and streaming step to bulk (non-boundary) cells
    void bulkCollideAndStream(Box2D domain);
private:
    /// Generic implementation of bulkCollideAndStream(domain).
    void linearBulkCollideAndStream(Box2D domain);
    /// Cache-efficient implementation of bulkCollideAndStream(domain)for
    ///   nearest-neighbor lattices.
    void blockwiseBulkCollideAndStream(Box2D domain);
private:
    /// Helper method for memory allocation
    void allocateAndInitialize();
    /// Helper method for memory de-allocation
    void releaseMemory();
    void implementPeriodicity();
private:
    void periodicDomain(Box2D domain);
private:
    Dynamics<T,Descriptor>* backgroundDynamics;
    /*
     * TODO: Inherit Cell from cudaManaged to enable Unified Memory;
     */
    Cell<T,Descriptor>     *rawData;
    Cell<T,Descriptor>    **grid;
    BlockLatticeDataTransfer2D<T,Descriptor> dataTransfer;
public:
    static CachePolicy2D& cachePolicy();
template<typename T_, template<typename U_> class Descriptor_>
    friend class ExternalRhoJcollideAndStream2D;

#ifndef PLB_CUDA_DISABLED
    /**
     * For CUDA_ENABLED, to enable \ref BlockLatticeDataTransfer2D access
     * \ref BlockLattice2D::rawData and \ref BlockLattice2D::grid directly;
     */
template<typename T_, template<typename U_> class Descriptor_>
    friend class BlockLatticeDataTransfer2D;
#endif

};

template<typename T, template<typename U> class Descriptor>
double getStoredAverageDensity(BlockLattice2D<T,Descriptor> const& blockLattice);

template<typename T, template<typename U> class Descriptor>
double getStoredAverageEnergy(BlockLattice2D<T,Descriptor> const& blockLattice);

template<typename T, template<typename U> class Descriptor>
double getStoredAverageVelocity(BlockLattice2D<T,Descriptor> const& blockLattice);

}  // namespace plb

#endif  // BLOCK_LATTICE_2D_H
