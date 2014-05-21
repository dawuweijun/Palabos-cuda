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

#ifndef FREE_SURFACE_ANALYSIS_2D_H
#define FREE_SURFACE_ANALYSIS_2D_H

#include "multiBlock/multiBlockLattice2D.h"

namespace plb {

template<typename T, template<typename U> class Descriptor>
T freeSurfaceAverageMass(std::vector<MultiBlock2D*> twoPhaseArgs, Box2D domain);
             
template<typename T, template<typename U> class Descriptor>
class FS_AverageMassFunctional2D : public PlainReductiveBoxProcessingFunctional2D
{
public:
    FS_AverageMassFunctional2D();
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D*> atomicBlocks);
    virtual FS_AverageMassFunctional2D<T,Descriptor>* clone() const;
    T getAverageMass() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        for (pluint i=0; i<modified.size(); ++i) {
            modified[i] = modif::nothing;
        }
    }
private:
    plint averageMassId;
};


template<typename T, template<typename U> class Descriptor>
T freeSurfaceTotalMass(std::vector<MultiBlock2D*> twoPhaseArgs, Box2D domain);
             
template<typename T, template<typename U> class Descriptor>
class FS_TotalMassFunctional2D : public PlainReductiveBoxProcessingFunctional2D
{
public:
    FS_TotalMassFunctional2D();
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D*> atomicBlocks);
    virtual FS_TotalMassFunctional2D<T,Descriptor>* clone() const;
    T getTotalMass() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        for (pluint i=0; i<modified.size(); ++i) {
            modified[i] = modif::nothing;
        }
    }
private:
    plint totalMassId;
};


template<typename T, template<typename U> class Descriptor>
T freeSurfaceAverageDensity(std::vector<MultiBlock2D*> twoPhaseArgs, Box2D domain);
             
template<typename T, template<typename U> class Descriptor>
class FS_AverageDensityFunctional2D : public PlainReductiveBoxProcessingFunctional2D
{
public:
    FS_AverageDensityFunctional2D();
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D*> atomicBlocks);
    virtual FS_AverageDensityFunctional2D<T,Descriptor>* clone() const;
    T getAverageDensity() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        for (pluint i=0; i<modified.size(); ++i) {
            modified[i] = modif::nothing;
        }
    }
private:
    plint averageDensityId;
};


template<typename T, template<typename U> class Descriptor>
T freeSurfaceAverageVolumeFraction(std::vector<MultiBlock2D*> twoPhaseArgs, Box2D domain);
             
template<typename T, template<typename U> class Descriptor>
class FS_AverageVolumeFractionFunctional2D : public PlainReductiveBoxProcessingFunctional2D
{
public:
    FS_AverageVolumeFractionFunctional2D();
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D*> atomicBlocks);
    virtual FS_AverageVolumeFractionFunctional2D<T,Descriptor>* clone() const;
    T getAverageVolumeFraction() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        for (pluint i=0; i<modified.size(); ++i) {
            modified[i] = modif::nothing;
        }
    }
private:
    plint averageVfId;
};


template<typename T, template<typename U> class Descriptor>
plint countFreeSurfaceElements(std::vector<MultiBlock2D*> twoPhaseArgs, plint flagToLookFor, Box2D domain);

template<typename T, template<typename U> class Descriptor>
class CountFreeSurfaceElementsFunctional2D : public PlainReductiveBoxProcessingFunctional2D
{
public:
    CountFreeSurfaceElementsFunctional2D(plint flagToLookFor_);
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D*> atomicBlocks);
    virtual CountFreeSurfaceElementsFunctional2D<T,Descriptor>* clone() const;
    plint getNumInterfaceCells() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        for (pluint i=0; i<modified.size(); ++i) {
            modified[i] = modif::nothing;
        }
    }
private:
    plint numCellsId;
    plint flagToLookFor;
};
    

template<typename T, template<typename U> class Descriptor>
Array<T,3> freeSurfaceAverageMomentum(std::vector<MultiBlock2D*> twoPhaseArgs, Box2D domain);
             
template<typename T, template<typename U> class Descriptor>
class FS_AverageMomentumFunctional2D : public PlainReductiveBoxProcessingFunctional2D
{
public:
    FS_AverageMomentumFunctional2D();
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D*> atomicBlocks);
    virtual FS_AverageMomentumFunctional2D<T,Descriptor>* clone() const;
    Array<T,3> getAverageMomentum() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        for (pluint i=0; i<modified.size(); ++i) {
            modified[i] = modif::nothing;
        }
    }
private:
    Array<plint,3> averageMomentumId;
};


template<typename T, template<typename U> class Descriptor>
T freeSurfaceAverageHeight(std::vector<MultiBlock2D*> twoPhaseArgs, Box2D domain);
             
template<typename T, template<typename U> class Descriptor>
class FS_AverageHeightFunctional2D : public PlainReductiveBoxProcessingFunctional2D
{
public:
    FS_AverageHeightFunctional2D();
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D*> atomicBlocks);
    virtual FS_AverageHeightFunctional2D<T,Descriptor>* clone() const;
    T getAverageHeight() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        for (pluint i=0; i<modified.size(); ++i) {
            modified[i] = modif::nothing;
        }
    }
private:
    plint averageHeightId;
};


template<typename T, template<typename U> class Descriptor>
T getAverageHeightAtXY(std::vector<MultiBlock2D*> twoPhaseArgs, plint N, Box2D domain);

template<typename T, template<typename U> class Descriptor>
class GetWaterLevelAtxyFunctional2D : public PlainReductiveBoxProcessingFunctional2D
{
public:
    GetWaterLevelAtxyFunctional2D();
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D*> atomicBlocks);
    virtual GetWaterLevelAtxyFunctional2D<T,Descriptor>* clone() const;
    plint getNumFluidCellsAtXY() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        for (pluint i=0; i<modified.size(); ++i) {
            modified[i] = modif::nothing;
        }
    }
private:
    plint numFluidOccupiedCellId;
};

}  // namespace plb

#endif  // FREE_SURFACE_ANALYSIS_2D_H

