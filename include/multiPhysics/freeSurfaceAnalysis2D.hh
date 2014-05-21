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

#ifndef FREE_SURFACE_ANALYSIS_2D_HH
#define FREE_SURFACE_ANALYSIS_2D_HH

#include "multiPhysics/freeSurfaceAnalysis2D.h"
#include "multiPhysics/freeSurfaceModel2D.h"

namespace plb
{

template<typename T, template<typename U> class Descriptor>
FS_AverageMassFunctional2D<T,Descriptor>::FS_AverageMassFunctional2D()
    : averageMassId ( this->getStatistics().subscribeAverage() )
{ }

template<typename T, template<typename U> class Descriptor>
void FS_AverageMassFunctional2D<T,Descriptor>::processGenericBlocks (
    Box2D domain, std::vector<AtomicBlock2D*> atomicBlocks )
{
//     using namespace twoPhaseFlag;
    FreeSurfaceProcessorParam2D<T,Descriptor> param ( atomicBlocks );
    BlockStatistics& statistics = this->getStatistics();

    for ( plint iX=domain.x0; iX<=domain.x1; ++iX )
    {
        for ( plint iY=domain.y0; iY<=domain.y1; ++iY )
        {

            statistics.gatherAverage ( averageMassId, param.mass ( iX,iY) );
            statistics.incrementStats();
        }
    }
}

template<typename T, template<typename U> class Descriptor>
FS_AverageMassFunctional2D<T,Descriptor>* FS_AverageMassFunctional2D<T,Descriptor>::clone() const
{
    return new FS_AverageMassFunctional2D<T,Descriptor> ( *this );
}

template<typename T, template<typename U> class Descriptor>
T FS_AverageMassFunctional2D<T,Descriptor>::getAverageMass() const
{
    return this->getStatistics().getAverage ( averageMassId );
}

template<typename T, template<typename U> class Descriptor>
T freeSurfaceAverageMass ( std::vector<MultiBlock2D*> twoPhaseArgs, Box2D domain )
{
    FS_AverageMassFunctional2D<T,Descriptor> functional;
    applyProcessingFunctional ( functional, domain, twoPhaseArgs );
    return functional.getAverageMass();
}


template<typename T, template<typename U> class Descriptor>
FS_TotalMassFunctional2D<T,Descriptor>::FS_TotalMassFunctional2D()
    : totalMassId ( this->getStatistics().subscribeSum() )
{ }

template<typename T, template<typename U> class Descriptor>
void FS_TotalMassFunctional2D<T,Descriptor>::processGenericBlocks (
    Box2D domain, std::vector<AtomicBlock2D*> atomicBlocks )
{
    using namespace twoPhaseFlag;
    FreeSurfaceProcessorParam2D<T,Descriptor> param ( atomicBlocks );
    BlockStatistics& statistics = this->getStatistics();

    for ( plint iX=domain.x0; iX<=domain.x1; ++iX )
    {
        for ( plint iY=domain.y0; iY<=domain.y1; ++iY )
        {

            statistics.gatherSum ( totalMassId, param.mass ( iX,iY) );
            statistics.incrementStats();
        }
    }
}

template<typename T, template<typename U> class Descriptor>
FS_TotalMassFunctional2D<T,Descriptor>* FS_TotalMassFunctional2D<T,Descriptor>::clone() const
{
    return new FS_TotalMassFunctional2D<T,Descriptor> ( *this );
}

template<typename T, template<typename U> class Descriptor>
T FS_TotalMassFunctional2D<T,Descriptor>::getTotalMass() const
{
    return this->getStatistics().getSum ( totalMassId );
}


template<typename T, template<typename U> class Descriptor>
T freeSurfaceTotalMass ( std::vector<MultiBlock2D*> twoPhaseArgs, Box2D domain )
{
    FS_TotalMassFunctional2D<T,Descriptor> functional;
    applyProcessingFunctional ( functional, domain, twoPhaseArgs );
    return functional.getTotalMass();
}


template<typename T, template<typename U> class Descriptor>
FS_AverageDensityFunctional2D<T,Descriptor>::FS_AverageDensityFunctional2D()
    : averageDensityId ( this->getStatistics().subscribeAverage() )
{ }

template<typename T, template<typename U> class Descriptor>
void FS_AverageDensityFunctional2D<T,Descriptor>::processGenericBlocks (
    Box2D domain, std::vector<AtomicBlock2D*> atomicBlocks )
{
    using namespace twoPhaseFlag;
    FreeSurfaceProcessorParam2D<T,Descriptor> param ( atomicBlocks );
    BlockStatistics& statistics = this->getStatistics();
    BounceBack<T,Descriptor> BBdynamics;

    for ( plint iX=domain.x0; iX<=domain.x1; ++iX )
    {
        for ( plint iY=domain.y0; iY<=domain.y1; ++iY )
        {

            if ( param.cell ( iX,iY).getDynamics().getId() != BBdynamics.getId() )
            {
                statistics.gatherAverage ( averageDensityId, param.getDensity ( iX,iY) );
                statistics.incrementStats();
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
FS_AverageDensityFunctional2D<T,Descriptor>* FS_AverageDensityFunctional2D<T,Descriptor>::clone() const
{
    return new FS_AverageDensityFunctional2D<T,Descriptor> ( *this );
}

template<typename T, template<typename U> class Descriptor>
T FS_AverageDensityFunctional2D<T,Descriptor>::getAverageDensity() const
{
    return this->getStatistics().getAverage ( averageDensityId );
}


template<typename T, template<typename U> class Descriptor>
T freeSurfaceAverageDensity ( std::vector<MultiBlock2D*> twoPhaseArgs, Box2D domain )
{
    FS_AverageDensityFunctional2D<T,Descriptor> functional;
    applyProcessingFunctional ( functional, domain, twoPhaseArgs );
    return functional.getAverageDensity();
}


template<typename T, template<typename U> class Descriptor>
FS_AverageVolumeFractionFunctional2D<T,Descriptor>::FS_AverageVolumeFractionFunctional2D()
    : averageVfId ( this->getStatistics().subscribeAverage() )
{ }

template<typename T, template<typename U> class Descriptor>
void FS_AverageVolumeFractionFunctional2D<T,Descriptor>::processGenericBlocks (
    Box2D domain, std::vector<AtomicBlock2D*> atomicBlocks )
{
    using namespace twoPhaseFlag;
    FreeSurfaceProcessorParam2D<T,Descriptor> param ( atomicBlocks );
    BlockStatistics& statistics = this->getStatistics();

    BounceBack<T,Descriptor> BBdynamics;
    for ( plint iX=domain.x0; iX<=domain.x1; ++iX )
    {
        for ( plint iY=domain.y0; iY<=domain.y1; ++iY )
        {

            if ( param.cell ( iX,iY).getDynamics().getId() != BBdynamics.getId() )
            {
                statistics.gatherAverage ( averageVfId, param.volumeFraction ( iX,iY) );
                statistics.incrementStats();
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
FS_AverageVolumeFractionFunctional2D<T,Descriptor>* FS_AverageVolumeFractionFunctional2D<T,Descriptor>::clone() const
{
    return new FS_AverageVolumeFractionFunctional2D<T,Descriptor> ( *this );
}

template<typename T, template<typename U> class Descriptor>
T FS_AverageVolumeFractionFunctional2D<T,Descriptor>::getAverageVolumeFraction() const
{
    return this->getStatistics().getAverage ( averageVfId );
}


template<typename T, template<typename U> class Descriptor>
T freeSurfaceAverageVolumeFraction ( std::vector<MultiBlock2D*> twoPhaseArgs, Box2D domain )
{
    FS_AverageVolumeFractionFunctional2D<T,Descriptor> functional;
    applyProcessingFunctional ( functional, domain, twoPhaseArgs );
    return functional.getAverageVolumeFraction();
}



template<typename T, template<typename U> class Descriptor>
CountFreeSurfaceElementsFunctional2D<T,Descriptor>::CountFreeSurfaceElementsFunctional2D ( plint flagToLookFor_ )
    : numCellsId ( this->getStatistics().subscribeIntSum() ),
      flagToLookFor ( flagToLookFor_ )
{ }

template<typename T, template<typename U> class Descriptor>
void CountFreeSurfaceElementsFunctional2D<T,Descriptor>::processGenericBlocks (
    Box2D domain, std::vector<AtomicBlock2D*> atomicBlocks )
{
    using namespace twoPhaseFlag;
    FreeSurfaceProcessorParam2D<T,Descriptor> param ( atomicBlocks );
    BlockStatistics& statistics = this->getStatistics();

    for ( plint iX=domain.x0; iX<=domain.x1; ++iX )
    {
        for ( plint iY=domain.y0; iY<=domain.y1; ++iY )
        {

            int materialIndex = param.flag ( iX,iY);
            if ( materialIndex==flagToLookFor )  // Fluid Cell
            {
                statistics.gatherIntSum ( numCellsId, 1 );
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
CountFreeSurfaceElementsFunctional2D<T,Descriptor>* CountFreeSurfaceElementsFunctional2D<T,Descriptor>::clone() const
{
    return new CountFreeSurfaceElementsFunctional2D<T,Descriptor> ( *this );
}

template<typename T, template<typename U> class Descriptor>
plint CountFreeSurfaceElementsFunctional2D<T,Descriptor>::getNumInterfaceCells() const
{
    return this->getStatistics().getIntSum ( numCellsId );
}

template<typename T, template<typename U> class Descriptor>
plint countFreeSurfaceElements ( std::vector<MultiBlock2D*> twoPhaseArgs, plint flagToLookFor, Box2D domain )
{
    CountFreeSurfaceElementsFunctional2D<T,Descriptor> functional ( flagToLookFor );
    applyProcessingFunctional ( functional, domain, twoPhaseArgs );
    return functional.getNumInterfaceCells();
}


template<typename T, template<typename U> class Descriptor>
FS_AverageMomentumFunctional2D<T,Descriptor>::FS_AverageMomentumFunctional2D()
    : averageMomentumId (
        this->getStatistics().subscribeAverage(),
        this->getStatistics().subscribeAverage(),
        this->getStatistics().subscribeAverage() )
{ }

template<typename T, template<typename U> class Descriptor>
void FS_AverageMomentumFunctional2D<T,Descriptor>::processGenericBlocks (
    Box2D domain, std::vector<AtomicBlock2D*> atomicBlocks )
{
    using namespace twoPhaseFlag;
    FreeSurfaceProcessorParam2D<T,Descriptor> param ( atomicBlocks );
    BlockStatistics& statistics = this->getStatistics();
    BounceBack<T,Descriptor> BBdynamics;
    NoDynamics<T,Descriptor> NNdynamics;

    for ( plint iX=domain.x0; iX<=domain.x1; ++iX )
    {
        for ( plint iY=domain.y0; iY<=domain.y1; ++iY )
        {

            if ( param.cell ( iX,iY).getDynamics().getId() != BBdynamics.getId() &&
                    param.cell ( iX,iY).getDynamics().getId() != NNdynamics.getId() )
            {
                Array<T,Descriptor<T>::d> j = param.getMomentum ( iX,iY);
                statistics.gatherAverage ( averageMomentumId[0], j[0] );
                statistics.gatherAverage ( averageMomentumId[1], j[1] );
                statistics.gatherAverage ( averageMomentumId[2], j[2] );
                statistics.incrementStats();
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
FS_AverageMomentumFunctional2D<T,Descriptor>* FS_AverageMomentumFunctional2D<T,Descriptor>::clone() const
{
    return new FS_AverageMomentumFunctional2D<T,Descriptor> ( *this );
}

template<typename T, template<typename U> class Descriptor>
Array<T,3> FS_AverageMomentumFunctional2D<T,Descriptor>::getAverageMomentum() const
{
    return Array<T,3> (
               this->getStatistics().getAverage ( averageMomentumId[0] ),
               this->getStatistics().getAverage ( averageMomentumId[1] ),
               this->getStatistics().getAverage ( averageMomentumId[2] ) );
}

template<typename T, template<typename U> class Descriptor>
Array<T,3> freeSurfaceAverageMomentum ( std::vector<MultiBlock2D*> twoPhaseArgs, Box2D domain )
{
    FS_AverageMomentumFunctional2D<T,Descriptor> functional;
    applyProcessingFunctional ( functional, domain, twoPhaseArgs );
    return functional.getAverageMomentum();
}


template<typename T, template<typename U> class Descriptor>
FS_AverageHeightFunctional2D<T,Descriptor>::FS_AverageHeightFunctional2D()
    : averageHeightId ( this->getStatistics().subscribeAverage() )
{ }

template<typename T, template<typename U> class Descriptor>
void FS_AverageHeightFunctional2D<T,Descriptor>::processGenericBlocks (
    Box2D domain, std::vector<AtomicBlock2D*> atomicBlocks )
{
    using namespace twoPhaseFlag;
    FreeSurfaceProcessorParam2D<T,Descriptor> param ( atomicBlocks );
    BlockStatistics& statistics = this->getStatistics();
    BounceBack<T,Descriptor> BBdynamics;
    Dot2D absOffset = param.absOffset();

    for ( plint iX=domain.x0; iX<=domain.x1; ++iX )
    {
        for ( plint iY=domain.y0; iY<=domain.y1; ++iY )
        {

            if ( param.cell ( iX,iY).getDynamics().getId() != BBdynamics.getId() )
            {
                T localHeight = T ( 0 );
                if ( param.volumeFraction ( iX,iY) == 1 )
                {
					//TODO FIX THIS
//                     localHeight = T ( absOffset.z );
                }
                statistics.gatherAverage ( averageHeightId, localHeight );
                statistics.incrementStats();
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
FS_AverageHeightFunctional2D<T,Descriptor>* FS_AverageHeightFunctional2D<T,Descriptor>::clone() const
{
    return new FS_AverageHeightFunctional2D<T,Descriptor> ( *this );
}

template<typename T, template<typename U> class Descriptor>
T FS_AverageHeightFunctional2D<T,Descriptor>::getAverageHeight() const
{
    return this->getStatistics().getAverage ( averageHeightId );
}

template<typename T, template<typename U> class Descriptor>
T freeSurfaceAverageHeight ( std::vector<MultiBlock2D*> twoPhaseArgs, Box2D domain )
{
    FS_AverageHeightFunctional2D<T,Descriptor> functional;
    applyProcessingFunctional ( functional, domain, twoPhaseArgs );
    return functional.getAverageHeight();
}


template<typename T, template<typename U> class Descriptor>
GetWaterLevelAtxyFunctional2D<T,Descriptor>::GetWaterLevelAtxyFunctional2D()
    : numFluidOccupiedCellId ( this->getStatistics().subscribeIntSum() )
{ }

template<typename T, template<typename U> class Descriptor>
void GetWaterLevelAtxyFunctional2D<T,Descriptor>::processGenericBlocks (
    Box2D domain, std::vector<AtomicBlock2D*> atomicBlocks )
{
    using namespace twoPhaseFlag;
    FreeSurfaceProcessorParam2D<T,Descriptor> param ( atomicBlocks );
    BlockStatistics& statistics = this->getStatistics();
    plint bbDynamics = BounceBack<T,Descriptor>().getId();

    for ( plint iX=domain.x0; iX<=domain.x1; ++iX )
    {
        for ( plint iY=domain.y0; iY<=domain.y1; ++iY )
        {

            if ( param.cell ( iX,iY).getDynamics().getId() != bbDynamics )
            {
                if ( param.volumeFraction ( iX,iY) >= 0.5 )
                {
                    statistics.gatherIntSum ( numFluidOccupiedCellId, 1 );
                }
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
GetWaterLevelAtxyFunctional2D<T,Descriptor>* GetWaterLevelAtxyFunctional2D<T,Descriptor>::clone() const
{
    return new GetWaterLevelAtxyFunctional2D<T,Descriptor> ( *this );
}

template<typename T, template<typename U> class Descriptor>
plint GetWaterLevelAtxyFunctional2D<T,Descriptor>::getNumFluidCellsAtXY() const
{
    return this->getStatistics().getIntSum ( numFluidOccupiedCellId );
}

template<typename T, template<typename U> class Descriptor>
T getAverageHeightAtXY ( std::vector<MultiBlock2D*> twoPhaseArgs, plint N, Box2D domain )
{
    GetWaterLevelAtxyFunctional2D<T,Descriptor> functional;
    applyProcessingFunctional ( functional, domain, twoPhaseArgs );
    plint length_domain = domain.x1-domain.x0 ; // number of cell along y direction
    if ( length_domain==0 )
        length_domain =1;
    plint width_domain = domain.y1-domain.y0 ; // number of cell along y direction
    if ( width_domain==0 )
        width_domain =1;
    T heightAtXY = functional.getNumFluidCellsAtXY() / ( T ( N ) *length_domain*width_domain );
    return heightAtXY;
}

}  // namespace plb

#endif  // FREE_SURFACE_ANALYSIS_2D_HH

