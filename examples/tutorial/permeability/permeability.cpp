/* This file is part of the Palabos library.
 *
 * Copyright (C) 2011-2012 FlowKit Sarl
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

/* Main author: Wim Degruyter */

#include "palabos3D.h"
#include "palabos3D.hh"
#include <vector>
#include <cmath>

using namespace plb;
using namespace std;

typedef double T;
#define DESCRIPTOR descriptors::D3Q19Descriptor

// This function object returns a zero velocity, and a pressure which decreases
//   linearly in x-direction. It is used to initialize the particle populations.
class PressureGradient {
public:
        PressureGradient(T deltaP_, plint nx_) : deltaP(deltaP_), nx(nx_)
        { }
        void operator() (plint iX, plint iY, plint iZ, T& density, Array<T,3>& velocity) const
        {
                velocity.resetToZero();
                density = 1. - deltaP*DESCRIPTOR<T>::invCs2 / (T)(nx-1) * (T)iX;

        }
private:
        T deltaP;
        plint nx;
};

void porousMediaSetup( MultiBlockLattice3D<T,DESCRIPTOR>& lattice,
                       OnLatticeBoundaryCondition3D<T,DESCRIPTOR>* boundaryCondition,
                       MultiScalarField3D<int>& geometry, T deltaP)
{
        const plint nx = lattice.getNx();
        const plint ny = lattice.getNy();
        const plint nz = lattice.getNz();

        pcout << "Definition of inlet/outlet." << endl;
        Box3D inlet (0,0,      1,ny-2, 1,nz-2);
        boundaryCondition->addPressureBoundary0N(inlet, lattice);
        setBoundaryDensity(lattice, inlet, 1.);

        Box3D outlet(nx-1,nx-1, 1,ny-2, 1,nz-2);
        boundaryCondition->addPressureBoundary0P(outlet, lattice);
        setBoundaryDensity(lattice, outlet, 1. - deltaP*DESCRIPTOR<T>::invCs2);

        pcout << "Definition of the geometry." << endl;
        // Where "geometry" evaluates to 1, use bounce-back.
        defineDynamics(lattice, geometry, new BounceBack<T,DESCRIPTOR>(), 1);
        // Where "geometry" evaluates to 2, use no-dynamics (which does nothing).
        defineDynamics(lattice, geometry, new NoDynamics<T,DESCRIPTOR>(), 2);

        pcout << "Initilization of rho and u." << endl;
        initializeAtEquilibrium( lattice, lattice.getBoundingBox(), PressureGradient(deltaP, nx) );

        lattice.initialize();
        delete boundaryCondition;
}

void writeGifs(MultiBlockLattice3D<T,DESCRIPTOR>& lattice, plint iter)
{
        const plint nx = lattice.getNx();
        const plint ny = lattice.getNy();
        const plint nz = lattice.getNz();

        const plint imSize = 600;
        ImageWriter<T> imageWriter("leeloo");

        // Write velocity-norm at x=0.
        imageWriter.writeScaledGif( createFileName("ux_inlet", iter, 6),
                        *computeVelocityNorm(lattice, Box3D(0,0, 0,ny-1, 0,nz-1)),
                        imSize, imSize );

        // Write velocity-norm at x=nx/2.
        imageWriter.writeScaledGif( createFileName("ux_half", iter, 6),
                        *computeVelocityNorm(lattice, Box3D(nx/2,nx/2, 0,ny-1, 0,nz-1)),
                        imSize, imSize );

}


void writeVTK(MultiBlockLattice3D<T,DESCRIPTOR>& lattice, plint iter)
{
        VtkImageOutput3D<T> vtkOut(createFileName("vtk", iter, 6), 1.);
        vtkOut.writeData<float>(*computeVelocityNorm(lattice), "velocityNorm", 1.);
        vtkOut.writeData<3,float>(*computeVelocity(lattice), "velocity", 1.);
}

T computePermeability (
                MultiBlockLattice3D<T,DESCRIPTOR>& lattice, T nu, T deltaP, Box3D domain )
{
        pcout << "Computing the permeability." << endl;

        // Compute only the x-direction of the velocity (direction of the flow).
        plint xComponent = 0;
        plint nx = lattice.getNx();

        T meanU = computeAverage (
                                *computeVelocityComponent (lattice, domain, xComponent ) );


        pcout << "Average velocity     = " << meanU             << endl;
        pcout << "Lattice viscosity nu = " << nu                << endl;
        pcout << "Grad P               = " << deltaP/(T)(nx-1)            << endl;
        pcout << "Permeability         = " << nu*meanU / (deltaP/(T)(nx-1)) << endl;

        return meanU;
}

int main(int argc, char **argv)
{
        plbInit(&argc, &argv);

        if (argc!=7)
        {
                pcout << "Error missing some input parameter\n";
                pcout << "The structure is :\n";
                pcout << "1. Input file name.\n";
                pcout << "2. Output directory name.\n";
                pcout << "3. number of cells in X direction.\n";
                pcout << "4. number of cells in Y direction.\n";
                pcout << "5. number of cells in Z direction.\n";
                pcout << "6. Delta P .\n";
                pcout << "Example: " << argv[0] << " twoSpheres.dat tmp/ 48 64 64 0.00005\n";
                exit (EXIT_FAILURE);
        }
        std::string fNameIn  = argv[1];
        std::string fNameOut = argv[2];

        const plint nx = atoi(argv[3]);
        const plint ny = atoi(argv[4]);
        const plint nz = atoi(argv[5]);
        const T   deltaP = atof(argv[6]);

        global::directories().setOutputDir(fNameOut+"/");

        const T omega = 1.0;
        const T nu    = ((T)1/omega-0.5)/DESCRIPTOR<T>::invCs2;

        pcout << "Creation of the lattice." << endl;
        MultiBlockLattice3D<T,DESCRIPTOR> lattice(nx,ny,nz, new BGKdynamics<T,DESCRIPTOR>(omega));
        // Switch off periodicity.
        lattice.periodicity().toggleAll(false);

        ///////
        
        MultiScalarField3D<int> geometry(nx,ny,nz);
        MultiScalarField3D<int> slice(1,ny,nz);
        for (iX=0; iX<nx-1; ++iX) {
            string fname = createFileName("slice_", iX, 4)+"_truc.dat";
            pcout << "Reading slice " << fname;
            plb_ifstream geometryFile(fNameIn.c_str());
            if (!geometryFile.is_open()) {
                    pcout << "Error: could not open geometry file " << fNameIn << endl;
                    return -1;
            }
            geometryFile >> slice;
            copy(slice, slice.getBoundingBox(), geometry, Box3D(iX,iX,0,ny-1,0,nz-1));
        }
       
        ////////

        pcout << "Reading the geometry file." << endl;
        MultiScalarField3D<int> geometry(nx,ny,nz);
        plb_ifstream geometryFile(fNameIn.c_str());
        if (!geometryFile.is_open()) {
                pcout << "Error: could not open geometry file " << fNameIn << endl;
                return -1;
        }
        geometryFile >> geometry;

        pcout << "nu = " << nu << endl;
        pcout << "deltaP = " << deltaP << endl;
        pcout << "omega = " << omega << endl;
        pcout << "nx = " << lattice.getNx() << endl;
        pcout << "ny = " << lattice.getNy() << endl;
        pcout << "nz = " << lattice.getNz() << endl;

        porousMediaSetup(lattice, createLocalBoundaryCondition3D<T,DESCRIPTOR>(), geometry, deltaP);

        // The value-tracer is used to stop the simulation once is has converged.
        // 1st parameter:velocity
        // 2nd parameter:size
        // 3rd parameters:threshold
        // 1st and second parameters ae used for the length of the time average (size/velocity)
        util::ValueTracer<T> converge(1.0,1000.0,1.0e-4);

        pcout << "Simulation begins" << endl;
        plint iT=0;

        const plint maxT = 30000;
        for (;iT<maxT; ++iT)
        {
                if (iT % 20 == 0) {
                        pcout << "Iteration " << iT << endl;
                }
                if (iT % 500 == 0 && iT>0) {
                        writeGifs(lattice,iT);
                        }

                lattice.collideAndStream();
                converge.takeValue(getStoredAverageEnergy(lattice),true);

                if (converge.hasConverged())
                {
                        break;
                }
        }

        pcout << "End of simulation at iteration " << iT << endl;

        pcout << "Permeability:" << endl << endl;
        computePermeability(lattice, nu, deltaP, lattice.getBoundingBox());
        pcout << endl;

        pcout << "Writing VTK file ..." << endl << endl;
        writeVTK(lattice,iT);
        pcout << "Finished!" << endl << endl;
}
