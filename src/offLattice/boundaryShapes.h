#ifndef BOUNDARY_SHAPES_H
#define BOUNDARY_SHAPES_H
namespace plb{
struct OffBoundary {
    enum Type {dirichlet, neumann, freeSlip, constRhoInlet, densityNeumann, flux, isolation};
};
}
#endif