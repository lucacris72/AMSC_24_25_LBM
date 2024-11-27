#ifndef __LBM_H
#define __LBM_H

#include <vector>

using namespace std;

const unsigned int scale = 1;
const unsigned int NX = 32 * scale;
const unsigned int NY = NX;
const unsigned int ndir = 9;
const size_t size0 = NX * NY;
const size_t size1 = NX * NY * ndir;
const size_t mem_size_0dir   = sizeof(double)*NX*NY;
const size_t mem_size_ndir = sizeof(double) * NX * NY * ndir;
const size_t mem_size_scalar = sizeof(double) * NX * NY;
const size_t total_mem_bytes = mem_size_0dir + 2*mem_size_ndir + 3*mem_size_scalar;

const unsigned int NSTEPS = 200*scale*scale;
const unsigned int NSAVE  =  50*scale*scale;
const unsigned int NMSG   =  50*scale*scale;
const bool computeFlowProperties = true;
const bool quiet = true;

const double w0 = 4.0 / 9.0;  // zero weight
const double ws = 1.0 / 9.0;  // adjacent weight
const double wd = 1.0 / 36.0; // diagonal weight

// Arrays of the lattice weights and direction components
const double wi[] = {w0, ws, ws, ws, ws, wd, wd, wd, wd};
const int dirx[] = {0, 1, 0, -1, 0, 1, -1, -1, 1};
const int diry[] = {0, 0, 1, 0, -1, 1, 1, -1, -1};

// The kinematic viscosity  and the corresponding relaxation parameter
const double nu = 1.0 / 6.0;
const double tau = 3.0 * nu + 0.5;

// The maximum flow speed
const double u_max = 0.04 / scale;

// The fluid density
const double rho0 = 1.0;

// The number of time steps in the simulation
const unsigned int NSTEPS = 200 * scale * scale;

inline size_t scalar_index(unsigned int x, unsigned int y)
{
    return NX * y + x;
}

inline size_t field0_index(unsigned int x, unsigned int y)
{
    return NX * y + x;
}
inline size_t fieldn_index(unsigned int x, int y, unsigned int d)
{
    return (ndir - 1) * (NX * (y + 1) + x) + (d - 1);
}

inline size_t field_index(unsigned int x, unsigned int y, unsigned int d)
{
    return NX * (NY * d + y) + x;
}

#endif