#ifndef __LBM_H
#define __LBM_H

#include <vector>
#include <cmath>
#include <string>

using namespace std;

void init_equilibrium(std::vector<double>&, std::vector<double>&, std::vector<double>&, std::vector<double>&, std::vector<double>&);
void compute_rho_u(std::vector<double>&, std::vector<double>&, std::vector<double>&, std::vector<double>&);
void collide(std::vector<double>&, std::vector<double>&, std::vector<double>&, std::vector<double>&);
void stream_collide_save(std::vector<double>&, std::vector<double>&, std::vector<double>&, std::vector<double>&, std::vector<double>&, std::vector<double>&, bool);
void compute_flow_properties(unsigned int, std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>);
void save_scalar(const char*, std::vector<double>&, unsigned int);

const unsigned int scale = 1;
const unsigned int NX = 32 * scale;
const unsigned int NY = NX;
const unsigned int ndir = 9;
const size_t size0 = NX * NY;
const size_t size1 = NX * NY * (ndir-1);
const size_t mem_size_0dir   = sizeof(double)*NX*NY;
const size_t mem_size_ndir = sizeof(double) * NX * NY * (ndir-1);
const size_t mem_size_scalar = sizeof(double) * NX * NY;
const size_t total_mem_bytes = mem_size_0dir + 2*mem_size_ndir + 3*mem_size_scalar;

const unsigned int NSTEPS = 40*scale*scale;
const unsigned int NSAVE  =  20*scale*scale;
const unsigned int NMSG   =  20*scale*scale;
const bool computeFlowProperties = true;
const bool quiet = true;

const double w0 = 4.0 / 9.0;  // zero weight
const double wst = 1.0 / 9.0;  // adjacent weight
const double wd = 1.0 / 36.0; // diagonal weight

const double cs = 1.0/sqrt(3.0); // sound velocity


const double vl[] = {0.0, 0.0}; //speed of left wall
const double vr[] = {0.0, 0.0}; // " " right wall
const double vu[] = {0.04, 0.0}; // " " upper wall
const double vd[] = {0.0, 0.0}; // " " lower wall




// Arrays of the lattice weights and direction components
const double wi[] = {w0, wst, wst, wst, wst, wd, wd, wd, wd};
const int dirx[] = {0, 1, 0, -1, 0, 1, -1, -1, 1};
const int diry[] = {0, 0, 1, 0, -1, 1, 1, -1, -1};

const int index_opp[] = {0, 3, 4, 1, 2, 7, 8, 5, 6}; // index of opposite directions


// The kinematic viscosity  and the corresponding relaxation parameter
const double nu = 1.0 / 6.0;
const double tau = 3.0 * nu + 0.5;

// The maximum flow speed
const double u_max = 0.04 / scale;

// The fluid density
const double rho0 = 1.0;

// The number of time steps in the simulation

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