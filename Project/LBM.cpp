#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <vector>
#include "LBM.h"

using namespace std;

void init_equilibrium(vector<double> &f0, vector<double> &f1, vector<double> &r,
                      vector<double> &u, vector<double> &v)
{
    for (unsigned int y = 0; y < NY; ++y)
    {
        for (unsigned int x = 0; x < NX; ++x)
        {
            double rho = r[scalar_index(x, y)];
            double ux = u[scalar_index(x, y)];
            double uy = v[scalar_index(x, y)];
            f0[field0_index(x,y)] = w0*(1.0 - 1.5*(ux*ux+uy*uy));
            for (unsigned int i = 0; i < ndir; ++i)
            {
                double cidotu = dirx[i] * ux + diry[i] * uy;
                f1[field_index(x, y, i)] =
                    wi[i] * rho * (1.0 + 3.0 * cidotu + 4.5 * cidotu * cidotu - 1.5 * (ux * ux + uy * uy));
            }
        }
    }
}

void stream(vector<double> &f_src, vector<double> &f_dst)
{
    for (unsigned int y = 0; y < NY; ++y)
    {
        for (unsigned int x = 0; x < NX; ++x)
        {
            for (unsigned int i = 0; i < ndir; ++i)
            {
                // enforce periodicity
                // add NX to ensure that value is positive
                unsigned int xmd = (NX + x - dirx[i]) % NX;
                unsigned int ymd = (NY + y - diry[i]) % NY;
                f_dst[field_index(x, y, i)] = f_src[field_index(xmd, ymd, i)];
            }
        }
    }
}

void compute_rho_u(vector<double> &f, vector<double> &r,
                   vector<double> &u, vector<double> &v)
{
    for (unsigned int y = 0; y < NY; ++y)
    {
        for (unsigned int x = 0; x < NX; ++x)
        {
            double rho = 0.0;
            double ux = 0.0;
            double uy = 0.0;
            for (unsigned int i = 0; i < ndir; ++i)
            {
                rho += f[field_index(x, y, i)];
                ux += dirx[i] * f[field_index(x, y, i)];
                uy += diry[i] * f[field_index(x, y, i)];
            }
            r[scalar_index(x, y)] = rho;
            u[scalar_index(x, y)] = ux / rho;
            v[scalar_index(x, y)] = uy / rho;
        }
    }
}

void collide(vector<double> &f, vector<double> &r, vector<double> &u, vector<double> &v)
{
    // useful constants
    const double tauinv = 2.0 / (6.0 * nu + 1.0); // 1/tau
    const double omtauinv = 1.0 - tauinv;
    // 1 - 1/tau
    for (unsigned int y = 0; y < NY; ++y)
    {
        for (unsigned int x = 0; x < NX; ++x)
        {
            double rho = r[scalar_index(x, y)];
            double ux = u[scalar_index(x, y)];
            double uy = v[scalar_index(x, y)];
            for (unsigned int i = 0; i < ndir; ++i)
            {
                // calculate dot product
                double cidotu = dirx[i] * ux + diry[i] * uy;
                // calculate equilibrium
                double feq = wi[i] * rho * (1.0 + 3.0 * cidotu + 4.5 * cidotu * cidotu - 1.5 * (ux * ux + uy * uy));
                // relax to equilibrium
                f[field_index(x, y, i)] =
                    omtauinv * f[field_index(x, y, i)] + tauinv * feq;
            }
        }
    }
}

// Function that performs streaming, computation of moments, and collision in one step
void stream_collide_save(vector<double> &f0, vector<double> &f1, vector<double> &f2,
                         vector<double> &r, vector<double> &u, vector<double> &v,
                         bool save)
{
    // useful constants
    const double tauinv = 2.0 / (6.0 * nu + 1.0); // 1/tau
    const double omtauinv = 1.0 - tauinv;
    // 1 - 1/tau
    for (unsigned int y = 0; y < NY; ++y)
    {
        for (unsigned int x = 0; x < NX; ++x)
        {
            unsigned int xp1 = (x + 1) % NX;
            unsigned int yp1 = (y + 1) % NY;
            unsigned int xm1 = (NX + x - 1) % NX;
            unsigned int ym1 = (NY + y - 1) % NY;
            // direction numbering scheme
            // 6 2 5
            // 3 0 1
            // 7 4 8
            double ft0 = f0[field0_index(x, y)];
            // load populations from adjacent nodes
            double ft1 = f1[fieldn_index(xm1, y, 1)];
            double ft2 = f1[fieldn_index(x, ym1, 2)];
            double ft3 = f1[fieldn_index(xp1, y, 3)];
            double ft4 = f1[fieldn_index(x, yp1, 4)];
            double ft5 = f1[fieldn_index(xm1, ym1, 5)];
            double ft6 = f1[fieldn_index(xp1, ym1, 6)];
            double ft7 = f1[fieldn_index(xp1, yp1, 7)];
            double ft8 = f1[fieldn_index(xm1, yp1, 8)];
            // compute moments
            double rho = ft0 + ft1 + ft2 + ft3 + ft4 + ft5 + ft6 + ft7 + ft8;
            double rhoinv = 1.0 / rho;
            double ux = rhoinv * (ft1 + ft5 + ft8 - (ft3 + ft6 + ft7));
            double uy = rhoinv * (ft2 + ft5 + ft6 - (ft4 + ft7 + ft8));
            // only write to memory when needed
            if (save)
            {
                r[scalar_index(x, y)] = rho;
                u[scalar_index(x, y)] = ux;
                v[scalar_index(x, y)] = uy;
            }
            // now compute and relax to equilibrium
            // note that
            // feq_i = w_i rho [1 + 3(ci . u)
            // +(9/2) (ci . u)^2 - (3/2) (u.u)]
            //= w_i rho [1 - 3/2 (u.u)
            // +(ci . 3u) + (1/2) (ci . 3u)^2]
            // = w_i rho [1 - 3/2 (u.u)
            // +(ci . 3u)(1 + (1/2) (ci . 3u))]
            // temporary variables
            double tw0r = tauinv * w0 * rho;                // w[0]*rho/tau
            double twsr = tauinv * wst * rho;                // w[1-4]*rho/tau
            double twdr = tauinv * wd * rho;                // w[5-8]*rho/tau
            double omusq = 1.0 - 1.5 * (ux * ux + uy * uy); // 1-(3/2)u.u
            double tux = 3.0 * ux;
            double tuy = 3.0 * uy;
            f0[field0_index(x, y)] = omtauinv * ft0 + tw0r * (omusq);
            double cidot3u = tux;
            f2[fieldn_index(x, y, 1)] = omtauinv * ft1 + twsr * (omusq + cidot3u * (1.0 + 0.5 * cidot3u));
            // ... similar expressions for directions 2-4
            cidot3u = tux + tuy;
            f2[fieldn_index(x, y, 5)] = omtauinv * ft5 + twdr * (omusq + cidot3u * (1.0 + 0.5 * cidot3u));
            // ... similar expressions for directions 6-8
        }
    }
}

void compute_flow_properties(unsigned int t, vector<double> r, vector<double> u, vector<double> v, vector<double> prop)
{

    // prop must point to space for 4 doubles:
    // 0: energy
    // 1: L2 error in rho
    // 2: L2 error in ux
    // 3: L2 error in uy
    double E = 0.0;
    double sumrhoe2 = 0.0;
    double sumuxe2 = 0.0;
    double sumuye2 = 0.0;
    double sumrhoa2 = 0.0;
    double sumuxa2 = 0.0;
    double sumuya2 = 0.0;
    for (unsigned int y = 0; y < NY; ++y)
    {
        for (unsigned int x = 0; x < NX; ++x)
        {
            double rho = r[scalar_index(x, y)];
            double ux = u[scalar_index(x, y)];
            double uy = v[scalar_index(x, y)];
            E += rho * (ux * ux + uy * uy);
            double rhoa, uxa, uya;
            // taylor_green(t,x,y,&rhoa,&uxa,&uya); //maybe for the future?
            sumrhoe2 += (rho - rhoa) * (rho - rhoa);
            sumuxe2 += (ux - uxa) * (ux - uxa);
            sumuye2 += (uy - uya) * (uy - uya);
            sumrhoa2 += (rhoa - rho0) * (rhoa - rho0);
            sumuxa2 += uxa * uxa;
            sumuya2 += uya * uya;
        }
    }
    prop[0] = E;
    prop[1] = sqrt(sumrhoe2 / sumrhoa2);
    prop[2] = sqrt(sumuxe2 / sumuxa2);
    prop[3] = sqrt(sumuye2 / sumuya2);
}

void save_scalar(const char *name, vector<double> &scalar,
                 unsigned int n)
{
    // assume reasonably-sized file names
    char filename[128];
    char format[16];
    // compute maximum number of digits
    int ndigits = floor(log10((double)NSTEPS) + 1.0);
    // generate format string
    // file name format is name0000nnn.bin
    sprintf(format, "%%s%%0%dd.bin", ndigits);
    sprintf(filename, format, name, n);
    // open file for writing
    FILE *fout = fopen(filename, "wb+");
    // write data
    fwrite(&scalar, 1, mem_size_scalar, fout);
    // close file
    fclose(fout);
}