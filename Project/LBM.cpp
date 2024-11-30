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
                f1[fieldn_index(x, y, i)] =
                    wi[i] * rho * (1.0 + 3.0 * cidotu + 4.5 * cidotu * cidotu - 1.5 * (ux * ux + uy * uy));
            }
        }
    }
}

/*void stream(vector<double> &f_src, vector<double> &f_dst)
{
    for (unsigned int y = 0; y < NY; ++y)
    {
        for (unsigned int x = 0; x < NX; ++x)
        {
            for (unsigned int i = 0; i < ndir; ++i)
            {
                // enforce periodicity
                // add NX to ensure that value is positive
                unsigned int xmd = (x - dirx[i]);
                unsigned int ymd = (y - diry[i]);

                if (!(0 <= xmd < NX || 0 <= ymd < NY))
                {
                    corner (f_src, f_dst, xmd, ymd);
                }

                else if (!(0 <= xmd < NX && 0 <= ymd < NY))
                {
                    bounceback (f_src, f_dst, xmd, ymd);
                }
                else
                {
                f_dst[field_index(x, y, i)] = f_src[field_index(xmd, ymd, i)];
                }            
            }
        }
    }
}*/

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
                rho += f[fieldn_index(x, y, i)];
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
                f[fieldn_index(x, y, i)] =
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
            double ft[ndir];        
                       
            /*unsigned int xp1 = (x + 1);
            unsigned int yp1 = (y + 1);
            unsigned int xm1 = (x - 1);
            unsigned int ym1 = (y - 1);*/

            // direction numbering scheme
            // 6 2 5
            // 3 0 1
            // 7 4 8

            ft[0] = f0[field0_index(x, y)];

            for (unsigned int i = 1; i < ndir; ++i)
            {
                unsigned int xmd = (x - dirx[i]) ;
                unsigned int ymd = (y - diry[i]) ;

                bool right = 0;
                bool left = 0;
                bool up = 0;
                bool down = 0;

                
                if (xmd < 0)
                {
                    left = 1;
                }

                if (xmd == NX)
                {
                    right = 1;
                }

                if (ymd == -1)
                {
                    up = 1;
                }

                if (ymd == NY)
                {
                    down = 1;
                }

                
                if (!(right||left||up||down))

                {
                    double coeff = 2.0 * wi[i] * (1.0/cs) * (1.0/cs) * r[scalar_index(x, y)];
                    double tmpx = dirx[i] * (left * vl[0]  + right * vr[0] + up * vu[0] + down * vd[0]);
                    double tmpy = diry[i] * (left * vl[1]  + right * vr[1] + up * vu[1] + down * vd[1]);

                    ft[i] = f1[fieldn_index(x, y, index_opp[i])] + coeff * tmpx + coeff * tmpy;
                }

                else
                {
                    ft[i] = f1[fieldn_index(xmd, ymd, i)];
                }            
            }

            // load populations from adjacent nodes
            /*double ft1 = f1[fieldn_index(xm1, y, 1)];
            double ft2 = f1[fieldn_index(x, ym1, 2)];
            double ft3 = f1[fieldn_index(xp1, y, 3)];
            double ft4 = f1[fieldn_index(x, yp1, 4)];
            double ft5 = f1[fieldn_index(xm1, ym1, 5)];
            double ft6 = f1[fieldn_index(xp1, ym1, 6)];
            double ft7 = f1[fieldn_index(xp1, yp1, 7)];
            double ft8 = f1[fieldn_index(xm1, yp1, 8)];*/


            // compute moments
            double rho = ft[0] + ft[1] + ft[2] + ft[3] + ft[4] + ft[5] + ft[6] + ft[7] + ft[8];
            double rhoinv = 1.0 / rho;
            double ux = rhoinv * (ft[1] + ft[5] + ft[8] - (ft[3] + ft[6] + ft[7]));
            double uy = rhoinv * (ft[2] + ft[5] + ft[6] - (ft[4] + ft[7] + ft[8]));
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

            for (unsigned int i = 1; i < ndir; ++i)
            {
                // calculate dot product
                double cidotu = dirx[i] * ux + diry[i] * uy;
                // calculate equilibrium
                double feq = wi[i] * rho * (1.0 + 3.0 * cidotu + 4.5 * cidotu * cidotu - 1.5 * (ux * ux + uy * uy));
                // relax to equilibrium
                f2[fieldn_index(x, y, i)] =
                    omtauinv * ft[i] + tauinv * feq;
            }

            f0[field0_index(x, y)] = omtauinv * ft[0] + tw0r * (omusq);
            // double cidot3u = tux;
            // f2[fieldn_index(x, y, 1)] = omtauinv * ft[1] + twsr * (omusq + cidot3u * (1.0 + 0.5 * cidot3u));
            // ... similar expressions for directions 2-4
            //cidot3u = tux + tuy;
            // f2[fieldn_index(x, y, 5)] = omtauinv * ft[5] + twdr * (omusq + cidot3u * (1.0 + 0.5 * cidot3u));
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

/*void streamBBM(vector<double> &f_src, vector<double> &f_dst, vector<double> &r)    // stream with Bouncing Back Method BC
{
    for (unsigned int y = 0; y < NY; ++y)
    {
        for (unsigned int x = 0; x < NX; ++x)
        {
            for (unsigned int i = 0; i < ndir; ++i)
            {
                unsigned int xmd = (x - dirx[i]) ;
                unsigned int ymd = (y - diry[i]) ;
                
                if (!(0 <= xmd < NX && 0 <= ymd < NY))
                {
                    bounceback (f_src, f_dst, xmd, ymd, i, r);
                }

                else
                {
                    f_dst[field_index(x, y, i)] = f_src[field_index(xmd, ymd, i)];
                }            
            }
        }
    }
}

void bounceback (vector<double> &f_src, vector<double> &f_dst, unsigned int xmd, unsigned int ymd, unsigned int i, vector<double> &r)
{
    unsigned int x = xmd + dirx[i];
    unsigned int y = ymd + diry[i];

    unsigned int i_opp = index_opp[i]

    if (ymd != -1)
    {   
    f_dst[field_index(x, y, i)] = f_src[field_index(x, y, i_opp[i])];
    }

    else 
    {
    
           f_dst[field_index(x, y, i)] = f_src[field_index(x, y, i_opp)]- 2 * wi[i_opp] * r[scalar_index(x, y)] ;
    }



}*/