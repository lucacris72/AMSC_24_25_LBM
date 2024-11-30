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
            
            // initial value for the macroscopic proprieties
            
            double rho = rho0;
            double ux = 0.0;
            double uy = 0.0;

            if (y == 0)
            {
                ux = vu[0]; // apply boundary conditions on upper wall
            } 
            
            r[scalar_index(x, y)] = rho;
            u[scalar_index(x, y)] = ux;
            v[scalar_index(x, y)] = uy;

            f0[field0_index(x,y)] = w0*(1.0 - 1.5*(ux*ux+uy*uy)); // initialize the distributions f1
            for (unsigned int i = 1; i < ndir; ++i)
            {
                double cidotu = dirx[i] * ux + diry[i] * uy;
                f1[fieldn_index(x, y, i)] = wi[i] * rho * (1.0 + 3.0 * cidotu + 4.5 * cidotu * cidotu - 1.5 * (ux * ux + uy * uy));
            }

           // printf("x = %d, y = %d, NY = %d \n", x, y, NY);

            
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
            
            double rho_in = f0[field0_index(x, y)];
            
            for (unsigned int i = 1; i < ndir; ++i)
            {
                rho_in = rho_in + f1[fieldn_index(x, y, i)];
            }
            

            ft[0] = f0[field0_index(x, y)];

            for (unsigned int i = 1; i < ndir; ++i)
            {
                unsigned int xmd = (x - dirx[i]) ;
                unsigned int ymd = (y - diry[i]) ;

                bool right = 0;
                bool left = 0;
                bool up = 0;
                bool down = 0;

                              
                if (xmd < 0) // check if I need a value which is outside the border
                {
                    left = 1;
                }

                if (xmd > NX -1)
                {
                    right = 1;
                }

                if (ymd < 0)
                {
                    up = 1;
                }

                if (ymd > NY - 1)
                {
                    down = 1;
                }

                
                if ((right||left||up||down)) // apply boundary conditions if I'm going outside the border

                {
                    double coeff = 2.0 * wi[i] * (1.0/cs) * (1.0/cs) * rho_in;
                    double tmpx = dirx[i] * (left * vl[0]  + right * vr[0] + up * vu[0] + down * vd[0]);
                    double tmpy = diry[i] * (left * vl[1]  + right * vr[1] + up * vu[1] + down * vd[1]);

                    ft[i] = f1[fieldn_index(x, y, index_opp[i])] + coeff * tmpx + coeff * tmpy;
                }

                else // streaming step if when I'm not going outside the border
                {
                    ft[i] = f1[fieldn_index(xmd, ymd, i)];
                }            
            }

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
    fwrite(&scalar[0], 1, mem_size_scalar, fout);
    // close file
    fclose(fout);

    if(ferror(fout))
    {
        fprintf(stderr,"Error saving to %s\n",filename);
        perror("");
    }
    else
    {
        if(!quiet)
            printf("Saved to %s\n",filename);
    }
}