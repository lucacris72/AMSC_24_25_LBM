#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <memory>
#include <string>

#include "LBM.h"
#include "seconds.h"

using namespace std;

void report_flow_properties(unsigned int t, vector<double> rho, vector<double> ux, vector<double> uy);

int main(int argc, char* argv[])
{
    vector<double> f0;
    vector<double> f1;
    vector<double> f2;
    vector<double> rho;
    vector<double> ux;
    vector<double> uy;

    try {
        f0.reserve(size0);
        f1.reserve(size1);
        f2.reserve(size1);
        rho.reserve(size0);
        ux.reserve(size0);
        uy.reserve(size0);
    } catch (const std::bad_alloc& e) {
        std::cerr << "Memory allocation failed: " << e.what() << std::endl;
        exit(-1);
    }

    double bytesPerMiB = 1024.0*1024.0;
    double bytesPerGiB = 1024.0*1024.0*1024.0;
    
    
    // compute flow at t=0 
    // to initialise rho, ux, uy fields.
    //Initialization function not written yet
    
    // initialise f1 as equilibrium for rho, ux, uy
    init_equilibrium(f0,f1,rho,ux,uy);
    
    save_scalar("rho", rho, 0);
    save_scalar("ux", ux, 0);
    save_scalar("uy", uy, 0);

    double start = seconds();
    
    // main simulation loop; take NSTEPS time steps
    for(unsigned int n = 0; n < NSTEPS; ++n)
    {
        bool save = (n+1)%NSAVE == 0;
        bool msg  = (n+1)%NMSG == 0;
        bool need_scalars = save || (msg && computeFlowProperties);
        
        // stream and collide from f1 storing to f2
        // optionally compute and save moments
        stream_collide_save(f0,f1,f2,rho,ux,uy,need_scalars);
        
        if(save)
        {
            save_scalar("rho",rho,n+1);
            save_scalar("ux", ux, n+1);
            save_scalar("uy", uy, n+1);
        }
        
        // swap pointers
        f1.swap(f2);
        
        if(msg)
        {
            if(computeFlowProperties)
            {
                report_flow_properties(n+1, rho, ux, uy);
            }
            
            if(!quiet)
                printf("completed timestep %d\n",n+1);
        }
    }

    double end = seconds();
    double runtime = end-start;

    size_t doubles_read = ndir; // per node every time step
    size_t doubles_written = ndir;
    size_t doubles_saved = 3; // per node every NSAVE time steps
    
    // note NX*NY overflows when NX=NY=65536
    size_t nodes_updated = NSTEPS*size_t(NX*NY);
    size_t nodes_saved   = (NSTEPS/NSAVE)*size_t(NX*NY);
    double speed = nodes_updated/(1e6*runtime);
    
    double bandwidth = (nodes_updated*(doubles_read + doubles_written)+nodes_saved*(doubles_saved))*sizeof(double)/(runtime*bytesPerGiB);
    
    printf(" ----- performance information -----\n");
    printf(" memory allocated: %.1f (MiB)\n",total_mem_bytes/bytesPerMiB);
    printf("        timesteps: %u\n",NSTEPS);
    printf("          runtime: %.3f (s)\n",runtime);
    printf("            speed: %.2f (Mlups)\n",speed);
    printf("        bandwidth: %.1f (GiB/s)\n",bandwidth);
    
    return 0;
}

void report_flow_properties(unsigned int t, vector<double> rho, vector<double> ux, vector<double> uy)
{
    vector<double> prop;
    prop.reserve(4);
    compute_flow_properties(t,rho,ux,uy,prop);
    printf("%u,%g,%g,%g,%g\n",t,prop[0],prop[1],prop[2],prop[3]);
}

void save_scalar(const string name, vector<double> scalar, unsigned int n)
{
    // assume reasonably-sized file names
    char filename[128];
    char format[16];
    
    // compute maximum number of digits
    int ndigits = floor(log10((double)NSTEPS)+1.0);
    
    // generate format string
    // file name format is name0000nnn.bin
    sprintf(format,"%%s%%0%dd.bin",ndigits);
    sprintf(filename,format,name,n);
    
    // open file for writing
    FILE *fout = fopen(filename,"wb+");
    
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