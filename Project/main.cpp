#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <memory>

#include "LBM.cpp"

using namespace std;

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
    
    
    // compute flow at t=0 
    // to initialise rho, ux, uy fields.
    //Initialization function not written yet
    
    // initialise f1 as equilibrium for rho, ux, uy
    init_equilibrium(f0,f1,rho,ux,uy);
    
    save_scalar("rho", rho, 0);
    save_scalar("ux", ux, 0);
    save_scalar("uy", uy, 0);
    
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
                //I need to know how to report data for visualization before writing this function
                //report_flow_properties();
            }
            
            if(!quiet)
                printf("completed timestep %d\n",n+1);
        }
    }
    
    return 0;
}