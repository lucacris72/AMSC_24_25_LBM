#include "save_csv.h"
#include "LBM.h"

void save_csv(const std::string &filename, const std::vector<double> &rho, const std::vector<double> &ux, const std::vector<double> &uy)
{
    std::ofstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Errore opening file: " << filename << std::endl;
        return;
    }

    // Head of CSV
    file << "x,y,rho,ux,uy" << std::endl;

    // Save data for every point of the lattice
    for (unsigned int y = 0; y < NY; ++y)
    {
        for (unsigned int x = 0; x < NX; ++x)
        {
            size_t idx = scalar_index(x, y);
            file << x << "," << y << "," << rho[idx] << "," << ux[idx] << "," << uy[idx] << std::endl;
        }
    }

    file.close();
}