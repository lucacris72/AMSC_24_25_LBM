#ifndef SAVE_CSV_H
#define SAVE_CSV_H

#include <vector>
#include <string>
#include <iostream>
#include <fstream>

void save_csv(const std::string &filename, const std::vector<double> &rho, const std::vector<double> &ux, const std::vector<double> &uy);

#endif
