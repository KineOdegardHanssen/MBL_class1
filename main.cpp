#include <iostream>
#include <vector>
#include <set_hamiltonian.h>
#include <diagonalize.h>
#include <find_quantities.h>

using namespace std;

void system_total_hom(int systemsize, int maxit, double tolerance, double h, double J, bool armabool);
void system_sector_hom(int systemsize, int maxit, double tolerance, double h, double J, bool armabool);
void system_total_random(int systemsize, int maxit, double tolerance, double h, double J, bool armabool);
void system_sector_random(int systemsize, int maxit, double tolerance, double h, double J, bool armabool);

int main()
{
    //const bool TRACE = false;
    int systemsize = 2;
    int maxit = 1e7;
    double J = 1;
    double h = 1;
    double tolerance = 1e-10;
    bool armabool = true;

    system_total_hom(systemsize, maxit, tolerance, h, J, armabool);

}

void system_total_hom(int systemsize, int maxit, double tolerance, double h, double J, bool armabool)
{
    bool sectorbool = false;
}

void system_total_random(int systemsize, int maxit, double tolerance, double h, double J, bool armabool)
{
    bool sectorbool = false;
}

void system_sector_hom(int systemsize, int maxit, double tolerance, double h, double J, bool armabool)
{
    bool sectorbool = true;

}

void system_sector_random(int systemsize, int maxit, double tolerance, double h, double J, bool armabool)
{
    bool sectorbool = true;

}
