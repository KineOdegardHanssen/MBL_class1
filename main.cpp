#include <iostream>
#include <vector>
#include <set_hamiltonian.h>
#include <diagonalize.h>
#include <find_quantities.h>

using namespace std;

void test_sectormatrix_diag(int systemsize, vector<double> hs, double J);
void test_totalmatrix_diag(int systemsize, vector<double> hs, double J);
void test_totalmatrix_spinnotconserved(int systemsize, vector<double> hs, vector<double> hxs, double J);

void system_total_hom(int systemsize, int maxit, double tolerance, double h, double J, bool armabool, bool inftempbool);
void system_sector_hom(int systemsize, int maxit, double tolerance, double h, double J, bool armabool);
void system_total_random(int systemsize, int maxit, double tolerance, double h, double J, bool armabool);
void system_sector_random(int systemsize, int maxit, double tolerance, double h, double J, bool armabool);

void system_total_spinnotconserved_zhom_xalt(int systemsize, int maxit, double tolerance, double h, double hx, double J, bool armabool, bool inftempbool);

int factorial(int i);

int main()
{
    //const bool TRACE = false;
    int systemsize = 3;
    int maxit = 1e7;
    double J = 1;
    double h = 1;
    double hx = 0.5;
    double tolerance = 1e-10;
    bool armabool = true;
    bool inftempbool = false;

    vector<double> hs = vector<double>(systemsize);
    hs[0] = 0.2;
    hs[1] = 0.7;
    hs[2] = 0.3;

    vector<double> hxs = vector<double>(systemsize);
    hxs[0] = 0.3;
    hxs[1] = 0.1;
    hxs[2] = 0.9;

    //test_totalmatrix_diag(systemsize, hs, J);
    //system_total_hom(systemsize, maxit, tolerance, h, J, armabool);
    //test_totalmatrix_diag(systemsize, hs, J);
    test_totalmatrix_spinnotconserved(systemsize, hs, hxs, J);
    //system_total_spinnotconserved_zhom_xalt(systemsize, maxit, tolerance, h, hx, J, armabool, inftempbool);

    /*
    cout << "Playing with factorials! " << endl;
    cout << "i = " << "; factorial(i) = " << endl;
    cout << 0 << "; "  << factorial(0) << endl;
    cout << 1 << "; "  << factorial(1) << endl;
    cout << 2 << "; "  << factorial(2) << endl;
    cout << 3 << "; "  << factorial(3) << endl;
    cout << 4 << "; "  << factorial(4) << endl;
    cout << 5 << "; "  << factorial(5) << endl;
    cout << 6 << "; "  << factorial(6) << endl;
    cout << 7 << "; "  << factorial(7) << endl;
    cout << 8 << "; "  << factorial(8) << endl;
    cout << 9 << "; "  << factorial(9) << endl;
    cout << 10 << "; " << factorial(10) << endl;
    cout << -1 << "; " << factorial(-1) << endl;
    */




}

void test_sectormatrix_diag(int systemsize, vector<double> hs, double J)
{
    int sector;
    bool sectorbool = true;
    bool armabool = true;

    for(int i=0; i<systemsize; i++)    cout << hs[i] << " " << endl;
    cout << "Solving using armadillo: " << endl;
    Set_Hamiltonian system(systemsize, J, hs, armabool, sectorbool);

    // Setting the sector
    if(systemsize%2==0)    sector = systemsize/2;
    else                   sector = (systemsize+1)/2;  // Or possibly give the sector?
    system.set_mysector(sector);
    system.palhuse_interacting_sectorHamiltonian();
    system.palhuse_diagonal_sectorHamiltonian();

    Diagonalize diagon(system);

    diagon.using_armadillo();
    diagon.print_using_armadillo();

    for(int i=0; i<systemsize; i++)    cout << hs[i] << " " << endl;
    cout << "Solving using Eigen and LAPACK: " << endl;
    armabool = false;
    Set_Hamiltonian system2(systemsize, J, hs, armabool, sectorbool);

    // Setting the sector
    if(systemsize%2==0)    sector = systemsize/2;
    else                   sector = (systemsize+1)/2;  // Or possibly give the sector?
    system2.set_mysector(sector);
    system2.palhuse_interacting_sectorHamiltonian();
    system2.palhuse_diagonal_sectorHamiltonian();

    Diagonalize diagon2(system2);

    diagon2.lapack_directly();

}


void test_totalmatrix_diag(int systemsize, vector<double> hs, double J)
{
    bool sectorbool = false;
    bool armabool = true;

    for(int i=0; i<systemsize; i++)    cout << hs[i] << " " << endl;
    cout << "Solving using armadillo: " << endl;
    Set_Hamiltonian system(systemsize, J, hs, armabool, sectorbool);

    system.palhuse_interacting_totalHamiltonian();
    system.palhuse_diagonal_totalHamiltonian();

    Diagonalize diagon(system);

    diagon.using_armadillo();
    diagon.print_using_armadillo();

    for(int i=0; i<systemsize; i++)    cout << hs[i] << " " << endl;
    cout << "Solving using Eigen and LAPACK: " << endl;
    armabool = false;
    Set_Hamiltonian system2(systemsize, J, hs, armabool, sectorbool);

    system2.palhuse_interacting_totalHamiltonian();
    system2.palhuse_diagonal_totalHamiltonian();

    Diagonalize diagon2(system2);

    diagon2.lapack_directly();

    cout << "One eigenvector: " << endl;
    for(int i=0; i<diagon.N; i++)   cout << diagon.eigenvectors_armadillo(i,1) << " " << endl;
}



void test_totalmatrix_spinnotconserved(int systemsize, vector<double> hs, vector<double> hxs, double J)
{

    bool sectorbool = false;
    bool armabool = true;

    for(int i=0; i<systemsize; i++)    cout << hs[i] << " " << endl;
    cout << "Solving using armadillo: " << endl;
    Set_Hamiltonian system(systemsize, J, hs, armabool, sectorbool);
    system.give_hxs(hxs);

    cout << "This is what armabool is, according to system: " << system.armadillobool << endl;
    cout << "False is: " << false << endl;

    system.palhuselike_spinnotconserved_interacting_totalHamiltonian();
    system.palhuse_diagonal_totalHamiltonian();

    Diagonalize diagon(system);

    diagon.using_armadillo();
    diagon.print_using_armadillo();

    for(int i=0; i<systemsize; i++)    cout << hs[i] << " " << endl;
    cout << "Solving using Eigen and LAPACK: " << endl;

    armabool = false;
    Set_Hamiltonian system2(systemsize, J, hs, armabool, sectorbool);
    system2.give_hxs(hxs);

    cout << "This is what armabool is, according to system2: " << system2.armadillobool << endl;
    cout << "False is: " << false << endl;

    system2.palhuselike_spinnotconserved_interacting_totalHamiltonian();
    system2.palhuse_diagonal_totalHamiltonian();

    Diagonalize diagon2(system2);

    diagon2.lapack_directly();

    cout << "One eigenvector: " << endl;
    for(int i=0; i<diagon.N; i++)   cout << diagon.eigenvectors_armadillo(i,1) << " " << endl;
}


void system_total_hom(int systemsize, int maxit, double tolerance, double h, double J, bool armabool, bool inftempbool)
{
    bool sectorbool = false;

    char field_type = 'H';

    Find_Quantities zyztehm(field_type, maxit, systemsize, tolerance, J, h, armabool, sectorbool, inftempbool);
    //for(int i=0; i<zyztehm; i++)
    zyztehm.ETH(1);

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



void system_total_spinnotconserved_zhom_xalt(int systemsize, int maxit, double tolerance, double h, double hx, double J, bool armabool, bool inftempbool)
{   // Should probably vary the systemsize here...
    /*
    char field_type = 'H';
    char field_type_x = 'A';

    Find_Quantities zyztehm(field_type, field_type_x, maxit, systemsize, tolerance, J, h, hx, armabool, inftempbool);
    //zyztehm.spinnotconserved(field_type, field_type_x, maxit, systemsize, tolerance, J, h, hx, armabool, inftempbool);
    for(int i=0; i<zyztehm.N; i++)
    {   // Should I print to file or something?
        cout << zyztehm.testquantumside(i) << endl;
    }
    */
}

int factorial(int i)
{
    double j=1;
    if(i<0)
    {
        cout << "Error! factorial(int i) must take input i>0!" << endl;
    }
    else if(i>1)
    {
      for(int k=2; k<=i; k++)    j *= k;
    }

    return    j;
}
