#ifndef SET_HAMILTONIAN_H
#define SET_HAMILTONIAN_H
#include <iostream>
#include <vector>
#include <cmath>
#include <armadillo>
#include <Eigen/Dense>

using namespace std;

class Set_Hamiltonian
{
public:

    int systemsize, no_of_states, no_of_hits, matrixsize, mysector;
    double J;
    bool armadillobool, sectorbool, sectorboolchanged, palhuse, testupip1downi, testdownip1upi;

    std::vector<int> sectorlist;
    std::vector<double> hs;

    Eigen::MatrixXd eigenH;
    arma::mat       armaH;

    //Initializers
    Set_Hamiltonian();
    Set_Hamiltonian(int systemsize, double J, vector<double> hs, bool armadillobool, bool sectorbool);

    // Sector functions
    void set_mysector(int mysector);
    void find_sector();      // Could merge these
    void trim_sectorlist();  // Could merge these
    void sector0();            // Do I need these if I am going to call set_hamiltonian from find_quantities?
    void sector1_2();          // Odd name, perhaps...

    void create_armadillo_matrix();
    void create_armadillo_matrix(int size);     // This is intended if we consider sectors
    void create_dense_Eigen_matrix();
    void create_dense_Eigen_matrix(int size);

    // Basic spin operations
    int give_spin(int i, int a);
    int set_spin_up(int i, int a);
    int set_spin_down(int i, int a);
    int flip_spin(int i, int a);      // Do I really need this one?

    // Spin operators
    double szi(int i, int a);
    double szip1szi(int i, int a);
    int upip1downi(int i, int a);  // A simpler version of spip1smi (if-test outside of function)
    int downip1upi(int i, int a);  // A simpler version of smip1spi (if-test outside of function)

    // State-specific functions
    int number_of_up(int a);
    int number_of_down(int a);
    void checktestupdown(int j, int i);

    //The systems
    // Well, I don't need these, do I?
    void randomize();
    void set_hs(std::vector<double> hs_in); // Maybe just this one. Or have it in the initializer/ constructor.
    void set_hs_hom();
    void set_hs_alt();


    //Hamiltonians: Different kinds of systems
    void set_elements(int i, int b);           // unspecified.
    void palhuse_set_elements(int i, int b);
    // Sector Hamiltonians
    void palhuse_interacting_sectorHamiltonian();
    void palhuse_diagonal_sectorHamiltonian();

    // Total Hamiltonian
    void palhuse_interacting_totalHamiltonian();
    void palhuse_diagonal_totalHamiltonian();

    void flip_diagonal_sectorHamiltonian();
};

#endif // SET_HAMILTONIAN_H
