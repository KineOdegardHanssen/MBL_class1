#ifndef FIND_QUANTITIES_H
#define FIND_QUANTITIES_H
#include <iostream>
#include <set_hamiltonian.h>
#include <diagonalize.h>
#include <armadillo>
#include <Eigen/Dense>
#include<vector>
#include<complex>

using namespace std;

class Find_Quantities
{
public:

    // One class on top of this, like runnit for 2p?


    Set_Hamiltonian system;            // Should I do something like this and nest them?  // Worry about this later
    Diagonalize eig_all;

    char field_type, field_type_x;
    int N, n_all, maxit, systemsize, middlesector, li, lh;
    double J, h, hx, Z, beta, min_ev, tolerance; // Only change smallest_ev this for every new instance of quantities...
    bool armadillobool, sectorbool, inftempbool, field_type_fail, field_type_fail_x;

    std::vector<double> hs, hxs;

    // The eigenvalues and eigenvectors for the sector, or all of H if we choose not to work with sectors
    Eigen::VectorXd eigenvalues_all_Eigen;
    Eigen::MatrixXd eigenvectors_all_Eigen;

    // Every eigenvalue and eigenvector. For when we work with sectors.
    Eigen::VectorXd eigvals;
    Eigen::MatrixXd eigmat;

    // The eigenvalues and eigenvectors for the sector, or all of H if we choose not to work with sectors
    arma::vec eigvals_a;
    arma::mat eigmat_a;

    // Every eigenvalue and eigenvector. For when we work with sectors.
    arma::vec eigenvalues_all_arma;
    arma::mat eigenvectors_all_arma;




    // Initializer
    Find_Quantities();
    Find_Quantities(char field_type, int maxit, int systemsize, double tolerance, double J, double h, bool armadillobool, bool sectorbool, bool inftempbool);
    //Find_Quantities::Find_Quantities(char field_type, char field_type_x, int maxit, int systemsize, double tolerance, double J, double h, double hx, bool armadillobool, bool inftempbool);
    void spinnotconserved(char field_type, char field_type_x, int maxit, int systemsize, double tolerance, double J, double h, double hx, bool armadillobool, bool inftempbool);
    void initialize_sector(); // Are these two only needed for ETH? In that case, I guess they shouldn't be called from the constructor.
    void initialize_all();
    void initialize_all_withsx();
    double testquantumside(int i);
    double testquantumside_arma(int i);
    double testquantumside_Eigen(int i);
    void make_hs_random();
    void make_hs_homogenous();
    void make_hs_alternating();
    void set_hs_manually(vector<double> hs_in);
    void make_hxs_random();
    void make_hxs_homogenous();
    void make_hxs_alternating();
    void make_hxs_linspace();
    void set_hxs_manually(vector<double> hs_in);


    //Functions

    // Basic functions
    void sort_energies();
    int middle_sector();
    int factorial(int i);
    arma::mat initialize_matrix_arma(int size);        // armadillo has a function that does this.
    Eigen::MatrixXd initialize_matrix_Eigen(int size); // maybe the same goes for Eigen?
    int signcompare(double fa, double fc);
    void calculateZ();
    void calculateZ_arma();

    // Finding beta
    // To be run for each eigenstate/eigenvalue
    // I may be giving up some speed or accuracy here: can divide by e^(-beta E_lowest), but not set e^(-beta(E_lowest-E_lowest)) to 1 automatically
    void newtonsmethod(double eigenvalue);
    void bisectionmethod(double eigenvalue);
    double self_consistency_beta(double eigenvalue, double betatest);
    double self_consistency_beta_derivative(double eigenvalue, double betatest);
    // armadillo
    void newtonsmethod_arma(double eigenvalue);
    void bisectionmethod_arma(double eigenvalue);
    double self_consistency_beta_a(double eigenvalue, double betatest);
    double self_consistency_beta_derivative_a(double eigenvalue, double betatest);

    // Eigenstate Thermalization Hypothesis
    double ETH(int i);       // Should this really be a void?
    double ETH_arma(int i);
    double ETH_arma_sector(int i);
    double ETH_arma_maziero(int i);
    double ETH_Eigen(int i);
    double ETH_Eigen_sector(int i);
    double ETH_Eigen_maziero(int i);

    // Eigen
    Eigen::MatrixXd trace_Eigen(Eigen::MatrixXd A, int size);
    Eigen::MatrixXd trace_Eigen_maziero(Eigen::MatrixXd A, int size);
    Eigen::MatrixXd trace_Eigen_sector(Eigen::MatrixXd A);
    Eigen::MatrixXd thermalmat_Eigen();             // See if I change this a bit.
    Eigen::MatrixXd eigenstatemat_Eigen(int i);

    // Armadillo
    arma::mat trace_arma(arma::mat A, int size);
    arma::mat trace_arma_maziero(arma::mat A, int size);
    arma::mat trace_arma_sector(arma::mat A);
    arma::mat thermalmat_arma();                    // See if I change this a bit.
    arma::mat eigenstatemat_arma(int i);
};

#endif // FIND_QUANTITIES_H
