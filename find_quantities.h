#ifndef FIND_QUANTITIES_H
#define FIND_QUANTITIES_H
#include <iostream>
#include <set_hamiltonian.h>
#include <diagonalize.h>
#include <armadillo>
#include <Eigen/Dense>
#include<vector>
#include<complex>


class Find_Quantities
{
public:

    // One class on top of this, like runnit for 2p?


    Set_Hamiltonian system;            // Should I do something like this and nest them?  // Worry about this later
    Diagonalize eig_all;

    char field_type;
    int N, maxit, systemsize, li, lh;
    double J, h, Z, beta, min_ev, tolerance; // Only change smallest_ev this for every new instance of quantities...
    bool armadillobool, sectorbool, inftempbool, field_type_fail;

    std::vector<double> hs;

    Eigen::VectorXd eigenvalues_all_Eigen;
    Eigen::VectorXd eigvals;
    Eigen::MatrixXd eigmat;

    arma::vec eigenvalues_all_arma;
    arma::vec eigvals_a;  // This will probably be a problem...
    arma::mat eigmat_a;   // Could rename it at each step. Troublesome...


    // Initializer
    Find_Quantities(char field_type, int maxit, int systemsize, double tolerance, double J, double h, bool armadillobool, bool sectorbool, bool inftempbool);
    void initialize_sector(); // Are these two only needed for ETH? In that case, I guess they shouldn't be called from the constructor.
    void initialize_all();
    void make_hs_random();
    void make_hs_homogenous();
    void make_hs_alternating();




    //Functions

    // Basic functions
    void sort_energies();
    int middle_sector();
    arma::mat initialize_matrix_arma(int size);
    Eigen::MatrixXd initialize_matrix_Eigen(int size);
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
    Eigen::MatrixXd trace_Eigen(Eigen::MatrixXd A);
    Eigen::MatrixXd trace_Eigen_maziero(Eigen::MatrixXd A);
    Eigen::MatrixXd trace_Eigen_sector(Eigen::MatrixXd A);
    Eigen::MatrixXd thermalmat_Eigen();             // See if I change this a bit.
    Eigen::MatrixXd eigenstatemat_Eigen(int i);

    // Armadillo
    arma::mat trace_arma(arma::mat A);
    arma::mat trace_arma_maziero(arma::mat A);
    arma::mat trace_arma_sector(arma::mat A);
    arma::mat thermalmat_arma();                    // See if I change this a bit.
    arma::mat eigenstatemat_arma(int i);
};

#endif // FIND_QUANTITIES_H
