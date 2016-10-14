#ifndef DIAGONALIZE_H
#define DIAGONALIZE_H
#include <iostream>
#include <set_hamiltonian.h>
#include <armadillo>
#include <Eigen/Dense>
#include <Eigen/QR>
#include <Eigen/Eigenvalues>
#include<vector>
#include<complex>


using namespace std;
//using namespace arma;

class Diagonalize
{
public:

    Set_Hamiltonian given;

    // For linking directly to LAPACK
    char JOBS;
    char UPLO;  // upper triangular matrix is stored after diag.
    int N;      // size of matrix
    int LDA;


    // For armadillo solvers
    double armatime, eigen_time, lapacktime;
    arma::vec eigenvalues_armadillo; // Should I use arma::vec? Probably a good idea. Eigen probably has something similar
    arma::mat eigenvectors_armadillo;

    // For Eigen using dense matrices
    Eigen::VectorXd eigenvalues_H;      // Not using these (yet...)
    Eigen::MatrixXd eigenmatrix_H;

    //Eigen::Map<Eigen::VectorXd> eigenvalues_H;
    //Eigen::Map<Eigen::MatrixXd> eigenmatrix_H;


    //Initializers
    Diagonalize();
    Diagonalize(Systems given);



    // Functions for using LAPACK
    void lapack_directly();

    // Functions for armadillo matrices (dense)
    void using_armadillo();
    void print_using_armadillo();

    void using_dense_eigen();
    void print_dense_using_eigen();

    // Functions for Eigen::SparseMatrix
    void print_sparse_Eigen();
};

#endif // DIAGONALIZE_H
