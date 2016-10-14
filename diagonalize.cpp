#include "diagonalize.h"
#include<lapacke.h>
//#include <stdlib.h>
//#include <stdio.h>

Diagonalize::Diagonalize()
{
}

Diagonalize::Diagonalize(Systems given)
{
    this->given = given;
}

void Diagonalize::lapack_directly()
{
    const bool TRACE = false;

    //Tests if I want to use them
    //cout << "compute the LU factorization..." << endl << endl;

    /*
    cout << "What we send in to LAPACK_dysyev" << endl;

    for(int i=0; i<N*N; i++)
    {
        cout << A[i] << " ";
    }
    cout << endl;

    cout << "W" << endl;
    for(int i=0; i<N; i++)
    {
        cout << W[i] << " ";
    }
    cout << endl;

    cout << "JOBS: " << JOBS << "; UPLO: " << UPLO << endl;
    */

    N = given.matrixsize;

    // note, to understand this part take a look in the MAN pages, at section of parameters.


    JOBS = 'V';
    //if(TRACE)          cout << "JOBS OK," << endl;
    UPLO = 'L';
    //if(TRACE)          cout << "UPLO OK," << endl;
    LDA  = N;
    //if(TRACE)          cout << "LDA OK," << endl;
    int LWORK = 3*N-1;
    //if(TRACE)          cout << "LWORK OK," << endl;
    int INFO;
    //if(TRACE)          cout << "INFO declared," << endl;

    std::vector<double> WORK  = std::vector<double> (LWORK);
    //if(TRACE)          cout << "WORK OK," << endl;
    //std::vector<double> RWORK = std::vector<double>(3*N-2);
    std::vector<double> W     = std::vector<double>(N);
    //if(TRACE)          cout << "W OK," << endl;

    //if(TRACE)    cout << "Before declaration double *A1" << endl;
    double *A = given.eigenH.data(); // Must I rewrite this?
    //if(TRACE)    cout << "After declaration double *A1" << endl;

    // end of declarations

    double start_time_lapack_directly = clock();
    LAPACK_dsyev(&JOBS, &UPLO, &N, A, &LDA, &W[0],&WORK[0],&LWORK, &INFO);
    double end_time_lapack_directly = clock();

    /*
    if(TRACE)    cout << "Have run LAPACK_dsyev" << endl;

    cout << "What we got out of LAPACK_dysyev" << endl;

    for(int i=0; i<N*N; i++)
    {
        cout << A[i] << " ";
    }
    cout << endl;

    cout << "W" << endl;
    for(int i=0; i<N; i++)
    {
        cout << W[i] << " ";
    }
    cout << endl;
    */


    // checks INFO, if INFO != 0 something goes wrong, for more information see the MAN page of dgetrf.
    if(INFO!=0)
    {
        cout << "an error occured : "<< INFO << endl << endl;
    }else{
        cout << "system solved..."<< endl << endl;
        lapacktime = (end_time_lapack_directly - start_time_lapack_directly)/CLOCKS_PER_SEC;
        if(INFO!=0)
        {
            // checks INFO, if INFO != 0 something goes wrong, for more information see the MAN page of dgetrs.
            cout << "an error occured : "<< INFO << endl << endl;
            lapacktime = 1e16;
        }else{
            /*
            cout << "Using LAPACK directly" << endl;
            cout << "Eigenvectors : " << endl;
            for (int i=0;i<N;i++)
            {
                for(int j=0; j<N; j++)
                {
                    cout << A[i+j*LDA] << " ";
                }
                cout << endl;
            }
            cout << endl;
            */

            cout << "Corresponding eigenvalues: " << endl;
            for(int k=0; k<N; k++)    cout << W[k] << ", ";
            cout << endl;
        }
    }


    /*
    cout << "Comparing matrices found by LAPACK: " << endl;
    cout << "Output from LAPACK: " << endl;
    for (int i=0;i<N;i++)
    {
        for(int j=0; j<N; j++)
        {
            cout << A[i+j*LDA] << " ";
        }
        cout << endl;
    }
    cout << endl;

    cout << "Corresponding eigenvalues: " << endl;
    for(int k=0; k<N; k++)    cout << W[k] << ", ";
    cout << endl;
    */

    // Setting things manually because it is a big noisy mess otherwise
    // Okay, time both...
    eigenmatrix_H = Eigen::MatrixXd(N,N);
    eigenvalues_H = Eigen::VectorXd(N);



    for (int i=0;i<N;i++)
    {
        for(int j=0; j<N; j++)            eigenmatrix_H(i,j) = A[i+j*LDA];
    }

    for(int k=0; k<N; k++)                eigenvalues_H(k) = W[k];

    //cout << "After putting data manually into an Eigen matrix: " << endl;
    //cout << eigenmatrix_H << endl;

    /*
    cout << "Manually: " << endl;
    cout << "Eigenmatrix: " << endl;
    cout << eigenmatrix_H << endl;

    cout << "Eigenvalues" << endl;
    cout << eigenvalues_H << endl;
    */




    //cout << "After using Eigen::Map: " << endl;
    Eigen::Map<Eigen::MatrixXd> B(A, N, N);
    //Eigen::Map<Eigen::VectorXd> evs(W);

    //eigenvalues_H = evs;
    //eigenmatrix_H = B;

    //eigenmatrix_H(A, N, N);

    //cout << "After using the Map function in Eigen: " << endl;
    cout << B << endl;

    if(TRACE)    cout << "Everything should be running just fine" << endl;

}




void Diagonalize::using_armadillo()
{
    N = given.matrixsize;

    eigenvalues_armadillo = arma::vec(N);
    eigenvectors_armadillo = arma::mat(N,N);

    //string method = "std";
    double start_time_arma = clock();        // If I want to take the time
    arma::eig_sym(eigenvalues_armadillo,eigenvectors_armadillo, given.armaH);
    double end_time_arma = clock();
    armatime = (end_time_arma - start_time_arma)/CLOCKS_PER_SEC;
}


void Diagonalize::print_using_armadillo()
{
    for(int i= 0; i<N; i++)
    {
        cout << "Eigenvalue = " << eigenvalues_armadillo(i) << endl;
        cout << "Its corresponding eigenvector:" << endl;
        for(int j=0; j<N; j++)       cout << eigenvectors_armadillo(j,i) << " ";
        cout << endl;
    }   // End for-loop over i

    cout << "Eigenvalues = " << eigenvalues_armadillo << endl;
    cout << "Eigenvectors : " << endl << eigenvectors_armadillo << endl;
}


void Diagonalize::using_dense_eigen()
{
    N = given.matrixsize;

    Eigen::EigenSolver<Eigen::MatrixXd> es(given.eigenH);

    //eigenvalues_H = es.eigenvalues();
    //eigenmatrix_H = es.eigenvectors();
}


void Diagonalize::print_dense_using_eigen()
{
    cout << "In print_dense_using_eigen" << endl;
    N = given.matrixsize;

    //cout << "The matrix is: " << endl;
    //cout << given.eigenH << endl;

    double eigen_time_start = clock();
    Eigen::EigenSolver<Eigen::MatrixXd> es(given.eigenH);
    double eigen_time_end = clock();

    eigen_time = (eigen_time_end - eigen_time_start)/CLOCKS_PER_SEC;

    cout << "The eigenvalues of H are:" << endl << es.eigenvalues() << endl;
    //cout << "The matrix of eigenvectors, V, is:" << endl << es.eigenvectors() << endl << endl;

    cout << "Exiting print_dense_using_eigen" << endl;
}



