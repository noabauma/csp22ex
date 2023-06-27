#include<iostream>
#include<fstream>
#include<string>
#include<cmath>
#include<random>
#include<Eigen/Dense>
#include <chrono>
#include"lattice.hpp"

const static Eigen::IOFormat CSVFormat(Eigen::StreamPrecision, Eigen::DontAlignCols, ",", "\n");

//global defined constants
int L;
const int J = 1;            //coupling constant
double p;                   //prob to add to cluster

std::mt19937 gen(42);
std::uniform_int_distribution<> uniform;
std::uniform_real_distribution<double> uniform_real(0.0, 1.0);

/*
Args:

Updates:
    Recursivly visits sites and updates them if they are part of the cluster
*/
template<typename T>
void visit(lattice_3D<T>& x, int& M, int& E, lattice_3D<bool>& visited, const T state, int i, int j, int k){

    //add to cluster with prob
    if(uniform_real(gen) < p){
        visited(i,j,k) = true; //site marked as visited (being part of the cluster!)

        const int x_old = x(i,j,k);
        const int nn = x_old*x.nn_sum(i,j,k);

        x(i,j,k) *= -1; //flip spin

        M -= 2*x_old;
        E += 2*J*nn;

        //now that the top site is part of the cluster we go visit other sites
        const int i_p = (i+1)%L;
        const int i_n = (i-1+L)%L;
        const int j_p = (j+1)%L;
        const int j_n = (j-1+L)%L;
        const int k_p = (k+1)%L;
        const int k_n = (k-1+L)%L;

        if(!visited(i_p,j,k) && x(i_p,j,k) == state){
            visit(x, M, E, visited, state, i_p, j, k);
        }

        if(!visited(i_n,j,k) && x(i_n,j,k) == state){
            visit(x, M, E, visited, state, i_n, j, k);
        }

        if(!visited(i,j_p,k) && x(i,j_p,k) == state){
            visit(x, M, E, visited, state, i, j_p, k);
        }

        if(!visited(i,j_n,k) && x(i,j_n,k) == state){
            visit(x, M, E, visited, state, i, j_n, k);
        }

        if(!visited(i,j,k_p) && x(i,j,k_p) == state){
            visit(x, M, E, visited, state, i, j, k_p);
        }

        if(!visited(i,j,k_n) && x(i,j,k_n) == state){
            visit(x, M, E, visited, state, i, j, k_n);
        }
        
    }
}

/*
Args:
    x: Spin configuration
    M: Magnetization of x
    E: Energy of x

Updates:
    x, M and E after one Monte Carlo move
*/
template<typename T>
void move(lattice_3D<T>& x, int& M, int& E){
    lattice_3D<bool> visited(L);    //all automatically false (pls check)

    //pick one site at random
    const int i = uniform(gen);
    const int j = uniform(gen);
    const int k = uniform(gen);

    //check for sites with same state
    const T state = x(i,j,k);

    visited(i,j,k) = true;

    const int x_old = x(i,j,k);
    const int nn = x_old*x.nn_sum(i,j,k);

    x(i,j,k) *= -1; //flip with 100% prob (because it is the first)

    M -= 2*x_old;
    E += 2*J*nn;

    const int i_p = (i+1)%L;
    const int i_n = (i-1+L)%L;
    const int j_p = (j+1)%L;
    const int j_n = (j-1+L)%L;
    const int k_p = (k+1)%L;
    const int k_n = (k-1+L)%L;

    if(x(i_p,j,k) == state){
        visit(x, M, E, visited, state, i_p, j, k);
    }

    if(x(i_n,j,k) == state){
        visit(x, M, E, visited, state, i_n, j, k);
    }

    if(x(i,j_p,k) == state){
        visit(x, M, E, visited, state, i, j_p, k);
    }

    if(x(i,j_n,k) == state){
        visit(x, M, E, visited, state, i, j_n, k);
    }

    if(x(i,j,k_p) == state){
        visit(x, M, E, visited, state, i, j, k_p);
    }

    if(x(i,j,k_n) == state){
        visit(x, M, E, visited, state, i, j, k_n);
    }
}



int main(int argc, char *argv[]){
    if(argc > 3 || argc < 3){
        printf("Give temperature T in range (0.0, 10.0), system size L\n");
        return 1;
    }

    const double T = atof(argv[1]);
    const double beta = 1.0/T;              //NOTE: beta = 1/(kbT) with kb = 1.0, hence inverse of Temp
    p = 1.0 - std::exp(-2.0*beta*J);

    L = atoi(argv[2]);
    const unsigned N     = L*L*L;
    const double   N_inv = 1.0/double(N);
    uniform = std::uniform_int_distribution<>(0,L-1);

    const unsigned Nsample = 20*N;          //number of samples (= size of the Markov chain)
    
    //initialize lattice matrix
    lattice_3D<int16_t> x(L);               //calloc (initiallized with zeros)
    x.fill(1);                              //fill with ones
    int M = N;
    int E = -3*J*N;

    /*
    //not needed for correlation function
    for(int i = 0; i < L; ++i){
        for(int j = 0; j < L; ++j){
            for(int k = 0; k < L; ++k){
                if(uniform_real(gen) < 0.5){
                    x(i,j,k) = -1;
                    M -= 2;
                    E += 2*J*x.nn_sum(i, j, k);
                }
            }
        }
    }
    */

    std::cout << "Wolff run with: L = " << L << ", T = " << T << "\n";

    std::cout << "(start nonlinear relaxation run...)";

    //nonlinear relaxation analysis
    Eigen::ArrayXd corr_nonlinear_magnetization = Eigen::ArrayXd::Zero(Nsample);

    const int num_chains = 50;  //avg over 50 markov chains
    double runtime = 0.0;
    for(int k = 0; k < num_chains; ++k){
        //reset configuration
        x.fill(1);                          
        M = N;
        E = -3*J*N;
        
        auto start = std::chrono::steady_clock::now();  //also stop time 
        for(unsigned n = 0; n < Nsample; ++n){
            move(x, M, E);
            corr_nonlinear_magnetization[n] += double(M)*N_inv;
        }
        auto end = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed_seconds = end-start;
        runtime += elapsed_seconds.count();
    }
    runtime /= double(Nsample*num_chains);
    

    corr_nonlinear_magnetization /= double(num_chains);

    std::cout << "(start linear relaxation run...)\n\n";
    
    //linear relaxation analysis
    Eigen::ArrayXd corr_linear_magnetization = Eigen::ArrayXd::Zero(Nsample);

    lattice_3D<int16_t> corr_t0_conf = x;
    const double corr_t0_magn_2 = double(M*M)*N_inv*N_inv;
    const double denominator = 1.0/(1.0 - corr_t0_magn_2);

    
    for(unsigned n = 0; n < Nsample; ++n){
        move(x, M, E);
        const double phi = (x.compute_magnetization_correlation_(corr_t0_conf) - corr_t0_magn_2) * denominator;
        corr_linear_magnetization[n] = phi;
    }


    //save E(t) into file
    std::ofstream myfile;
    myfile.open ("outputs/wolff_" + std::to_string(L) + "_" + std::to_string(T) + ".csv");
    myfile << L << "," << T << "," << runtime << "\n";
    myfile << corr_nonlinear_magnetization.transpose().format(CSVFormat) << "\n";
    myfile << corr_linear_magnetization.transpose().format(CSVFormat) << "\n";
    myfile.close();

    return 0;
}