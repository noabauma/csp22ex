#include<iostream>
#include<cmath>
#include<random>
#include<Eigen/Dense>
#include"lattice.hpp"

//global defined constants
int L;
int N;
double N_inv;
double N2_inv;

const int J = 1;            //coupling constant
double p;                   //prob to add to cluster

const int Nthermalization = int(10e5); //number of thermalization steps
const int Nsample = 5000;              //number of samples (= size of the Markov chain)
int Nsubsweep;            //number of subsweeps (to generate better samples)

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
    const double beta = 1.0/T;  //NOTE: beta = 1/(kbT) with kb = 1.0, hence inverse of Temp
    p = 1.0 - std::exp(-2.0*beta*J);

    L = atoi(argv[2]);
    N = L*L*L;
    N_inv = 1.0/double(N);
    N2_inv = N_inv*N_inv;
    uniform = std::uniform_int_distribution<>(0,L-1);
    Nsubsweep = 10*N;

    //initialize lattice matrix
    lattice_3D<int> x(L);       //calloc (initiallized to zero)
    x.increment();               //element-wise ++x
    int M = N;
    int E = -3*J*N;


    
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


    //thermalisation loop
    for(int n = 0; n < Nthermalization; ++n){
        move(x, M, E);
    }
    
    
    
    //move(x, M, E);

    //std::cout << x;

    //measurement of M and E

    Eigen::ArrayXd M_data = Eigen::ArrayXd::Zero(Nsample);
    Eigen::ArrayXd E_data = Eigen::ArrayXd::Zero(Nsample);

    M_data(0) = std::abs(double(M))*N_inv;
    E_data(0) = double(E)*N_inv;

    for(int n = 1; n < Nsample; ++n){
        for(int N_ = 0; N_ < Nsubsweep; ++N_){
            move(x, M, E);
        }
        M_data(n) = std::abs(double(M))*N_inv;
        E_data(n) = double(E)*N_inv;
    }


    const double M_4 = M_data.square().square().mean();
    const double M_2 = M_data.square().mean();

    const double U_L = 1.0 - M_4/(3.0*M_2*M_2);

    //const double M_mean = M_data.mean();
    //const double M_var  = M_data.square().mean() - M_mean*M_mean;


    //save into file for python ploting
    //L, T, U_L, chi
    std::cout << L << "," << T << "," << U_L << "\n";

    return 0;
}