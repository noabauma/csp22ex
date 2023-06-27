#include<iostream>
#include<cmath>
#include<random>
#include<Eigen/Dense>
#include"lattice.hpp"

//global defined constants
int L = 16;
int N = L*L*L;
double N_inv = 1.0/double(N);
const int    J = 1;               //coupling constant

const int Nsample = 5000;              //number of samples (= size of the Markov chain)
const int Nsubsweep = 3*N;             //number of subsweeps (to generate better samples)

const int Emax = 40;    //maximum energy difference 
int Ed   = 4;           //intial demon energy

std::mt19937 gen(42);
std::uniform_int_distribution<> uniform(0,L-1);

/*
Args:
    x: Spin configuration
    M: Magnetization of x
    E: Energy of x

Updates:
    x, M and E after one Monte Carlo move
*/
template<typename T>
void increase_energy(lattice_3D<T>& x, int& M, int& E){
	
    //pick one site at random
    const int i = uniform(gen);
    const int j = uniform(gen);
    const int k = uniform(gen);

    //flip the spin at site (i,j,k) if and only if the energy increases
    //update magnetisation and energy after the flip
    const int x_old = x(i,j,k);
    const int nn = x_old*x.nn_sum(i,j,k);
    const int dE = 2*J*nn;

    if(dE > 0.0){
        //flip the spin
        x(i,j,k) *= -1;

        //update the magnetisation and energy
        M -= 2*x_old;
        E += dE;
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
void move(lattice_3D<T>& x, int& M, int& E, int& Ed){
    //pick one site at random
    const int i = uniform(gen);
    const int j = uniform(gen);
    const int k = uniform(gen);

    //flip the spin according to the Creutz algorithm
    //compute energy difference dE due to flipping spin at site (i,j,k)
    //perform the Creutz step by comparing `dE` to the demon energy `Ed`
    //update `M`, `E` and also `Ed` if the spin flip is accepted
    const int x_old = x(i,j,k);
    const int nn = x_old*x.nn_sum(i,j,k);
    const int dE = 2*J*nn;

    if(Ed - dE >= 0 && Emax >= Ed - dE){
        //flip the spin
        x(i,j,k) *= -1;

        //update the magnetisation and energy
        M  -= 2*x_old;
        E  += dE;
        Ed -= dE;
    }
}



int main(int argc, char *argv[]){
    int Esys;

    if(argc == 2){
        Esys  = atoi(argv[1]);
    }else if(argc == 3){
        Esys  = atoi(argv[1]);
        L = atoi(argv[2]);
        N = L*L*L;
        N_inv = 1.0/double(N);
        uniform = std::uniform_int_distribution<>(0,L-1);
    }else{
        printf("Give system energies E in range [-3N, 0]\n");
        return 1;
    }
    
    
    
    //initialize lattice matrix
    lattice_3D<int8_t> x(L);        //calloc (initiallized to zero)
    x.increment();               //element-wise ++x
    int M = N;
    int E = -3*J*N;

    while(Esys > E){
        increase_energy(x, M, E);
    }
    

    Eigen::ArrayXd M_data = Eigen::ArrayXd::Zero(Nsample);
    Eigen::ArrayXd E_data = Eigen::ArrayXd::Zero(Nsample);
    Eigen::ArrayXd Ed_arr = Eigen::ArrayXd::Zero((Nsample-1)*Nsubsweep);

    M_data(0) = std::abs(double(M))*N_inv;
    E_data(0) = double(E)*N_inv;

    for(int n = 1; n < Nsample; ++n){
        for(int t = 0; t < Nsubsweep; ++t){
            move(x, M, E, Ed);
            Ed_arr((n-1)*Nsubsweep+t) = Ed;
        }
        M_data(n) = std::abs(double(M))*N_inv;
        E_data(n) = double(E)*N_inv;
    }

    const double M_mean = M_data.mean();
    const double E_mean = E_data.mean();
    const double M_var  = M_data.square().mean() - M_mean*M_mean;
    const double E_var  = E_data.square().mean() - E_mean*E_mean;

    //estimate the temperature of the system from the statistics of demon energies using data `Ed_arr`
    const double beta = 0.25 * std::log(1.0 + 4.0/Ed_arr.mean());
    const double T = 1.0/beta;

    //save into file for python ploting
    //L, Esys, T, M, M_std, chi, E, E_std, Cv
    std::cout << L << "," << Esys << "," << T << "," << M_mean << "," << std::sqrt(M_var) << "," << N*beta*M_var << "," << E_mean << "," << std::sqrt(E_var) << "," << beta*beta*E_var << "\n";
    
    return 0;
}