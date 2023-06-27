#include<iostream>
#include<fstream>
#include<cmath>
#include<random>
#include<Eigen/Dense>
#include"lattice.hpp"

//global defined constants
int L, N;
double N_inv;
const int    J = 1;               //coupling constant
//const double H = 0.0;
//const double kB = 1.0;              //Boltzmann constant [m^2 kg s^-2 K^-1]

const int Nthermalization = int(10e8); //number of thermalization steps
const int Nsample = 5000;              //number of samples (= size of the Markov chain)
const int Nsubsweep = 10000;            //number of subsweeps (to generate better samples)

std::mt19937 gen(42);
std::uniform_int_distribution<> uniform;
std::uniform_real_distribution<double> uniform_real(0.0, 1.0);


/*
Args:
    x: Spin configuration
    M: Magnetization of x
    E: Energy of x

Updates:
    x, M and E after one Monte Carlo move
*/
template<typename T>
void move(lattice_3D<T>& x, int& M, int& E, const Eigen::ArrayXd& R_vec){
    // Probability look-up tables
	// TODO: optionally use probability lookup tables
	
    //pick one site at random
    const int i = uniform(gen);
    const int j = uniform(gen);
    const int k = uniform(gen);

    //Flip the spin of that site according to the Metropolis algorithm
    //Compute the local magnetic field at site (i,j) due to nearest-neighbours
    const int x_old = x(i,j,k);
    const int nn = x_old*x.nn_sum(i,j,k);
    const int dE = 2*J*nn;

    if(dE <= 0){
        //flip the spin
        x(i,j,k) *= -1;

        //update the magnetisation and energy
        M -= 2*x_old;
        E += 2*J*nn;
    }
    else{
        //Compute the Metropolis acceptance probability `R` and compare it to a random number in [0,1)
        const double R = R_vec(int(0.5*(nn+6)));
        const double eta = uniform_real(gen);

        if(R > eta){
            //flip the spin
            x(i,j,k) *= -1;

            //update the magnetisation and energy
            M -= 2*x_old;
            E += 2*J*nn;
        }
    }
}



int main(int argc, char *argv[]){
    if(argc > 3){
        printf("Give temperature T in range (0.0, 10.0) and Lattice size\n");
        return 1;
    }
    const double temp = atof(argv[1]);
    const double beta = 1.0/temp;  //NOTE: beta = 1/(kbT) with kb = 1.0, hence inverse of Temp

    L = atoi(argv[2]);
    N = L*L*L;
    N_inv = 1.0/double(N);
    uniform = std::uniform_int_distribution<>(0,L-1);

    //fill lookup table of acceptance probability
    Eigen::ArrayXd R_vec(7);
    R_vec << std::exp(-beta*2.0*J*(-6)), std::exp(-beta*2.0*J*(-4)), std::exp(-beta*2.0*J*(-2)), 1.0, std::exp(-beta*2.0*J*(2)), std::exp(-beta*2.0*J*(4)), std::exp(-beta*2.0*J*(6));

    //initialize lattice matrix
    lattice_3D<int> x(L);        //calloc (initiallized to zero)
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
        move(x, M, E, R_vec);
    }

    //measurement of M and E
    Eigen::ArrayXd M_data = Eigen::ArrayXd::Zero(Nsample);

    M_data(0) = std::abs(double(M))*N_inv;

    for(int n = 1; n < Nsample; ++n){
        for(int N_ = 0; N_ < Nsubsweep; ++N_){
            move(x, M, E, R_vec);
        }
        M_data(n) = std::abs(double(M))*N_inv;
    }

    const double M_mean = M_data.mean();
    const double M_var  = M_data.square().mean() - M_mean*M_mean;
    
    //save into file for python ploting
    //L, T, chi
    std::cout << L << "," << temp << "," << N*beta*M_var << "\n";
    
    return 0;
}