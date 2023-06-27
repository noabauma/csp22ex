#include<iostream>
#include<fstream>
#include<cmath>
#include<random>
#include<Eigen/Dense>
#include"lattice.hpp"

const static Eigen::IOFormat CSVFormat(Eigen::StreamPrecision, Eigen::DontAlignCols, ",", "\n");

//global defined constants
const int L = 16;
const int N = L*L*L;
const double N_inv = 1.0/double(N);
const int    J = 1;               //coupling constant
//const double H = 0.0;
//const double kB = 1.0;              //Boltzmann constant [m^2 kg s^-2 K^-1]

const int Nthermalization = int(10e8); //number of thermalization steps
const int Nsample = 5000;              //number of samples (= size of the Markov chain)
const int Nsubsweep = 1;            //number of subsweeps (to generate better samples)

std::mt19937 gen(42);
std::uniform_int_distribution<> uniform(0,L-1);
std::uniform_real_distribution<double> uniform_real(0.0, 1.0);


/*
Args:
    x: Spin configuration

Returns:
    Total energy of configuration x.
*/
template<typename T>
int total_energy(const lattice_3D<T>& x){
    int energy = 0.0;
    for(int i = 0; i < L; ++i){
        for(int j = 0; j < L; ++j){
            for(int k = 0; k < L; ++k){
                energy += x(i,j,k)*x.nn_sum(i,j,k);
            }
        }
    }
    energy *= -J*0.5;   //is it 0.5 or 1/3 due to 3dim?
			
    return energy;  //Note J is an int due to performance. If J is no longer an int. change this formula
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
    if(argc > 2){
        printf("Give temperature T in range (0.0, 10.0)\n");
        return 1;
    }
    const double temp = atof(argv[1]);
    const double beta = 1.0/temp;  //NOTE: beta = 1/(kbT) with kb = 1.0, hence inverse of Temp

    //fill lookup table of acceptance probability
    Eigen::ArrayXd R_vec(7);
    R_vec << std::exp(-beta*2.0*J*(-6)), std::exp(-beta*2.0*J*(-4)), std::exp(-beta*2.0*J*(-2)), 1.0, std::exp(-beta*2.0*J*(2)), std::exp(-beta*2.0*J*(4)), std::exp(-beta*2.0*J*(6));

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
        move(x, M, E, R_vec);
    }

    //measurement of M and E

    //printf("Sampling M and E ...\n");

    Eigen::ArrayXd avg_spins   = Eigen::ArrayXd::Zero(Nsample);  //Task1
    Eigen::ArrayXd avg_spins_0 = Eigen::ArrayXd::Zero(Nsample);  //Task1

    avg_spins(0)   = x.sum()*N_inv;
    avg_spins_0(0) = 1;

    lattice_3D<int> x_0 = x; 

    for(int n = 1; n < Nsample; ++n){
        for(int N_ = 0; N_ < Nsubsweep; ++N_){  //still good to make bigger leaps
            move(x, M, E, R_vec);
        }
        avg_spins(n)   = x.sum()*N_inv;

        unsigned count = 0;
        for(int i = 0; i < N; ++i){
            count += x.sum(x_0(i));
        }
        avg_spins_0(n) = count*N_inv*N_inv;
    }

    const double spins_mean_0    = avg_spins(0);
    const double spins_mean_inf  = avg_spins(Nsample-1);

    const Eigen::ArrayXd Phi    = (avg_spins - spins_mean_inf)/(spins_mean_0 - spins_mean_inf);
    const Eigen::ArrayXd Phi_nl = (avg_spins_0 - spins_mean_0*spins_mean_0)/(1.0 - spins_mean_0*spins_mean_0);

    //save into file
    std::ofstream myfile;
    myfile.open ("task1.csv");
    myfile << Phi.transpose().format(CSVFormat) << "\n" << Phi_nl.transpose().format(CSVFormat) << "\n";
    myfile.close();
    
    return 0;
}