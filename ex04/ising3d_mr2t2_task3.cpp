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
const int J = 1;               //coupling constant

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
    if(argc > 3 || argc < 3){
        printf("Give temperature T in range (0.0, 10.0), system size L\n");
        return 1;
    }
    
    //initialize variables
    const double T = atof(argv[1]);
    const double beta = 1.0/T;  //NOTE: beta = 1/(kbT) with kb = 1.0, hence inverse of Temp

    const unsigned L = atoi(argv[2]);
    const unsigned N = L*L*L;
    const double N_inv = 1.0/double(N);
    uniform = std::uniform_int_distribution<>(0,L-1);

    const unsigned Nsample = 20*N;              //number of samples (= size of the Markov chain)

    //fill lookup table of acceptance probability
    Eigen::ArrayXd R_vec(7);
    R_vec << std::exp(-beta*2.0*J*(-6)), std::exp(-beta*2.0*J*(-4)), std::exp(-beta*2.0*J*(-2)), 1.0, std::exp(-beta*2.0*J*(2)), std::exp(-beta*2.0*J*(4)), std::exp(-beta*2.0*J*(6));

    //initialize lattice matrix
    lattice_3D<int16_t> x(L);        //calloc (initiallized with zeros)
    x.fill(1);                   //fill with ones
    int M = N;
    int E = -3*J*N;

    /*
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

    std::cout << "MR2T2 run with: L = " << L << ", T = " << T << "\n";

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
            move(x, M, E, R_vec);
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
        move(x, M, E, R_vec);
        const double phi = (x.compute_magnetization_correlation_(corr_t0_conf) - corr_t0_magn_2) * denominator;
        corr_linear_magnetization[n] = phi;
    }


    //save E(t) into file
    std::ofstream myfile;
    myfile.open ("outputs/mr2t2_" + std::to_string(L) + "_" + std::to_string(T) + ".csv");
    myfile << L << "," << T << "," << runtime << "\n";
    myfile << corr_nonlinear_magnetization.transpose().format(CSVFormat) << "\n";
    myfile << corr_linear_magnetization.transpose().format(CSVFormat) << "\n";
    myfile.close();


    return 0;
}