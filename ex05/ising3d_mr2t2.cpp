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

const int J = 1;               //coupling constant
double beta;                    

int Nthermalization; //number of thermalization steps
int Nsample;              //number of samples (= size of the Markov chain)
int Nsubsweep;            //number of subsweeps (to generate better samples)

std::mt19937 gen(42);
std::uniform_int_distribution<> uniform;
std::uniform_real_distribution<double> uniform_real(0.0, 1.0);
std::normal_distribution<double> normal(0.0, 1.0);



/*
Args:
    x: Spin configuration
    M: Magnetization of x
    E: Energy of x

Updates:
    x, M and E after one Monte Carlo move
*/
template<typename T>
void move(lattice_3D<T>& x, Eigen::Vector3d& M, double& E){
	
    //pick one site at random
    const int i = uniform(gen);
    const int j = uniform(gen);
    const int k = uniform(gen);

    //generate a random flip coordinate
    Eigen::Vector3d s_flip(normal(gen), normal(gen), normal(gen));
    s_flip.normalize();


    //Flip the spin of that site according to the Metropolis algorithm
    //Compute the local magnetic field at site (i,j) due to nearest-neighbours
    const Eigen::Vector3d x_old = x(i,j,k);
    const double dE = (s_flip - x_old).dot(x.nn_sum(i,j,k));

    if(dE < 0.0){
        //flip the spin
        x(i,j,k) = s_flip;

        //update the magnetisation and energy
        M += s_flip - x_old;
        E += dE;
    }
    else{
        //Compute the Metropolis acceptance probability `R` and compare it to a random number in [0,1)
        const double R = std::exp(-beta*dE);
        const double eta = uniform_real(gen);

        if(R > eta){
            //flip the spin
            x(i,j,k) = s_flip;

            //update the magnetisation and energy
            M += s_flip - x_old;
            E += dE;
        }
    }
}



int main(int argc, char *argv[]){
    if(argc != 3){
        printf("Give temperature T in range (0.0, 10.0), system size L\n");
        return 1;
    }
    
    const double T = atof(argv[1]);
    beta = 1.0/T;  //NOTE: beta = 1/(kbT) with kb = 1.0, hence inverse of Temp

    L = atoi(argv[2]);
    N = L*L*L;
    N_inv = 1.0/double(N);
    N2_inv = N_inv*N_inv;
    uniform = std::uniform_int_distribution<>(0,L-1);
    Nthermalization = 30*N;
    Nsample = 4000;
    Nsubsweep = 10*N;


    //initialize lattice matrix
    lattice_3D<Eigen::Vector3d> x(L);     //calloc (initiallized to zero)
    const Eigen::Vector3d tmp(1.0, 0.0, 0.0);
    //x.fill(Eigen::Vector3d::Constant(1.0/std::sqrt(3)));      //fill
    x.fill(tmp);
    Eigen::Vector3d M(double(N), 0.0, 0.0);
    double E = -3*J*N;

    /*
    for(int i = 0; i < L; ++i){
        for(int j = 0; j < L; ++j){
            for(int k = 0; k < L; ++k){
                Eigen::Vector3d s_flip(normal(gen), normal(gen), normal(gen));
                s_flip.normalize();

                //update Magnetization
                M += s_flip - tmp;
                
                //update total Energy
                const double dE = (s_flip - tmp).dot(x.nn_sum(i,j,k));
                E += dE;


                x(i,j,k) = s_flip;
            }
        }
    }
    */

    //thermalisation loop
    for(int n = 0; n < Nthermalization; ++n){
        move(x, M, E);
    }

    //measurement of M and E

    Eigen::ArrayXd M_data = Eigen::ArrayXd::Zero(Nsample);
    Eigen::ArrayXd E_data = Eigen::ArrayXd::Zero(Nsample);

    M_data(0) = M.norm()*N_inv;
    E_data(0) = E*N_inv;

    for(int n = 1; n < Nsample; ++n){
        for(int N_ = 0; N_ < Nsubsweep; ++N_){
            move(x, M, E);
        }
        M_data(n) = M.norm()*N_inv;
        E_data(n) = E*N_inv;
    }


    const double M_4 = M_data.square().square().mean();
    const double M_2 = M_data.square().mean();

    const double U_L = 1.0 - M_4/(3.0*M_2*M_2);

    const double M_mean = M_data.mean();
    const double M_var  = M_2 - M_mean*M_mean;
    const double E_mean = E_data.mean();
    const double E_var  = E_data.square().mean() - E_mean*E_mean;


    //save into file for python ploting
    //L, T, U_L, <M>, chi, <E>, Cv
    std::cout << L << "," << T << "," << U_L << "," << M_mean << "," << N*beta*M_var << "," << E_mean << "," << beta*beta*E_var << "\n";


    return 0;
}