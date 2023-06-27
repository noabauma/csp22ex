#include<iostream>
#include<cmath>
#include<random>
#include<Eigen/Dense>

//THIS is version 2! Improvements and correcting mistakes (hopefully)


//global defined constants
const int L = 10;
const int N = L*L;
const double N_inv = 1.0/double(N);
const int    J = 1;               //coupling constant
//const double H = 0.0;
//const double kB = 1.0;              //Boltzmann constant [m^2 kg s^-2 K^-1]

const int Nthermalization = int(10e5); //number of thermalization steps
const int Nsample = 5000;              //number of samples (= size of the Markov chain)
const int Nsubsweep = 10*N;            //number of subsweeps (to generate better samples)

std::mt19937 gen(42);
std::uniform_int_distribution<> uniform(0,L-1);
std::uniform_real_distribution<double> uniform_real(0.0, 1.0);

/*
Args:
    x: Spin configuration
    i, j: Indices describing the position of one spin

Returns:
    Sum of the spins in x which are nearest neighbors of (i, j)
*/
int nn_sum(const Eigen::MatrixXi& x, const int i, const int j){
    return x((i+1)%L, j) + x((i-1+L)%L, j) + x(i, (j+1)%L) + x(i, (j-1+L)%L);  
}


/*
Args:
    x: Spin configuration

Returns:
    Total energy of configuration x.
*/
int total_energy(const Eigen::MatrixXi& x){
    int energy = 0.0;
    for(int i = 0; i < L; ++i){
        for(int j = 0; j < L; ++j){
            energy += x(i,j)*nn_sum(x,i,j);
        }
    }
    energy *= -J*0.5;
			
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
void move(Eigen::MatrixXi& x, int& M, int& E, const Eigen::VectorXd& R_vec){
    // Probability look-up tables
	// TODO: optionally use probability lookup tables
	
    //pick one site at random
    //const int i = rand_int[2*n_];
    //const int j = rand_int[2*n_+1];
    const int i = uniform(gen);
    const int j = uniform(gen);

    //Flip the spin of that site according to the Metropolis algorithm
    //Compute the local magnetic field at site (i,j) due to nearest-neighbours
    const int x_old = x(i,j);
    const int nn = x_old*nn_sum(x,i,j);
    const int dE = 2*J*nn;

    if(dE <= 0){
        //flip the spin
        //E += J*(nn_sum(x,i,j) + nn_sum(x,(i+1)%L,j) + nn_sum(x,(i-1+L)%L,j) + nn_sum(x,i,(j+1)%L) + nn_sum(x,i,(j-1+L)%L));
        x(i,j) *= -1;

        //update the magnetisation and energy
        M -= 2*x_old;
        E += dE;    // TODO: this is weird and looks wrong pls check
        //E += -J*(nn_sum(x,i,j) + nn_sum(x,(i+1)%L,j) + nn_sum(x,(i-1+L)%L,j) + nn_sum(x,i,(j+1)%L) + nn_sum(x,i,(j-1+L)%L));
    }
    else{
        //Compute the Metropolis acceptance probability `R` and compare it to a random number in [0,1)
        const double R = R_vec(int(0.5*(nn+4)));
        const double eta = uniform_real(gen);

        if(R > eta){
            //flip the spin
            //E += J*(nn_sum(x,i,j) + nn_sum(x,(i+1)%L,j) + nn_sum(x,(i-1+L)%L,j) + nn_sum(x,i,(j+1)%L) + nn_sum(x,i,(j-1+L)%L));
            x(i,j) *= -1;

            //update the magnetisation and energy
            M -= 2*x_old;
            E += dE;    // TODO: this is weird and looks wrong pls check
            //E += -J*(nn_sum(x,i,j) + nn_sum(x,(i+1)%L,j) + nn_sum(x,(i-1+L)%L,j) + nn_sum(x,i,(j+1)%L) + nn_sum(x,i,(j-1+L)%L));
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
    Eigen::VectorXd R_vec(5);
    R_vec << std::exp(-beta*2.0*J*(-4)), std::exp(-beta*2.0*J*(-2)), 1.0, std::exp(-beta*2.0*J*(2)), std::exp(-beta*2.0*J*(4));

    //initialize lattice matrix
    Eigen::MatrixXi x = Eigen::MatrixXi::Ones(L,L);
    int M = N;
    int E = -2*J*N;


    for(int i = 0; i < L; ++i){
        for(int j = 0; j < L; ++j){
            if(uniform_real(gen) < 0.5){
                x(i,j) = -1;
                M -= 2;
                E += 2*J*nn_sum(x, i, j);
            }
        }
    }
    

    //thermalisation loop
    for(int n = 0; n < Nthermalization; ++n){
        move(x, M, E, R_vec);
    }

    //measurement of M and E

    //printf("Sampling M and E ...\n");

    Eigen::ArrayXd M_data = Eigen::ArrayXd::Zero(Nsample);
    Eigen::ArrayXd E_data = Eigen::ArrayXd::Zero(Nsample);

    M_data(0) = std::abs(double(M))*N_inv;
    E_data(0) = double(E)*N_inv;            //if averaging energy wanted yes

    for(int n = 1; n < Nsample; ++n){
        for(int N_ = 0; N_ < Nsubsweep; ++N_){
            move(x, M, E, R_vec);
        }
        M_data(n) = std::abs(double(M))*N_inv;
        E_data(n) = double(E)*N_inv;
    }

    const double M_mean = M_data.mean();
    const double E_mean = E_data.mean();
    const double M_std  = M_data.square().mean() - M_mean*M_mean;
    const double E_std  = E_data.square().mean() - E_mean*E_mean;

    //save into file for python ploting
    //T, M, M_std, chi, E, E_std, Cv
    std::cout << temp << "," << M_mean << "," << std::sqrt(M_std) << "," << N*beta*M_std << "," << E_mean << "," << std::sqrt(E_std) << "," << beta*beta*E_std << "\n";

    return 0;
}