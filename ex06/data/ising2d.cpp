#include<iostream>
#include<fstream>
#include<cmath>
#include<random>
#include<Eigen/Dense>

const static Eigen::IOFormat CSVFormat(Eigen::StreamPrecision, Eigen::DontAlignCols, ",", "\n");

//global defined constants
const int L = 32;
const int N = L*L;
const double N_inv = 1.0/double(N);
const int    J = 1;               //coupling constant

const int Nthermalization = 30*N;      //number of thermalization steps
const int Nsample = 4000;              //number of samples (= size of the Markov chain)
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
    M: Magnetization of x
    E: Energy of x

Updates:
    x, M and E after one Monte Carlo move
*/
void move(Eigen::MatrixXi& x, int& M, int& E, const Eigen::VectorXd& R_vec){
	
    //pick one site at random
    const int i = uniform(gen);
    const int j = uniform(gen);

    //Flip the spin of that site according to the Metropolis algorithm
    //Compute the local magnetic field at site (i,j) due to nearest-neighbours
    const int x_old = x(i,j);
    const int nn = x_old*nn_sum(x,i,j);
    const int dE = 2*J*nn;

    if(dE <= 0){
        //flip the spin
        x(i,j) *= -1;

        //update the magnetisation and energy
        M -= 2*x_old;
        E += dE;
    }
    else{
        //Compute the Metropolis acceptance probability `R` and compare it to a random number in [0,1)
        const double R = R_vec(int(0.5*(nn+4)));
        const double eta = uniform_real(gen);

        if(R > eta){
            //flip the spin
            x(i,j) *= -1;

            //update the magnetisation and energy
            M -= 2*x_old;
            E += dE;
        }
    }
}



int main(int argc, char *argv[]){
    if(argc != 3){
        printf("Give temperature T in range (0.0, 10.0) & if train_data\n");
        return 1;
    }
    const double T = atof(argv[1]);
    const double beta = 1.0/T;  //NOTE: beta = 1/(kbT) with kb = 1.0, hence inverse of Temp

    const bool make_train_data = bool(atoi(argv[2]));
    const std::string folder = make_train_data ? "train_data/" : "test_data/";

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
    Eigen::MatrixXi container(Nsample, N);
    container.row(0) = Eigen::Map<Eigen::RowVectorXi>(x.data(), x.size());

    for(int n = 1; n < Nsample; ++n){
        for(int N_ = 0; N_ < Nsubsweep; ++N_){
            move(x, M, E, R_vec);
        }
        container.row(n) = Eigen::Map<Eigen::RowVectorXi>(x.data(), x.size());
    }

    //save into file for python ploting
    std::ofstream myfile;
    myfile.open (folder + std::to_string(T) + ".csv");
    myfile << T << "\n";
    myfile << container.format(CSVFormat) << "\n";
    myfile.close();


    return 0;
}