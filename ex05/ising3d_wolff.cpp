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

Returns:
    Total energy of configuration x.
*/
template<typename T>
double total_energy(lattice_3D<T>& x){
    double energy = 0.0;
    for(int i = 0; i < L; ++i){
        for(int j = 0; j < L; ++j){
            for(int k = 0; k < L; ++k){
                energy += x(i,j,k).dot(x.nn_sum(i,j,k));
            }
        }
    }
    energy *= -J*0.5;
			
    return energy;  //Note J is an int due to performance. If J is no longer an int. change this formula
}



/*
Args:

Updates:
    Recursivly visits sites and updates them if they are part of the cluster
*/
template<typename T>
void visit(lattice_3D<T>& x, Eigen::Vector3d& M, double& E, lattice_3D<bool>& visited, const Eigen::Vector3d n_hat, int i, int j, int k){
    
    visited(i,j,k) = true; //site marked as visited (being part of the cluster!)

    const Eigen::Vector3d x_old = x(i,j,k);

    x(i,j,k) -= 2.0*(x_old.dot(n_hat))*n_hat; //flip spin

    M += x(i,j,k) - x_old;
    //E += (x(i,j,k) - x_old).dot(x.nn_sum(i,j,k)); //TODO: pls check

    //now that the top site is part of the cluster we go visit other nn sites
    const int i_p = (i+1)%L;
    const int i_n = (i-1+L)%L;
    const int j_p = (j+1)%L;
    const int j_n = (j-1+L)%L;
    const int k_p = (k+1)%L;
    const int k_n = (k-1+L)%L;


    if(!visited(i_p,j,k)){
        //const Eigen::Vector3d R_Sj = x(i_p,j,k) - 2.0*n_hat*x(i_p,j,k).dot(n_hat);
        //const double dE_i_p = J*n_hat.dot(R_Sj)*n_hat.dot(x(i,j,k));
        const double dE_i_p = (x(i,j,k).dot(n_hat)) * (x(i_p,j,k).dot(n_hat));
        if(dE_i_p < 0.0){
            const double p_i_p = 1.0 - std::exp(2.0*beta*dE_i_p);
            if(uniform_real(gen) < p_i_p){
                visit(x, M, E, visited, n_hat, i_p, j, k);
            }
        }
    }

    if(!visited(i_n,j,k)){
        //const Eigen::Vector3d R_Sj = x(i_n,j,k) - 2.0*n_hat*x(i_n,j,k).dot(n_hat);
        //const double dE_i_n = J*n_hat.dot(R_Sj)*n_hat.dot(x(i,j,k));
        const double dE_i_n = (x(i,j,k).dot(n_hat)) * (x(i_n,j,k).dot(n_hat));
        if(dE_i_n < 0.0){
            const double p_i_n = 1.0 - std::exp(2.0*beta*dE_i_n);
            if(uniform_real(gen) < p_i_n){
                visit(x, M, E, visited, n_hat, i_n, j, k);
            }
        }
    }

    if(!visited(i,j_p,k)){
        //const Eigen::Vector3d R_Sj = x(i,j_p,k) - 2.0*n_hat*x(i,j_p,k).dot(n_hat);
        //const double dE_j_p = J*n_hat.dot(R_Sj)*n_hat.dot(x(i,j,k));
        const double dE_j_p = (x(i,j,k).dot(n_hat)) * (x(i,j_p,k).dot(n_hat));
        if(dE_j_p < 0.0){
            const double p_j_p = 1.0 - std::exp(2.0*beta*dE_j_p);
            if(uniform_real(gen) < p_j_p){
                visit(x, M, E, visited, n_hat, i, j_p, k);
            }
        }
    }

    if(!visited(i,j_n,k)){
        //const Eigen::Vector3d R_Sj = x(i,j_n,k) - 2.0*n_hat*x(i,j_n,k).dot(n_hat);
        //const double dE_j_n = J*n_hat.dot(R_Sj)*n_hat.dot(x(i,j,k));
        const double dE_j_n = (x(i,j,k).dot(n_hat)) * (x(i,j_n,k).dot(n_hat));
        if(dE_j_n < 0.0){
            const double p_j_n = 1.0 - std::exp(2.0*beta*dE_j_n);
            if(uniform_real(gen) < p_j_n){
                visit(x, M, E, visited, n_hat, i, j_n, k);
            }
        }
    }

    if(!visited(i,j,k_p)){
        //const Eigen::Vector3d R_Sj = x(i,j,k_p) - 2.0*n_hat*x(i,j,k_p).dot(n_hat);
        //const double dE_k_p = J*n_hat.dot(R_Sj)*n_hat.dot(x(i,j,k));
        const double dE_k_p = (x(i,j,k).dot(n_hat)) * (x(i,j,k_p).dot(n_hat));
        if(dE_k_p < 0.0){
            const double p_k_p = 1.0 - std::exp(2.0*beta*dE_k_p);
            if(uniform_real(gen) < p_k_p){
                visit(x, M, E, visited, n_hat, i, j, k_p);
            }
        }
    }

    if(!visited(i,j,k_n)){
        //const Eigen::Vector3d R_Sj = x(i,j,k_n) - 2.0*n_hat*x(i,j,k_n).dot(n_hat);
        //const double dE_k_n = J*n_hat.dot(R_Sj)*n_hat.dot(x(i,j,k));
        const double dE_k_n = (x(i,j,k).dot(n_hat)) * (x(i,j,k_n).dot(n_hat));
        if(dE_k_n < 0.0){
            const double p_k_n = 1.0 - std::exp(2.0*beta*dE_k_n);
            if(uniform_real(gen) < p_k_n){
                visit(x, M, E, visited, n_hat, i, j, k_n);
            }
        }
    }
}

/*
Args:
    x: Spin configuration
    M: Magnetization of x
    E: Energy of x

Updates:
    x, M and E after one Wolff move
*/
template<typename T>
void move(lattice_3D<T>& x, Eigen::Vector3d& M, double& E){
    lattice_3D<bool> visited(L);    //all automatically false (pls check)

    //pick one site at random
    const int i = uniform(gen);
    const int j = uniform(gen);
    const int k = uniform(gen);

    //generate a random flip coordinate
    Eigen::Vector3d n_hat(normal(gen), normal(gen), normal(gen));
    n_hat.normalize();

    visited(i,j,k) = true;

    //TODO: I think the first one have to be updated
    const Eigen::Vector3d x_old = x(i,j,k);

    x(i,j,k) -= 2.0*(x_old.dot(n_hat))*n_hat; //flip with 100% prob (because it is the first)

    M += x(i,j,k) - x_old;
    E += (x(i,j,k) - x_old).dot(x.nn_sum(i,j,k)); //TODO: pls check

    const int i_p = (i+1)%L;
    const int i_n = (i-1+L)%L;
    const int j_p = (j+1)%L;
    const int j_n = (j-1+L)%L;
    const int k_p = (k+1)%L;
    const int k_n = (k-1+L)%L;

    const double tmp = x(i,j,k).dot(n_hat);

    const double dE_i_p = tmp * (x(i_p,j,k).dot(n_hat));
    const double dE_i_n = tmp * (x(i_n,j,k).dot(n_hat));
    const double dE_j_p = tmp * (x(i,j_p,k).dot(n_hat));
    const double dE_j_n = tmp * (x(i,j_n,k).dot(n_hat));
    const double dE_k_p = tmp * (x(i,j,k_p).dot(n_hat));
    const double dE_k_n = tmp * (x(i,j,k_n).dot(n_hat));

    // const double dE_i_p = tmp * J*n_hat.dot(x(i_p,j,k) - 2.0*n_hat*x(i_p,j,k).dot(n_hat));
    // const double dE_i_n = tmp * J*n_hat.dot(x(i_n,j,k) - 2.0*n_hat*x(i_n,j,k).dot(n_hat));
    // const double dE_j_p = tmp * J*n_hat.dot(x(i,j_p,k) - 2.0*n_hat*x(i,j_p,k).dot(n_hat));
    // const double dE_j_n = tmp * J*n_hat.dot(x(i,j_n,k) - 2.0*n_hat*x(i,j_n,k).dot(n_hat));
    // const double dE_k_p = tmp * J*n_hat.dot(x(i,j,k_p) - 2.0*n_hat*x(i,j,k_p).dot(n_hat));
    // const double dE_k_n = tmp * J*n_hat.dot(x(i,j,k_n) - 2.0*n_hat*x(i,j,k_n).dot(n_hat));

    if(dE_i_p < 0.0){
        const double p_i_p = 1.0 - std::exp(2.0*beta*dE_i_p);
        if(uniform_real(gen) < p_i_p){
            visit(x, M, E, visited, n_hat, i_p, j, k);
        }
    }

    if(dE_i_n < 0.0){
        const double p_i_n = 1.0 - std::exp(2.0*beta*dE_i_n);
        if(uniform_real(gen) < p_i_n){
            visit(x, M, E, visited, n_hat, i_n, j, k);
        }
    }

    if(dE_j_p < 0.0){
        const double p_j_p = 1.0 - std::exp(2.0*beta*dE_j_p);
        if(uniform_real(gen) < p_j_p){
            visit(x, M, E, visited, n_hat, i, j_p, k);
        }
    }

    if(dE_j_n < 0.0){
        const double p_j_n = 1.0 - std::exp(2.0*beta*dE_j_n);
        if(uniform_real(gen) < p_j_n){
            visit(x, M, E, visited, n_hat, i, j_n, k);
        }
    }

    if(dE_k_p < 0.0){
        const double p_k_p = 1.0 - std::exp(2.0*beta*dE_k_p);
        if(uniform_real(gen) < p_k_p){
            visit(x, M, E, visited, n_hat, i, j, k_p);
        }
    }

    if(dE_k_n < 0.0){
        const double p_k_n = 1.0 - std::exp(2.0*beta*dE_k_n);
        if(uniform_real(gen) < p_k_n){
            visit(x, M, E, visited, n_hat, i, j, k_n);
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

    
    for(int i = 0; i < L; ++i){
        for(int j = 0; j < L; ++j){
            for(int k = 0; k < L; ++k){
                Eigen::Vector3d s_flip(normal(gen), normal(gen), normal(gen));
                s_flip.normalize();

                //update Magnetization
                M += s_flip - tmp;
                
                //update total Energy
                const double dE = (s_flip - tmp).transpose()*x.nn_sum(i,j,k);
                E += dE;


                x(i,j,k) = s_flip;
            }
        }
    }
    

    //thermalisation loop
    for(int n = 0; n < Nthermalization; ++n){
        move(x, M, E);
    }

    //measurement of M and E

    Eigen::ArrayXd M_data = Eigen::ArrayXd::Zero(Nsample);
    Eigen::ArrayXd E_data = Eigen::ArrayXd::Zero(Nsample);

    M_data(0) = M.norm()*N_inv;
    E_data(0) = total_energy(x)*N_inv;

    for(int n = 1; n < Nsample; ++n){
        for(int N_ = 0; N_ < Nsubsweep; ++N_){
            move(x, M, E);
        }
        M_data(n) = M.norm()*N_inv;
        E_data(n) = total_energy(x)*N_inv;
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