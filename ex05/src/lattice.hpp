#include<iostream>
#include<Eigen/Dense>


//This file consists of a class for constructing a 3D lattice

// For typename = Eigen::Vector3d!!

template<typename T>
class lattice_3D{
    public:
        //normal constructor
        lattice_3D(const unsigned L_){
            L = L_;
            LL = L*L;
            LLL = LL*L;
            N_inv = 1.0/double(LLL);
            lattice = (T*) calloc(LLL, sizeof(T));
        }
        
        //copy constructor
        lattice_3D(const lattice_3D<T>& x){
            L = x.L;
            LL = L*L;
            LLL = LL*L;
            N_inv = 1.0/double(LLL);
            lattice = (T*) calloc(LLL, sizeof(T));
            memcpy(lattice, x.lattice, LLL*sizeof(T));
        }

        //desctructor
        ~lattice_3D(){
            free(lattice);
        }

        void fill(const T val) const{
            for(unsigned i = 0; i < LLL; ++i){
                lattice[i] = val;
            }
        }

        void increment() const{
            for(unsigned i = 0; i < LLL; ++i){
                ++lattice[i];
            }
        }

        T& operator()(const unsigned i, const unsigned j, const unsigned k) const {
            return lattice[i*LL + j*L + k];
        }

        T& operator[](const unsigned i) const {
            return lattice[i];
        }

        
        T nn_sum(const unsigned i, const unsigned j, const unsigned k){
            return lattice[i*LL + j*L + (k+1)%L] + lattice[i*LL + j*L + (k-1+L)%L] + lattice[i*LL + (j+1)%L *L + k] + lattice[i*LL + (j-1+L)%L *L + k] + lattice[(i+1)%L *LL + j*L + k] + lattice[(i-1+L)%L *LL + j*L + k];  
        }

        T sum(){
            T count = T::Zero();
            for(unsigned i = 0; i < LLL; ++i){
                count += lattice[i];
            }
            return count;
        }

        double compute_magnetization_correlation_(const lattice_3D<T> &y){
            T result = T::Zero();
            for(unsigned i = 0; i < LLL; ++i){
                result += lattice[i]*(y[i]);
            }
            return double(result) * N_inv;
        }       

        unsigned shape(){
            return L;
        }

        
    private:
        unsigned L, LL, LLL;
        double N_inv;
        T *lattice;
};

template<typename T>
std::ostream& operator<<(std::ostream& os, lattice_3D<T>& x){
    const unsigned L = x.shape();

    for(unsigned i = 0; i < L; ++i){
        for(unsigned j = 0; j < L; ++j){
            for(unsigned k = 0; k < L; ++k){
                const T s = x(i, j, k);
                const std::string str = " "; 
                os << str << s.transpose();
            }
            os << "\n";
        }
        os << "\n";
    }
    return os;
}

/*
int main(){
    lattice_3D<Eigen::Vector3d> x(2);
    x.fill(Eigen::Vector3d::Constant(0.5));
    x(0,0,1) *= 3.0;
    x(0,0,1)(1) = 5.0;
    std::cout << "x =\n" << x << "\n";
    std::cout << "1. = " << x(0,0,0).transpose()*x(0,0,1) << "\n";
    std::cout << "2. = " << x(0,0,0).dot(x(0,0,1)) << "\n";
    std::cout << "3. = " << x(0,0,0).cross(x(0,0,1)) << "\n";
    std::cout << "4. = " << x(0,0,0).adjoint()*x(0,0,1) << "\n";

    std::cout << sizeof(Eigen::Vector3d) << "\n";
    std::cout << sizeof(double) << "\n";
}
*/