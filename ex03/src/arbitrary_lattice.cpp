#include<iostream>
#include<tuple>

//This file consists of a class for constructing a 3D lattice

template<typename T>
class lattice_3D{
    public:
        lattice_3D(const unsigned N, const unsigned M, const unsigned O){
            N_ = N;
            M_ = M;
            O_ = O;
            NMO_ = N*M*O;
            lattice = (T*) calloc(N_ * M_ * O_, sizeof(T));
        }

        void fill(const T val) const{
            for(unsigned i = 0; i < NMO_; ++i){
                lattice[i] = val;
            }
        }

        void increment() const{
            for(unsigned i = 0; i < N_*M_*O_; ++i){
                ++lattice[i];
            }
        }

        T operator()(const unsigned i, const unsigned j, const unsigned k){
            return lattice[i*M_*O_ + j*O_ + k];
        }

        T& operator=(const T& value){
            return *this = value;
        }

        
        T nn_sum(const unsigned i, const unsigned j, const unsigned k){
            const T s = lattice[i*M_*O_ + j*O_ + k];

            return s*(lattice[i*M_*O_ + j*O_ + (k+1)%O_]);  
        }

        T sum(){
            T count = 0;
            for(unsigned i = 0; i < NMO_; ++i){
                count += lattice[i];
            }
        }
        

        std::tuple<unsigned, unsigned, unsigned> shape(){
            return std::make_tuple(N_, M_, O_);
        }

        
    private:
        unsigned N_, M_, O_, NMO_;
        T *lattice;

};

template<typename T>
std::ostream& operator<<(std::ostream& os, lattice_3D<T>& lattice){
    std::tuple<unsigned, unsigned, unsigned> shapes = lattice.shape();
    const unsigned N = std::get<0>(shapes);
    const unsigned M = std::get<1>(shapes);
    const unsigned O = std::get<2>(shapes);

    for(unsigned k = 0; k < O; ++k){
        for(unsigned i = 0; i < N; ++i){
            for(unsigned j = 0; j < M; ++j){
                os << lattice(i, j, k) << " ";
            }
            os << "\n";
        }
        os << "\n";
    }
    return os;
}


int main(){
    lattice_3D<int> x(2,2,3);
    x.increment();
    x(1,1,1) = 3;
    std::cout << x << "\n";
}