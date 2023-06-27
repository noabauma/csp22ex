#include<iostream>


//This file consists of a class for constructing a 3D lattice

template<typename T>
class lattice_3D{
    public:
        //normal constructor
        lattice_3D(const unsigned L_){
            L = L_;
            LL = L*L;
            LLL = LL*L;
            lattice = (T*) calloc(LLL, sizeof(T));
        }
        
        //copy constructor
        lattice_3D(const lattice_3D<T>& x){
            L = x.L;
            LL = L*L;
            LLL = LL*L;
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

        T& operator()(const unsigned i, const unsigned j, const unsigned k){
            return lattice[i*LL + j*L + k];
        }

        T& operator()(const unsigned i){
            return lattice[i];
        }

        //Todo: not working atm. Goal: component-wise multiply two lattice and output a new copy.
        lattice_3D<T>& operator*(const lattice_3D<T> &y){  
            lattice_3D<T> ret = *this;
            for(unsigned i = 0; i < LLL; ++i){
                ret.lattice[i] *= y.lattice[i];
            }
            //std::cout << ret << "\n";
            return ret;
        }

        
        T nn_sum(const unsigned i, const unsigned j, const unsigned k){
            return lattice[i*LL + j*L + (k+1)%L] + lattice[i*LL + j*L + (k-1+L)%L] + lattice[i*LL + (j+1)%L *L + k] + lattice[i*LL + (j-1+L)%L *L + k] + lattice[(i+1)%L *LL + j*L + k] + lattice[(i-1+L)%L *LL + j*L + k];  
        }

        T sum(){
            T count = 0;
            for(unsigned i = 0; i < LLL; ++i){
                count += lattice[i];
            }
            return count;
        }

        T sum(const T val){
            T count = 0;
            for(unsigned i = 0; i < LLL; ++i){
                count += lattice[i]*val;
            }
            return count;
        }
        

        unsigned shape(){
            return L;
        }

        
    private:
        unsigned L, LL, LLL;
        T *lattice;
};

template<typename T>
std::ostream& operator<<(std::ostream& os, lattice_3D<T>& x){
    const unsigned L = x.shape();

    for(unsigned i = 0; i < L; ++i){
        for(unsigned j = 0; j < L; ++j){
            for(unsigned k = 0; k < L; ++k){
                const T s = x(i, j, k);
                const std::string str = s < 0 ? " " : "  "; 
                os << str << s;
            }
            os << "\n";
        }
        os << "\n";
    }
    return os;
}

/*
int main(){
    lattice_3D<int> x(3);
    x.fill(1);
    x(1,1,1) *= -1;
    std::cout << "x =\n" << x << "\n";

    lattice_3D<int> y = x;
    y(0,0,1) *= -1;
    std::cout << "y =\n" << y << "\n";
    std::cout << "x*y =\n" << x(13) << "\n";
    std::cout << "x =\n" << x << "\n";
}
*/