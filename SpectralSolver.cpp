#include <iostream>
#include <numbers>
#include <cmath>
#include <algorithm>
#include <vector>
#include <chrono>
#include <iterator>
#include <stdexcept>

#define N (2 << 4)
#define M (N << 1)

#define A 0.166666666666667
#define B 0.333333333333333


template <typename data_type> 
class marray
{
    //  DATA MEMBERS...  
    private:
        // std::vector demarcating the extents of the multi-array.
        std::vector<std::size_t> extents;
        // std::vector holding the position to be accessed.
        std::vector<std::size_t> position;
        // The number of extents... Can be interpreted as the order of a tensor.
        std::size_t order;
        
    public:
        // Linearly stored data of the multi-array.
        std::vector<data_type> lumen;
        // The number of elements. Must remain constant across all reshapings and transpositions.
        const std::size_t element_count;

    // MEMBER FUNCTIONS...
    private:

    // Checks for invalid access arguments.
        void check_valid_position(){
            if (extents.size() != position.size())
                throw std::invalid_argument("Improper index argument!");
            else
                for (std::size_t i = 0; i != order; i++){
                    if(extents[i] <= position[i]){
                        throw std::invalid_argument("Index argument exceeds the extents of the object!");
                        break;
                    }
                }
            }


    // Returns the linearized memory index of the referenced element: (a,b,c,...).
        std::size_t linearize(){

            // Initialize accumulator values. 
            std::size_t pos_product = 1;
            std::size_t linear_index = 0;

            // Calculates the linear index given the data stored in the 'position' vector.
            for(std::size_t k = 0; k != order; k++){
                for(std::size_t j = k + 1; j != order; j++){
                    pos_product *= extents[j];
                }
                linear_index += position[k] * pos_product;
                pos_product = 1;
            }
            return linear_index;
        }

    //Extents extraction from parameter list: "std::size_t... dimensions".
        template<typename... index_type>
        std::vector<std::size_t> extractor(const index_type&... args){
            return {std::forward<std::size_t>(args)...};
        }

    public:
        template<typename... index_type>
        data_type& operator()(const index_type&... coordinates){
            
            // Uses the extractor member function to unfold the pack and apply the coordinates to the positon vector.
            position = extractor(coordinates...);
            
            // Argument check...
            check_valid_position();

            // Returns the linear index of the coordinates.
            std::size_t offset = linearize();
            
            // Returns a pointer to the location in memory corresponding to the requested position coordinates.
            return *(lumen.begin() + offset);
        }

    // Constructors...
        template<typename... index_type>
        marray(const index_type&... dimensions): extents(extractor(dimensions...)), order(sizeof...(dimensions)), element_count((dimensions * ...)){
            // Set additional capacity for 'lumen' for potential accumulation operations.           
            std::size_t padding = *std::max_element(extents.begin(),extents.end());
            lumen.resize(element_count);
            lumen.reserve(element_count + padding);
            
            // Sets size of the position vector.
            position.resize(order);
        }
};

void discreteFourier(marray<double>& input, bool invert = 0)
{
    marray<double> accumulator(M);

    if (invert == true){

        // Allocate memory for the Inverse Fourier Operator.
        marray<double> inverseOperator(N, M);

        // Initialize Inverse Fourier Operator.
        for (unsigned p = 0; p < N ; ++p){
            for (unsigned q = 0; q < M ; q += 2){
            inverseOperator(p,q) = cos((std::numbers::pi*p*q)/N);
            inverseOperator(p,q+1) = -sin((std::numbers::pi*p*q)/N);
            }
        }

        // Carry out the Inverse Fourier Transform.
        for (unsigned p = 0; p < N ; ++p){
            for (unsigned q = 0; q < M ; q += 2){
            accumulator(2*p) += ((inverseOperator(p,q) * input(q)) - (inverseOperator(p,q+1) * input(q+1))) / N;
            accumulator(2*p+1) += ((inverseOperator(p,q) * input(q+1)) + (inverseOperator(p,q+1) * input(q))) / N;
            }
        }

        input.lumen = accumulator.lumen;

    }else{

        // Allocate memory for the Forward Fourier Operator.
        marray<double> forwardOperator (N, M);

        // Initialize Forward Fourier Operator
        for (unsigned p = 0; p < N ; ++p){
            for (unsigned q = 0; q < M ; q += 2){
            forwardOperator(p,q) = cos((std::numbers::pi*p*q)/N);
            forwardOperator(p,q+1) = sin((std::numbers::pi*p*q)/N);
            }
        }

        // Carry out the Inverse Fourier Transform.
        for (unsigned p = 0; p < N ; ++p){
            for (unsigned q = 0; q < M ; q += 2){
            accumulator(2*p) += ((forwardOperator(p,q) * input(q)) - (forwardOperator(p,q+1) * input(q+1)));
            accumulator(2*p+1) += ((forwardOperator(p,q) * input(q+1)) + (forwardOperator(p,q+1) * input(q)));
            }
        }

        input.lumen = accumulator.lumen;

    }
}

//  Performs a single forward timestep via RK4 method.

double RK4(double(*dydx)(const double, const double), double inity, double stepsize = 0.01)
{

    double next = inity;
    static double t = 0;

    double k1;
    double k2;
    double k3;
    double k4;

    k1 = stepsize * dydx(t, next);
    k2 = stepsize * dydx(t + 0.5 * stepsize, next + 0.5 * k1);
    k3 = stepsize * dydx(t + 0.5 * stepsize, next + 0.5 * k2);
    k4 = stepsize * dydx(t + stepsize, next + k3);
    next = next + (A * k1) + (B * k2) + (B * k3) + (A * k4);
    
    t += stepsize;

    return next;

}

// dydx = f(x,y) populator function.

double fofxy(const double x, const double y)
{
    double fofxy = exp(x);
    return fofxy;
}

int main()
{

}
