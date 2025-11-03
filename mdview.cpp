/* 
The following is an implementation of a multi-array class. The main functionality is to provide
a mapping from linearly stored memory into the space of a multi-dimensional
array. This is intended as the base class for a general tensor class.

The class object 'marray' containing a multi-array of doubles is defined with the statement:

marray<double> object_name(extent1, extent2, ... , extent'n');

Where the extent parameters are integer counts of the length of the array in each distinct 
'dimension'. The above object would be of dimension n. If we were implementing a tensor it
would be a tensor of order n.

There is a multi-indexed access and linearly indexed access via
the data() member function and the data member lumen respectively.
*/

#include <vector>
#include <algorithm>
#include <stdexcept>

template <typename data_type> 
class mdview
{
    //  DATA MEMBERS...  
    private:
        // std::vector demarcating the extents of the multi-array.
        std::vector<std::size_t> extents;
        // std::vector holding the position to be accessed.
        std::vector<std::size_t> position;
        // The number of extents... Should be interpreted as the dimensionality of the view.
        std::size_t order;
        
    public:
        // Reference to linear 'source' array.
        std::vector<data_type>& source;
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

    // Extents extraction from parameter list: "std::size_t... dimensions".
        template<typename... index_type>
        std::vector<std::size_t> extractor(index_type... args){
            return {std::forward<std::size_t>(args)...};
        }

    public:
        template<typename... index_type>
        data_type& operator()(index_type... coordinates){
            // Uses the extractor member function to unfold the pack and apply the coordinates to the positon vector.
            position = extractor(coordinates...);
            // Argument check...
            check_valid_position();
            // Returns the linear index of the coordinates.
            std::size_t offset = linearize();
            // Returns a pointer to the location in memory corresponding to the requested position coordinates.
            return *(source.begin() + offset);
        }
        
    // Constructors...
        template<typename... index_type>
        mdview(std::vector<data_type>& source, index_type... dimensions): source(source), extents(extractor(dimensions...)), order(sizeof...(dimensions)), element_count((dimensions * ...)){
            // Set additional capacity for 'source' for potential accumulation operations.           
            // Sets size of the position vector.
            position.resize(order);
        }
};
