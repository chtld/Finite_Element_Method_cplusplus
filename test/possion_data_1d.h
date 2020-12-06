#ifndef __possion_data_1d__
#define __possion_data_1d__

#include<cmath>

class possion_data_1d;


class possion_data_1d{
private:
    /* data */
public:
    possion_data_1d();
    ~possion_data_1d();
    static double coefficient_function(double x);
    static double exact_solution_function(double x);
    static double exact_solution_function_grad(double x);
    static double right_hand_side_function(double x);
    static double dirichlet_function(double x);
    static double neumann_function(double x);
};

possion_data_1d::possion_data_1d(){
}

possion_data_1d::~possion_data_1d(){
}

inline double 
possion_data_1d::coefficient_function(double x){
    double result = exp(x);
    return result;
}

inline double 
possion_data_1d::exact_solution_function(double x){
    double result = x * cos(x);
    return result;
}

inline double 
possion_data_1d::exact_solution_function_grad(double x){
    double result = cos(x) - x * sin(x);
    return result;
}

inline double 
possion_data_1d::right_hand_side_function(double x){
    double result = -exp(x) * (cos(x) - 2 * sin(x) - x * cos(x) - x * sin(x));
    return result;
}

inline double
possion_data_1d::dirichlet_function(double x){
    double result = x * cos(x);
    return result;
}

inline double
possion_data_1d::neumann_function(double){
    return 0.0;
}

#endif //__possion_data_1d__