#ifndef __error_compute__
#define __error_compute__

double compute_Hs_error(double (*exact)(double), finite_element fe, int der = 0) {
    double error = 0.0;
    int number_of_elements = fe.Tb.cols();
    for (int k = 0; k < number_of_elements; ++ k) {
        
    }
    return error;
}

#endif