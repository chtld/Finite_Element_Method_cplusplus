#ifndef __error_compute__
#define __error_compute__
#include<Eigen/Dense>
#include<gauss_formula.h>
#include<local_basis_function.h>
#include<cmath>

VectorXd get_local_solution(const VectorXd& solution, VectorXi& index);
double gauss_int_error_1d(VectorXd uhlocal, double (*exact)(double), double left, double right, int basis_type, int der, int gauss_type);
double FE_function_1d(VectorXd uh_local, double x, double left, double right, int basis_type, int der);
double compute_Hs_error(double (*exact)(double), VectorXd& solution, mesh msh, finite_element fe, int basis_type, int der, int gauss_type);

double compute_Hs_error(double (*exact)(double), VectorXd& solution, mesh msh, finite_element fe, int basis_type, int der, int gauss_type) {
    double error = 0.0;
    int number_of_elements = fe.Tb.cols();
    for (int k = 0; k < number_of_elements; ++ k) {
        int left_global_index = msh.T.col(k)(0);
        int right_global_index = msh.T.col(k)(1);
        double left = msh.P(left_global_index);
        double right = msh.P(right_global_index);
        VectorXi index = fe.Tb.col(k);
        // cout << "here" << endl;
        VectorXd uh_local = get_local_solution(solution, index);
        // cout << uh_local << endl;
        error += gauss_int_error_1d(uh_local, exact, left, right, basis_type, der, gauss_type);
    }
    return sqrt(error);
}


VectorXd get_local_solution(const VectorXd& solution, VectorXi& index) {
    VectorXd uh_local(index.size());
    for (int i = 0; i < index.size(); ++ i) {
        uh_local(i) = solution(index(i));
    }
    return uh_local;
}

double gauss_int_error_1d(VectorXd uh_local, double (*exact)(double), double left, double right, int basis_type, int der, int gauss_type){
    double result = 0.0;
    VectorXd gauss_weights(gauss_type);
    VectorXd gauss_points(gauss_type);
    generate_gauss_formula_1d(left, right, gauss_type, gauss_weights, gauss_points);
    for (int i = 0; i < gauss_type; ++ i) {
        result += gauss_weights(i) * (FE_function_1d(uh_local, gauss_points(i), left, right, basis_type, der) - exact(gauss_points(i)))
                                   * (FE_function_1d(uh_local, gauss_points(i), left, right, basis_type, der) - exact(gauss_points(i)));
    }
    return result;
}

double FE_function_1d(VectorXd uh_local, double x, double left, double right, int basis_type, int der) {
    double result = 0.0;
    int number_of_basis_function = uh_local.size();
    for (int i = 0; i < number_of_basis_function; ++ i) {
        result += uh_local(i) * local_basis_function_1d(x, left, right, basis_type, i, der);
    }
    return result;
}

#endif