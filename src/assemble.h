#ifndef __assemble__
#define __assemble__

#include<Eigen/Dense>
#include<mesh.h>
#include<finite_element.h>
using namespace Eigen;

double gauss_int_trial_test_1d(double (*coef)(double), double left, double right,
                                int basis_type_trial, int basis_index_trial, int der_trial ,
                                int basis_type_test, int basis_index_test, int der_test, int gauss_type);

double gauss_int_test_1d(double (*righthandside)(double), double left, double right,
                            int basis_type_test, int basis_index_test, int der_test, int gauss_type);

MatrixXd assemble_matrix_1d(double (*coef)(double), mesh msh, finite_element fe, int basis_type_trial,
        int der_trial, int basis_type_test, int der_test, int gauss_type){
    int number_of_elements = msh.T.cols();
    int number_of_local_basis = fe.Tb.rows();
    int number_of_finite_element_nodes = fe.Pb.cols();
    MatrixXd A = MatrixXd::Zero(number_of_finite_element_nodes, number_of_finite_element_nodes);
    for (int k = 0; k < number_of_elements; ++ k) {
        int left_global_index = msh.T.col(k)(0);
        int right_global_index = msh.T.col(k)(1);
        double left = msh.P(left_global_index);
        double right = msh.P(right_global_index);
        for (int j = 0; j < number_of_local_basis; ++ j) {
            for (int i = 0; i < number_of_local_basis; ++ i) {
                int row_global_index = fe.Tb.col(k)(i);
                int col_global_index = fe.Tb.col(k)(j);
                A(row_global_index, col_global_index) += gauss_int_trial_test_1d(coef, left, right, basis_type_trial, j, der_trial, basis_type_test, i, der_test, gauss_type);
            }
        }
    }
    return A;
}

VectorXd assemble_vector_1d(double (*righthandside)(double), mesh msh, finite_element fe, int basis_type_test, int der_test, int gauss_type) {
    int number_of_elements = msh.T.cols();
    int number_of_local_basis = fe.Tb.rows();
    int number_of_finite_element_nodes = fe.Pb.cols();
    VectorXd b = VectorXd::Zero(number_of_finite_element_nodes);
    for (int k = 0; k < number_of_elements; ++ k) {
        int left_global_index = msh.T.col(k)(0);
        int right_global_index = msh.T.col(k)(1);
        double left = msh.P(left_global_index);
        double right = msh.P(right_global_index);
        for (int i = 0; i < number_of_local_basis; ++ i) {
            int row_global_index = fe.Tb.col(k)(i);
            b(row_global_index) += gauss_int_test_1d(righthandside, left, right, basis_type_test, i, der_test, gauss_type);
        }
    }
    return b;
}


void generate_gauss_formula_1d(double left, double right, int gauss_type, VectorXd& gauss_weights, VectorXd& gauss_points) {
    switch (gauss_type)
    {
    case 3:
        gauss_weights << 0.5555555556, 0.8888888889, 0.5555555556;
        gauss_points << -0.7745966692, 0, 0.7745966692;
        gauss_weights.array() *= (right - left) / 2.0;
        gauss_points.array() *= (right - left) / 2.0;
        gauss_points.array() += (left + right) / 2.0; 
        break;
    case 4:
        gauss_points << -0.8611363116, -0.3399810436, 0.3399810436, 0.8611363116;
        gauss_weights << 0.3478548461,  0.6521451549, 0.6521451549, 0.3478548461;
        gauss_weights.array() *= (right - left) / 2.0;
        gauss_points.array() *= (right - left) / 2.0;
        gauss_points.array() += (left + right) / 2.0;
    default:
        break;
    }
}

double local_basis_function_1d(double x, double left, double right, int basis_type, int basis_index, int der) {
    double result = 0.0;
    if (basis_type == 101) {
        if (der == 0) {
            if (basis_index == 0) {
                result = (right - x) / (right - left);
            } else if (basis_index == 1) {
                result = (x - left) / (right - left);
            }
        } else if (der == 1) {
            if (basis_index == 0) {
                result = -1.0 / (right - left);
            } else if (basis_index == 1) {
                result = 1.0 / (right - left);
            }
        }
    }
    return result;
}

double gauss_int_test_1d(double (*righthandside)(double), double left, double right,
                            int basis_type_test, int basis_index_test, int der_test, int gauss_type) {
    double result = 0.0;
    VectorXd gauss_weights(gauss_type);
    VectorXd gauss_points(gauss_type);
    generate_gauss_formula_1d(left, right, gauss_type, gauss_weights, gauss_points);
    for (int i = 0; i < gauss_weights.size(); ++ i) {
        result += gauss_weights(i) * righthandside(gauss_points(i)) * local_basis_function_1d(gauss_points(i), left, right, basis_type_test, basis_index_test, der_test);
    }
    return result;
}

double gauss_int_trial_test_1d(double (*coef)(double), double left, double right,
                                int basis_type_trial, int basis_index_trial, int der_trial ,
                                int basis_type_test, int basis_index_test, int der_test, int gauss_type){
    double result = 0.0;
    VectorXd gauss_weights(gauss_type);
    VectorXd gauss_points(gauss_type);
    generate_gauss_formula_1d(left, right, gauss_type, gauss_weights, gauss_points);
    for (int i = 0; i < gauss_weights.size(); ++ i) {
        result += gauss_weights(i) * coef(gauss_points(i)) * local_basis_function_1d(gauss_points(i), left, right, basis_type_trial, basis_index_trial, der_trial)
                                                           * local_basis_function_1d(gauss_points(i), left, right, basis_type_test, basis_index_test, der_test); 
    }
    return result;
}

#endif