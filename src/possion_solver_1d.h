#ifndef __possion_solver_1d__
#define __possion_solver_1d__

#include<mesh.h>
#include<finite_element.h>
#include<assemble.h>
#include<Eigen/Dense>

using namespace Eigen;

void treat_dirichlet_boundary_condition(MatrixXd& A, VectorXd& b,finite_element fe, double (*dirichlet_function)(double));
VectorXd compute_exact_solution(double (*exact)(double), MatrixXd Pb);


void possion_solver_1d(possion_data_1d pde, mesh msh, int basis_type, int gauss_type){
    int basis_type_test = basis_type;
    int basis_type_trial = basis_type;
    // Todo 生成有限元节点和网格
    finite_element fe(msh, basis_type, "dirichlet");
    // Todo 组装刚度矩阵
    MatrixXd A = assemble_matrix_1d(pde.coefficient_function, msh, fe, basis_type_trial, 1, basis_type_test, 1, gauss_type);
    // cout << A << endl;
    // Todo 组装右端项
    VectorXd b = assemble_vector_1d(pde.right_hand_side_function, msh, fe, basis_type_test, 0, gauss_type);
    // cout << b << endl;
    // Todo 边界条件处理
    treat_dirichlet_boundary_condition(A, b, fe, pde.dirichlet_function);
    cout << A << endl << b << endl;
    // Todo 求解线性方程组
    VectorXd x = A.lu().solve(b);
    cout << x << endl;
    // x = A.ldlt().solve(b);
    // cout << x << endl;
    // Todo 误差计算
    VectorXd exact_solution = compute_exact_solution(pde.exact_solution_function, fe.Pb);
    cout << exact_solution-x << endl;
    double l2_error = compute_Hs_error(pde.exact_solution_function, fe.Pb);
    // Todo 输出结果到文件

}

VectorXd compute_exact_solution(double (*exact)(double), MatrixXd Pb){
    VectorXd exact_solution = VectorXd::Zero(Pb.cols());
    for (int i = 0; i < Pb.cols(); ++ i) {
        exact_solution(i) = exact(Pb(i));
    }
    return exact_solution;
}



void treat_dirichlet_boundary_condition(MatrixXd& A, VectorXd& b,finite_element fe, double (*dirichlet_function)(double)) {
    int number_of_boundary_nodes = fe.boundary_nodes.cols();
    for (int i = 0; i < number_of_boundary_nodes; ++ i) {
        if (fe.boundary_nodes(0, i) == -1) {
            for (int j = 0; j < A.cols(); ++ j) {
                A(fe.boundary_nodes(1, i), j) = 0.0;
            }
            A(fe.boundary_nodes(1, i), fe.boundary_nodes(1, i)) = 1.0;
            b(fe.boundary_nodes(1, i)) = dirichlet_function(fe.Pb(fe.boundary_nodes(1, i)));
        }
    }
}

#endif