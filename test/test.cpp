#include<iostream>
#include<possion_data_1d.h>
#include<Eigen/Dense>
#include<../src/mesh.h>
#include<../src/possion_solver_1d.h>
#include<cmath>

using namespace std;
using namespace Eigen;

int main(){
    possion_data_1d pde;
    for (int i = 3; i < 10; ++ i) {
        mesh msh;
        int number_of_nodes = pow(2, i) + 1;
        msh.generate_pt_1d(0, 1, number_of_nodes);
        int basis_type = 101;
        possion_solver_1d(pde, msh, 101, 4);
    }
    return 0;
}