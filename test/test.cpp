#include<iostream>
#include<possion_data_1d.h>
#include<Eigen/Dense>
#include<../src/mesh.h>
#include<../src/possion_solver_1d.h>

using namespace std;
using namespace Eigen;

int main(){
    possion_data_1d pde;
    mesh msh;
    msh.generate_pt_1d(0, 1, 10);
    cout << msh.P << endl;
    cout << msh.T << endl;
    int basis_type = 101;
    possion_solver_1d(pde, msh, 101, 4);
    return 0;
}