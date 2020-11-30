#ifndef __finite_element__
#define __finite_element__

#include<mesh.h>
#include<string>
#include<Eigen/Dense>

using namespace Eigen;
using namespace std;

class finite_element
{
private:
    /* data */
public:
    MatrixXd Pb;
    MatrixXi Tb;
    MatrixXi boundary_nodes;
    finite_element(mesh msh, int basis_type, string flag);
    ~finite_element();
};

finite_element::finite_element(mesh msh, int basis_type, string flag)
{
    switch (basis_type)
    {
    case 101:
        Pb = msh.P;
        Tb = msh.T;
        if (flag == "dirichlet") {
            boundary_nodes = MatrixXi::Zero(2, 2);
            boundary_nodes << -1, -1, 0, Pb.cols() - 1;
        }        
        break;
    
    default:
        break;
    }
}

finite_element::~finite_element()
{
}


#endif