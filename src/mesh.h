#ifndef __mesh__
#define __mesh__

#include<Eigen/Dense>
using namespace Eigen;

class mesh
{
public:
    void generate_pt_1d(double left, double right, int number_of_nodes);
    mesh();
    ~mesh();
    MatrixXd P;
    MatrixXi T;
};

mesh::mesh()
{
}

mesh::~mesh()
{
}

void mesh::generate_pt_1d(double left, double right, int number_of_nodes){
    P = VectorXd::LinSpaced(number_of_nodes, left, right).transpose();    
    T = MatrixXi::Zero(2, number_of_nodes - 1);
    for (int i = 0; i < T.rows(); ++ i) {
        for (int j = 0; j < T.cols(); ++ j) {
            T(i, j) = j + i;
        }
    }
    return;
}



#endif