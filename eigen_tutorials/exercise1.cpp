#include<iostream>
#include<Eigen/Dense>

using namespace std;
using namespace Eigen;

int main(){
    MatrixXd m1(2, 2);
    m1(0, 0) = 3;
    m1(1, 0) = 2.5;
    m1(0, 1) = -1;
    m1(1, 1) = m1(1, 0) + m1(0, 1);
    cout << m1 << endl;
    
    MatrixXd m = MatrixXd::Random(3, 3);
    m = (m + MatrixXd::Constant(3, 3, 1.2)) * 50;
    cout << "m = " << endl << m << endl;
    VectorXd v(3);
    v << 1, 2, 3;
    cout << "m * v = " << endl << m * v << endl;
    return 0;
}