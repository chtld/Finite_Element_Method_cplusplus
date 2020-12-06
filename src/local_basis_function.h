#ifndef __local_basis_function__
#define __local_basis_function__

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

#endif