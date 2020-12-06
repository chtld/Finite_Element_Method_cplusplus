#ifndef __gauss_formula__
#define __gauss_formula__

void generate_gauss_formula_1d(double left, double right, int gauss_type, VectorXd& gauss_weights, VectorXd& gauss_points);

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

#endif