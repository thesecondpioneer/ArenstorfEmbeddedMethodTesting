#include <iostream>
#include "RK45Solver.h"
#include "RK54Solver.h"
#include "DOPRI54Solver.h"
#include "SRK64Solver.h"

const double mu = 0.012277471, mus = 1.0 - mu;

Eigen::VectorXd f(double x, Eigen::VectorXd y) {
    Eigen::VectorXd result(4);
    double d1 = std::pow((y[0] + mu) * (y[0] + mu) + y[2] * y[2], 1.5),
            d2 = std::pow((y[0] - mus) * (y[0] - mus) + y[2] * y[2], 1.5);
    result(0) = y[1];
    result(1) = y[0] + 2 * y[3] - mus * (y[0] + mu) / d1 - mu * (y[0] - mus) / d2;
    result(2) = y[3];
    result(3) = y[2] - 2 * y[1] - mus * y[2] / d1 - mu * y[2] / d2;
    return result;
}

double F(double x, Eigen::VectorXd y) {
    double d1 = std::pow((y[0] + mu) * (y[0] + mu) + y[2] * y[2], 1.5),
            d2 = std::pow((y[0] - mus) * (y[0] - mus) + y[2] * y[2], 1.5);
    return y[0] + 2 * y[3] - mus * (y[0] + mu) / d1 - mu * (y[0] - mus) / d2;
}

double G(double x, Eigen::VectorXd y) {
    double d1 = std::pow((y[0] + mu) * (y[0] + mu) + y[2] * y[2], 1.5),
            d2 = std::pow((y[0] - mus) * (y[0] - mus) + y[2] * y[2], 1.5);
    return y[2] - 2 * y[1] - mus * y[2] / d1 - mu * y[2] / d2;
}

int main() {
    double t0 = 0, tfin = 17.065216556015796;
    Eigen::VectorXd y(4);
    y(0) = 0.994;
    y(1) = 0;
    y(2) = 0;
    y(3) = -2.001585106379083;
    RK45Solver<double> solver1(t0, tfin, 1e-6, y,
                               std::function<Eigen::VectorX<double>(double, Eigen::VectorX<double>)>(f),
                               true, "../Outputs/output.csv", "../Outputs/output1.csv", "../Outputs/Trust/rkf45.csv");
    RK54Solver<double> solver2(t0, tfin, 1e-6, y,
                               std::function<Eigen::VectorX<double>(double, Eigen::VectorX<double>)>(f),
                               true, "../Outputs/output4.csv", "../Outputs/output5.csv", "../Outputs/Trust/rkf54.csv");
    DOPRI54Solver<double> solver3(t0, tfin, 1e-6, y,
                                  std::function<Eigen::VectorX<double>(double, Eigen::VectorX<double>)>(f),
                                  true, "../Outputs/output2.csv", "../Outputs/output3.csv",
                                  "../Outputs/Trust/dopri54.csv");
    SRK64Solver<double> solver4(t0, tfin, 1e-6, y, std::function<double(double, Eigen::VectorX<double>)>(F),
                                std::function<double(double, Eigen::VectorX<double>)>(G),
                                true, "../Outputs/output6.csv", "../Outputs/output7.csv",
                                "../Outputs/Trust/srk64.csv");
    solver1.solve();
    solver2.solve();
    solver3.solve();
    solver4.solve();
    return 0;
}
