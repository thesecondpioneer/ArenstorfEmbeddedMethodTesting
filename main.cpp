#include <iostream>
#include <cmath>
#include <vector>

const double mu = 0.012277471, mus = 1.0 - mu;

//method coefficients
const double a2 = 1.0 / 4, a3 = 3.0 / 8, a4 = 12.0 / 13, a5 = 1.0, a6 = 1.0 / 2,
             b21 = 1.0 / 4, b31 = 3.0 / 32, b41 = 1932.0 / 2197, b51 = 439.0 / 216, b61 = -8.0 / 27,
             b32 = 9.0 / 32, b42 = -7200.0 / 2197, b52 = -8.0, b62 = 2.0,
             b43 = 7296.0 / 2197, b53 = 3680.0 / 513, b63 = -3544.0 / 2565,
             b54 = -845.0 / 4104, b64 = 1859.0 / 4104,
             b65 = -11.0 / 40,
             c1 = 25.0 / 216, c3 = 1408.0 / 2565, c4 = 2197.0 / 4104, c5 =  -1.0 / 5,
             ch1 = 16.0 / 135, ch3 = 6656.0 / 12825, ch4 = 28561.0 / 56430, ch5 = -9.0 / 50, ch6 = 2.0 / 55,
             ct1 = -1.0 / 360, ct3 = 128.0 / 4275, ct4 = 2197.0 / 75240, ct5 = -1.0 / 50, ct6 = -2.0 / 55;

std::vector<double> operator+(std::vector<double> a, std::vector<double> b) {
    for (int i = 0; i < a.size(); i++) {
        a[i] += b[i];
    }
    return a;
}

std::vector<double> operator-(std::vector<double> a, std::vector<double> b) {
    for (int i = 0; i < a.size(); i++) {
        a[i] -= b[i];
    }
    return a;
}

std::vector<double> operator*(double a, std::vector<double> b) {
    for (int i = 0; i < b.size(); i++) {
        b[i] *= a;
    }
    return b;
}

double eunorm(std::vector<double> a) {
    double sum(0);
    for (int i = 0; i < a.size(); i++) {
        sum += a[i] * a[i];
    }
    return std::sqrt(sum);
}

std::vector<double> f(double x, std::vector<double> y) {
    double d1 = std::pow((y[0] + mu) * (y[0] + mu) + y[2] * y[2], 1.5),
            d2 = std::pow((y[0] - mus) * (y[0] - mus) + y[2] * y[2], 1.5);
    return std::vector<double>{y[1],
                               y[0] + 2 * y[3] - mus * (y[0] + mu) / d1 - mu * (y[0] - mu) / d2,
                               y[3],
                               y[2] - 2 * y[1] - mus * y[2] / d1 - mu * y[2] / d2
    };
}

std::vector<double> euler_step(double x0, const std::vector<double> &y0, double h) {
    return y0 + h * f(x0, y0);
}

void rk45_step(double &x0, std::vector<double> &y0, double &h, const double tol) {
    std::vector<double> k1 = h * f(x0, y0),
                        k2 = h * f(x0 + a2 * h, y0 + b21 * k1),
                        k3 = h * f(x0 + a3 * h, y0 + b31 * k1 + b32 * k2),
                        k4 = h * f(x0 + a4 * h, y0 + b41 * k1 + b42 * k2 + b43 * k3),
                        k5 = h * f(x0 + a5 * h, y0 + b51 * k1 + b52 * k2 + b53 * k3 + b54 * k4),
                        k6 = h * f(x0 + a6 * h, y0 + b61 * k1 + b62 * k2 + b63 * k3 + b64 * k4 + b65 * k5);
    double te = eunorm(ct1 * k1 + ct3 * k3 + ct4 * k4 + ct5 * k5 + ct6 * k6);
    h = 0.9 * h * std::pow((tol / te), 1.0 / 5);
    if (te <= tol){
        y0 = y0 + ch1 * k1 + ch3 * k3 + ch4 * k4 + ch5 * k5 + ch6 * k6;
        x0 += h;
    }
}

int main() {
    std::vector<double> y0{0.994, 0, 0, -2.001585106379083}, y_h1, y_prev;
    double t0 = 0, tfin = 17.065216556015796, h = 1, tol, delta;
    std::cout.precision(16);
    std::cin >> tol;
    while (t0 < tfin) {
        h = std::min(tfin - t0, h);
        y_prev = y0;
        rk45_step(t0, y0, h, tol);
        if (y_prev != y0){
            std::cout << std::fixed << y0[0] << ' ' << y0[2] << std::endl;
        }
    }
    return 0;
}
