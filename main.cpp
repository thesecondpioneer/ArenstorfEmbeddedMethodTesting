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
        c1 = 25.0 / 216, c3 = 1408.0 / 2565, c4 = 2197.0 / 4104, c5 = -1.0 / 5,
        ch1 = 16.0 / 135, ch3 = 6656.0 / 12825, ch4 = 28561.0 / 56430, ch5 = -9.0 / 50, ch6 = 2.0 / 55,
        ct1 = -1.0 / 360, ct3 = 128.0 / 4275, ct4 = 2197.0 / 75240, ct5 = -1.0 / 50, ct6 = -2.0 / 55;

double c2 = 0.45, A = -2, B = 2, C = -2;

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

double infnorm(std::vector<double> a){
    double result = -1;
    for (int i = 0; i < a.size(); i++){
        result = std::max(result, std::abs(a[i]));
    }
    return result;
}

std::vector<double> f(double x, std::vector<double> y) {
    double d1 = std::pow((y[0] + mu) * (y[0] + mu) + y[2] * y[2], 1.5),
            d2 = std::pow((y[0] - mus) * (y[0] - mus) + y[2] * y[2], 1.5);
    return std::vector<double>{y[1],
                               y[0] + 2 * y[3] - mus * (y[0] + mu) / d1 - mu * (y[0] - mus) / d2,
                               y[3],
                               y[2] - 2 * y[1] - mus * y[2] / d1 - mu * y[2] / d2
    };
}

std::vector<double> f1(double x, std::vector<double> y) {
    if ((y[1] < 0) or (y[0] < 0)) {
        std::cout << "Illegal argument, can't use this stepsize" << std::endl;
        return std::vector<double>(4, std::numeric_limits<double>::quiet_NaN());
    }
    return std::vector<double>{2.0 * x * std::pow(y[1], 1.0 / B) * y[3],
                               2.0 * B * x * std::exp(B / C * (y[2] - A)) * y[3],
                               2.0 * C * x * y[3],
                               -2.0 * x * std::log(y[0])
    };
}

std::vector<double> y(double x) {
    double x2 = std::pow(x, 2), sx2 = std::sin(x2), bsx2 = B * sx2, csx2 = C * sx2;
    return std::vector<double>{std::exp(sx2),
                               std::exp(bsx2),
                               csx2 + A,
                               std::cos(x2),
    };
}

std::vector<double> euler_step(std::vector<double> (*f)(double, std::vector<double>), double x0, const std::vector<double> &y0, double h) {
    return y0 + h * f(x0, y0);
}

void rk45_step(std::vector<double> (*f)(double, std::vector<double>), double &x0, std::vector<double> &y0, double &h,
               const double tol) {
    std::vector<double> k1 = h * f(x0, y0),
            k2 = h * f(x0 + a2 * h, y0 + b21 * k1),
            k3 = h * f(x0 + a3 * h, y0 + b31 * k1 + b32 * k2),
            k4 = h * f(x0 + a4 * h, y0 + b41 * k1 + b42 * k2 + b43 * k3),
            k5 = h * f(x0 + a5 * h, y0 + b51 * k1 + b52 * k2 + b53 * k3 + b54 * k4),
            k6 = h * f(x0 + a6 * h, y0 + b61 * k1 + b62 * k2 + b63 * k3 + b64 * k4 + b65 * k5);
    double te = infnorm(ct1 * k1 + ct3 * k3 + ct4 * k4 + ct5 * k5 + ct6 * k6);
    h = 0.9 * h * std::pow((tol / te), 1.0 / 5);
    if (te <= tol) {
        y0 = y0 + ch1 * k1 + ch3 * k3 + ch4 * k4 + ch5 * k5 + ch6 * k6;
        x0 += h;
    }
}

int main() {
    std::vector<double> y0{0.994, 0, 0, -2.001585106379083}, y0t = y(0), y_h1, y_prev, yr;
    double t0 = 0, tfin = 17.065216556015796, h = 0.00001, tol, delta, x0 = 0;
    y0t = y(x0);
    std::cout.precision(16);
    std::cin >> tol;
    while (t0 < tfin) {
        h = std::min(tfin - t0, h);
        y_prev = y0;
        rk45_step(f, t0,y0,h,tol);
        if (y_prev != y0t) {
            std::cout << std::fixed << y0[0] << ' ' << y0[2] << std::endl;
        }
    }
    /*while (x0 < 5.0) {
        h = std::min(5.0 - x0, h);
        y_prev = y0t;
        rk45_step(f1, x0,y0t,h,tol);
        if (y_prev != y0t) {
            std::cout << std::fixed << x0 << ' ' << y0t[3] << std::endl;
        }
    }*/
    /*while (x0 < 5.0) {
        h = std::min(5.0 - x0, h);
        y_prev = y0t;
        y0t = euler_step(f1, x0, y0t, h);
        x0 += h;
        if (y_prev != y0t) {
            std::cout << std::fixed << x0 << ' ' << y0t[0] << ' ' << eunorm(y(x0) - y0t) << std::endl;
        }
    }*/
    return 0;
}
