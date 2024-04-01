#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <iomanip>

const double mu = 0.012277471, mus = 1.0 - mu;
bool step_made = false;
int calls = 0;

//rk45 method coefficients
const double rk45_a2 = 1.0 / 4, rk45_a3 = 3.0 / 8, rk45_a4 = 12.0 / 13, rk45_a5 = 1.0, rk45_a6 = 1.0 / 2,
        rk45_b21 = 1.0 / 4, rk45_b31 = 3.0 / 32, rk45_b41 = 1932.0 / 2197, rk45_b51 = 439.0 / 216, rk45_b61 = -8.0 / 27,
        rk45_b32 = 9.0 / 32, rk45_b42 = -7200.0 / 2197, rk45_b52 = -8.0, rk45_b62 = 2.0,
        rk45_b43 = 7296.0 / 2197, rk45_b53 = 3680.0 / 513, rk45_b63 = -3544.0 / 2565,
        rk45_b54 = -845.0 / 4104, rk45_b64 = 1859.0 / 4104,
        rk45_b65 = -11.0 / 40,
        rk45_c1 = 25.0 / 216, rk45_c3 = 1408.0 / 2565, rk45_c4 = 2197.0 / 4104, rk45_c5 = -1.0 / 5,
        rk45_ch1 = 16.0 / 135, rk45_ch3 = 6656.0 / 12825, rk45_ch4 = 28561.0 / 56430, rk45_ch5 = -9.0 / 50, rk45_ch6 =
        2.0 / 55,
        rk45_ct1 = -1.0 / 360, rk45_ct3 = 128.0 / 4275, rk45_ct4 = 2197.0 / 75240, rk45_ct5 = -1.0 / 50, rk45_ct6 =
        -2.0 / 55;
//dopri54 method coefficients
const double dopri54_a2 = 1.0 / 5, dopri54_a3 = 3.0 / 10, dopri54_a4 = 4.0 / 5, dopri54_a5 = 8.0 / 9, dopri54_a6 = 1.0,
        dopri54_a7 = 1.0,
        dopri54_b21 = 1.0 / 5, dopri54_b31 = 3.0 / 40, dopri54_b41 = 44.0 / 45, dopri54_b51 =
        19372.0 / 6561, dopri54_b61 = 9017.0 / 3168,
        dopri54_b71 = 35.0 / 384,
        dopri54_b32 = 9.0 / 40, dopri54_b42 = -56.0 / 15, dopri54_b52 = -25360.0 / 2187, dopri54_b62 = -355.0 / 33,
        dopri54_b72 = 0,
        dopri54_b43 = 32.0 / 9, dopri54_b53 = 64448.0 / 6561, dopri54_b63 = 46732.0 / 5247, dopri54_b73 = 500.0 / 1113,
        dopri54_b54 = -212.0 / 729, dopri54_b64 = 49.0 / 176, dopri54_b74 = 125.0 / 192,
        dopri54_b65 = -5103.0 / 18656, dopri54_b75 = -2187.0 / 6784,
        dopri54_b76 = 11.0 / 84,
        dopri54_c1 = 35.0 / 384, dopri54_c3 = 500.0 / 1113, dopri54_c4 = 125.0 / 192, dopri54_c5 =
        -2187.0 / 6784, dopri54_c6 = 11.0 / 84,
        dopri54_ch1 = 5179.0 / 57600, dopri54_ch3 = 7571.0 / 16695, dopri54_ch4 = 393.0 / 640, dopri54_ch5 =
        -92097.0 / 339200,
        dopri54_ch6 = 187.0 / 2100, dopri54_ch7 = 1.0 / 40,
        dopri54_ct1 = 71.0 / 57600, dopri54_ct3 = -71.0 / 16695, dopri54_ct4 = 71.0 / 1920, dopri54_ct5 =
        -17253.0 / 339200,
        dopri54_ct6 = 22.0 / 525, dopri54_ct7 = -1.0 / 40;


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

double infnorm(std::vector<double> a) {
    double result = -1;
    for (int i = 0; i < a.size(); i++) {
        result = std::max(result, std::abs(a[i]));
    }
    return result;
}

std::vector<double> f(double x, std::vector<double> y) {
    double d1 = std::pow((y[0] + mu) * (y[0] + mu) + y[2] * y[2], 1.5),
            d2 = std::pow((y[0] - mus) * (y[0] - mus) + y[2] * y[2], 1.5);
    calls++;
    return std::vector<double>{y[1],
                               y[0] + 2 * y[3] - mus * (y[0] + mu) / d1 - mu * (y[0] - mus) / d2,
                               y[3],
                               y[2] - 2 * y[1] - mus * y[2] / d1 - mu * y[2] / d2
    };
}

std::vector<double>
euler_step(std::vector<double> (*f)(double, std::vector<double>), double x0, const std::vector<double> &y0, double h) {
    return y0 + h * f(x0, y0);
}

void rk45_step(std::vector<double> (*f)(double, std::vector<double>), double &x0, std::vector<double> &y0, double &h,
               const double &tol, const double &x_fin) {
    std::vector<double> k1 = h * f(x0, y0),
            k2 = h * f(x0 + rk45_a2 * h, y0 + rk45_b21 * k1),
            k3 = h * f(x0 + rk45_a3 * h, y0 + rk45_b31 * k1 + rk45_b32 * k2),
            k4 = h * f(x0 + rk45_a4 * h, y0 + rk45_b41 * k1 + rk45_b42 * k2 + rk45_b43 * k3),
            k5 = h * f(x0 + rk45_a5 * h, y0 + rk45_b51 * k1 + rk45_b52 * k2 + rk45_b53 * k3 + rk45_b54 * k4),
            k6 =
            h * f(x0 + rk45_a6 * h, y0 + rk45_b61 * k1 + rk45_b62 * k2 + rk45_b63 * k3 + rk45_b64 * k4 + rk45_b65 * k5);
    double te = infnorm(rk45_ct1 * k1 + rk45_ct3 * k3 + rk45_ct4 * k4 + rk45_ct5 * k5 + rk45_ct6 * k6);
    h = std::min(0.9 * h * std::pow((tol / te), 1.0 / 5), x_fin - x0);
    if (te <= tol) {
        y0 = y0 + rk45_c1 * k1 + rk45_c3 * k3 + rk45_c4 * k4 + rk45_c5 * k5;
        x0 += h;
        step_made = true;
    }
}

void rkf54_step(std::vector<double> (*f)(double, std::vector<double>), double &x0, std::vector<double> &y0, double &h,
                const double &tol, const double &x_fin) {
    std::vector<double> k1 = h * f(x0, y0),
            k2 = h * f(x0 + rk45_a2 * h, y0 + rk45_b21 * k1),
            k3 = h * f(x0 + rk45_a3 * h, y0 + rk45_b31 * k1 + rk45_b32 * k2),
            k4 = h * f(x0 + rk45_a4 * h, y0 + rk45_b41 * k1 + rk45_b42 * k2 + rk45_b43 * k3),
            k5 = h * f(x0 + rk45_a5 * h, y0 + rk45_b51 * k1 + rk45_b52 * k2 + rk45_b53 * k3 + rk45_b54 * k4),
            k6 =
            h * f(x0 + rk45_a6 * h, y0 + rk45_b61 * k1 + rk45_b62 * k2 + rk45_b63 * k3 + rk45_b64 * k4 + rk45_b65 * k5);
    double te = infnorm(rk45_ct1 * k1 + rk45_ct3 * k3 + rk45_ct4 * k4 + rk45_ct5 * k5 + rk45_ct6 * k6);
    h = std::min(0.9 * h * std::pow((tol / te), 1.0 / 5), x_fin - x0);
    if (te <= tol) {
        y0 = y0 + rk45_ch1 * k1 + rk45_ch3 * k3 + rk45_ch4 * k4 + rk45_ch5 * k5 + rk45_ch6 * k6;
        x0 += h;
        step_made = true;
    }
}

void dopri54_step(std::vector<double> (*f)(double, std::vector<double>), double &x0, std::vector<double> &y0, double &h,
                  const double &tol, const double &x_fin, std::vector<double> &f_last) {
    std::vector<double> k1;
    if (std::fpclassify(f_last[0]) == FP_NAN) {
        k1 = h * f(x0, y0);
    } else {
        k1 = h * f_last;
    }
    std::vector<double>
            k2 = h * f(x0 + dopri54_a2 * h, y0 + dopri54_b21 * k1),
            k3 = h * f(x0 + dopri54_a3 * h, y0 + dopri54_b31 * k1 + dopri54_b32 * k2),
            k4 = h * f(x0 + dopri54_a4 * h, y0 + dopri54_b41 * k1 + dopri54_b42 * k2 + dopri54_b43 * k3),
            k5 =
            h * f(x0 + dopri54_a5 * h, y0 + dopri54_b51 * k1 + dopri54_b52 * k2 + dopri54_b53 * k3 + dopri54_b54 * k4),
            k6 = h * f(x0 + dopri54_a6 * h,
                       y0 + dopri54_b61 * k1 + dopri54_b62 * k2 + dopri54_b63 * k3 + dopri54_b64 * k4 +
                       dopri54_b65 * k5);
    f_last = f(x0 + dopri54_a7 * h,
               y0 + dopri54_b71 * k1 + dopri54_b72 * k2 + dopri54_b73 * k3 + dopri54_b74 * k4 + dopri54_b75 * k5 +
               dopri54_b76 * k6);
    std::vector<double> k7 = h * f_last;
    double te = infnorm(dopri54_ct1 * k1 + dopri54_ct3 * k3 + dopri54_ct4 * k4 + dopri54_ct5 * k5 + dopri54_ct6 * k6 +
                        dopri54_ct7 * k7);
    h = std::min(0.9 * h * std::pow((tol / te), 1.0 / 5), x_fin - x0);
    if (te <= tol) {
        y0 = y0 + dopri54_c1 * k1 + dopri54_c3 * k3 + dopri54_c4 * k4 + dopri54_c5 * k5 + dopri54_c6 * k6;
        x0 += h;
        step_made = true;
    }
}

std::vector<std::vector<double>> A11(7), A12(7), A21(7), A22(7), B12(7), B21(7);
std::vector<double> C1, C2, D1, Q1, Di1, Qi1, D2, Q2, Di2, Qi2;

int main() {
    std::vector<double> y0{0.994, 0, 0, -2.001585106379083}, y_h1, yr, f_last(4,
                                                                              std::numeric_limits<double>::quiet_NaN());
    int p = 4;
    double t0 = 0, tfin = 17.065216556015796, h, tol, delta, x0 = 0;
    std::cin >> tol;
    //Выбор начального шага
    delta = std::pow(1.0 / std::max(std::abs(t0), std::abs(tfin)), p + 1) + std::pow(eunorm(f(t0, y0)), p + 1);
    h = std::pow(tol / delta, 1.0 / (p + 1));
    y_h1 = euler_step(f, x0, y0, h);
    delta = std::pow(1.0 / std::max(std::abs(t0), std::abs(tfin)), p + 1) + std::pow(eunorm(f(x0 + h, y_h1)), p + 1);
    h = std::pow(tol / delta, 1.0 / (p + 1));

    std::string mode;
    std::cin >> mode;
    if (mode == "rkf45") {
        std::ofstream fout1("E:\\All\\Work\\Code\\SPbU\\NIR\\Arenstorf\\output.csv"),
                fout2("E:\\All\\Work\\Code\\SPbU\\NIR\\Arenstorf\\output3.csv");

        fout1 << std::setprecision(16) << std::fixed << y0[0] << ',' << y0[2] << ',' << t0 << std::endl;
        fout2 << std::setprecision(16) << std::fixed << h << ',' << t0 << std::endl;
        while (t0 < tfin) {
            step_made = false;
            rk45_step(f, t0, y0, h, tol, tfin);
            fout2 << std::setprecision(16) << std::fixed << h << ',' << t0 << std::endl;
            if (step_made) {
                fout1 << std::setprecision(16) << std::fixed << y0[0] << ',' << y0[2] << ',' << t0 << std::endl;
            }
        }
    } else if (mode == "dopri") {
        std::ofstream fout1("E:\\All\\Work\\Code\\SPbU\\NIR\\Arenstorf\\output1.csv"),
                fout2("E:\\All\\Work\\Code\\SPbU\\NIR\\Arenstorf\\output4.csv");
        fout1 << std::setprecision(16) << std::fixed << y0[0] << ',' << y0[2] << ',' << t0 << std::endl;
        fout2 << std::setprecision(16) << std::fixed << h << ',' << t0 << std::endl;
        while (t0 < tfin) {
            step_made = false;
            dopri54_step(f, t0, y0, h, tol, tfin, f_last);
            fout2 << std::setprecision(16) << std::fixed << h << ',' << t0 << std::endl;
            if (step_made) {
                fout1 << std::setprecision(16) << std::fixed << y0[0] << ',' << y0[2] << ',' << t0 << std::endl;
            }
        }
    } else if (mode == "rkf54") {
        std::ofstream fout1("E:\\All\\Work\\Code\\SPbU\\NIR\\Arenstorf\\output2.csv"),
                fout2("E:\\All\\Work\\Code\\SPbU\\NIR\\Arenstorf\\output5.csv");
        fout1 << std::setprecision(16) << std::fixed << y0[0] << ',' << y0[2] << ',' << t0 << std::endl;
        fout2 << std::setprecision(16) << std::fixed << h << ',' << t0 << std::endl;
        while (t0 < tfin) {
            step_made = false;
            rkf54_step(f, t0, y0, h, tol, tfin);
            fout2 << std::setprecision(16) << std::fixed << h << ',' << t0 << std::endl;
            if (step_made) {
                fout1 << std::setprecision(16) << std::fixed << y0[0] << ',' << y0[2] << ',' << t0 << std::endl;
            }
        }
    } else if (mode == "struct") {
        C1 = {0, 1.0 / 5, 1.0 / 5, 3.0 / 10, 8.0 / 11, 1, 1};
        C2 = {0, 1.0 / 5, 1.0 / 3, 1.0 / 4, 8.0 / 11, 1, 1};
        D1 = {23.0 / 288, 0, 25.0 / 348, 100.0 / 423, 14641.0 / 130848, 0, 0};
        Q1 = {23.0 / 288, 0, 125.0 / 1392, 1000.0 / 2961, 161051.0 / 392544, 83.0 / 1008, 0};
        Di1 = {343.0 / 3456, 0, -185.0 / 5568, 5995.0 / 17766, 422411.0 / 4710528, 83.0 / 12096, 0};
        Qi1 = {1.0 / 12, 0, 25.0 / 348, 50.0 / 141, 6655.0 / 16356, 0, 1.0 / 12};
        D2 = {13.0 / 160, 0, 27.0 / 260, 64.0 / 315, 14641.0 / 131040, 0, 0};
        Q2 = {13.0 / 160, 0, 81.0 / 520, 256.0 / 945, 161051.0 / 393120, 89.0 / 1080, 0};
        Di2 = {-4817.0 / 21120, 975.0 / 176, 65947.0 / 22880, -241216.0 / 31185, 161051.0 / 4717440, 89.0 / 12960, 0};
        Qi2 = {1.0 / 12, 0, 9.0 / 52, 16.0 / 63, 1331.0 / 3276, 0, 1.0 / 12};
        A11[0] = {0};
        A11[1] = {0};
        A11[2] = {1.0 / 100, 1.0 / 100};
        A11[3] = {19.0 / 800, 61.0 / 2400, -1.0 / 240};
        A11[4] = {14350.0 / 131769, -7922.0 / 131769, -15196.0 / 131769, 43616.0 / 131769};
        A11[5] = {-124.0 / 747, 551.0 / 2988, 20378.0 / 21663, -22624.0 / 35109, 83853.0 / 452516};
        A11[6] = {23.0 / 288, 0, 25.0 / 348, 100.0 / 423, 14641.0 / 130848, 0};
        A12[0] = {0};
        A12[1] = {1.0 / 50};
        A12[2] = {1.0 / 100, 1.0 / 100};
        A12[3] = {191.0 / 8000, 67.0 / 3200, 3.0 / 16000};
        A12[4] = {-34016.0 / 483153, 956885.0 / 483153, 768413.0 / 805255, -2093568.0 / 805255};
        A12[5] = {-1819.0 / 4980, 8975.0 / 996, 3991.0 / 830, -5376.0 / 415, 0};
        A12[6] = {13.0 / 160, 0, 27.0 / 260, 64.0 / 315, 14641.0 / 131040, 0};
        A21[0] = {0};
        A21[1] = {1.0 / 100, 1.0 / 100};
        A21[2] = {1.0 / 36, 1.0 / 36, 0};
        A21[3] = {239.0 / 13824, 59.0 / 3072, -55.0 / 9216, 5.0 / 6912};
        A21[4] = {60979.0 / 1581228, -2668.0 / 43923, 61100.0 / 347391, 1711840.0 / 18579429, 12285.0 / 659692};
        A21[5] = {7055.0 / 38448, 199.0 / 1068, -11755.0 / 23229, 61400.0 / 112941, 179685.0 / 1940912};
        A21[6] = {23.0 / 288, 0, 25.0 / 348, 100.0 / 423, 14641.0 / 130848, 0};
        A22[0] = {0};
        A22[1] = {1.0 / 50};
        A22[2] = {1.0 / 54, 1.0 / 27};
        A22[3] = {161.0 / 7680, 31.0 / 3072, 1.0 / 5120};
        A22[4] = {-166432.0 / 2415765, 957302.0 / 483153, 769574.0 / 805255, -419328.0 / 161051};
        A22[5] = {-983.0 / 2670, 2405.0 / 267, 2138.0 / 445, -1152.0 / 89, 0, 0};
        A22[6] = {13.0 / 160, 0, 27.0 / 260, 64.0 / 315, 14641.0 / 131040, 0};
        B12[0] = {0};
        B12[1] = {1.0 / 5};
        B12[2] = {1.0 / 10, 1.0 / 10};
        B12[3] = {27.0 / 400, 39.0 / 160, -9.0 / 800};
        B12[4] = {-47852.0 / 73205, 195970.0 / 14641, 516954.0 / 73205, -1395712.0 / 7305};
        B12[5] = {1601.0 / 415, -11255.0 / 166, -36555.0 / 1079, 40448.0 / 415, 14641.0 / 10790};
        B12[6] = {13.0 / 160, 0, 81.0 / 520, 256.0 / 945, 161051.0 / 393120, 89.0 / 1080};
        B21[0] = {0};
        B21[1] = {1.0 / 10, 1.0 / 10};
        B21[2] = {1.0 / 18, -5.0 / 54, 10.0 / 27};
        B21[3] = {49.0 / 576, 5.0 / 128, 55.0 / 384, -5.0 / 288};
        B21[4] = {329.0 / 2178, 20.0 / 1331, -40240.0 / 115797, 39520.0 / 51183, 4095.0 / 29986};
        B21[5] = {-1067.0 / 6408, -5.0 / 178, 11060.0 / 7743, -34360.0 / 37647, 658845.0 / 970456};
        B21[6] = {23.0 / 288, 0, 125.0 / 1392, 1000.0 / 2961, 161051.0 / 392544, 83.0 / 1008};
    }
    return 0;
}