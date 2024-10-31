#ifndef EMBEDDED_METHODS_CPP_SRK64SOLVER_H
#define EMBEDDED_METHODS_CPP_SRK64SOLVER_H

#include "EmbeddedRungeKutta.h"
#include <vector>

//This method is purposed for solving structural ODE systems with cross structure relative to the first derivatives
//Which means, a system like
//x'' = F(t, x, y, y')
//y'' = G(t, x, y, x')

//The input vector y for this method should be structured like this:
//y[0]: x
//y[1]: x'
//y[2]: y
//y[3]: y'

template<typename Scalar>
class SRK64Solver : public EmbeddedRungeKutta<Scalar> {
protected:
    //method coefficients
    const std::vector<Scalar> C1 = {0, 1.0 / 5, 1.0 / 5, 3.0 / 10, 8.0 / 11, 1, 1},
            C2 = {0, 1.0 / 5, 1.0 / 3, 1.0 / 4, 8.0 / 11, 1, 1},
            D1 = {23.0 / 288, 0, 25.0 / 348, 100.0 / 423, 14641.0 / 130848, 0, 0},
            Q1 = {23.0 / 288, 0, 125.0 / 1392, 1000.0 / 2961, 161051.0 / 392544, 83.0 / 1008, 0},
            Di1 = {343.0 / 3456, 0, -185.0 / 5568, 5995.0 / 17766, 422411.0 / 4710528, 83.0 / 12096, 0},
            Qi1 = {1.0 / 12, 0, 25.0 / 348, 50.0 / 141, 6655.0 / 16356, 0, 1.0 / 12},
            D2 = {13.0 / 160, 0, 27.0 / 260, 64.0 / 315, 14641.0 / 131040, 0, 0},
            Q2 = {13.0 / 160, 0, 81.0 / 520, 256.0 / 945, 161051.0 / 393120, 89.0 / 1080, 0},
            Di2 = {-4817.0 / 21120, 975.0 / 176, 65947.0 / 22880, -241216.0 / 31185, 161051.0 / 4717440, 89.0 / 12960, 0},
            Qi2 = {1.0 / 12, 0, 9.0 / 52, 16.0 / 63, 1331.0 / 3276, 0, 1.0 / 12};
    const std::vector<std::vector<Scalar>> A11 = {
            {0},
            {0},
            {1.0 / 100, 1.0 / 100},
            {19.0 / 800, 61.0 / 2400, -1.0 / 240},
            {14350.0 / 131769, -7922.0 / 131769, -15196.0 / 131769, 43616.0 / 131769},
            {-124.0 / 747, 551.0 / 2988, 20378.0 / 21663, -22624.0 / 35109, 83853.0 / 452516},
            {23.0 / 288, 0, 25.0 / 348, 100.0 / 423, 14641.0 / 130848, 0}},
            A12 = {
            {0},
            {1.0 / 50},
            {1.0 / 100, 1.0 / 100},
            {191.0 / 8000, 67.0 / 3200, 3.0 / 16000},
            {-34016.0 / 483153, 956885.0 / 483153, 768413.0 / 805255, -2093568.0 / 805255},
            {-1819.0 / 4980, 8975.0 / 996, 3991.0 / 830, -5376.0 / 415, 0},
            {13.0 / 160, 0, 27.0 / 260, 64.0 / 315, 14641.0 / 131040, 0}},
            A21 = {
            {0},
            {1.0 / 100, 1.0 / 100},
            {1.0 / 36, 1.0 / 36, 0},
            {239.0 / 13824, 59.0 / 3072, -55.0 / 9216, 5.0 / 6912},
            {60979.0 / 1581228, -2668.0 / 43923, 61100.0 / 347391, 1711840.0 / 18579429, 12285.0 / 659692},
            {7055.0 / 38448, 199.0 / 1068, -11755.0 / 23229, 61400.0 / 112941, 179685.0 / 1940912, 0},
            {23.0 / 288, 0, 25.0 / 348, 100.0 / 423, 14641.0 / 130848, 0, 0}},
            A22 = {
            {0},
            {1.0 / 50},
            {1.0 / 54, 1.0 / 27},
            {161.0 / 7680, 31.0 / 3072, 1.0 / 5120},
            {-166432.0 / 2415765, 957302.0 / 483153, 769574.0 / 805255, -419328.0 / 161051},
            {-983.0 / 2670, 2405.0 / 267, 2138.0 / 445, -1152.0 / 89, 0},
            {13.0 / 160, 0, 27.0 / 260, 64.0 / 315, 14641.0 / 131040, 0}},
            B12 = {
            {0},
            {1.0 / 5},
            {1.0 / 10, 1.0 / 10},
            {27.0 / 400, 39.0 / 160, -9.0 / 800},
            {-47852.0 / 73205, 195970.0 / 14641, 516954.0 / 73205, -1395712.0 / 73205},
            {1601.0 / 415, -11255.0 / 166, -36555.0 / 1079, 40448.0 / 415, 14641.0 / 10790},
            {13.0 / 160, 0, 81.0 / 520, 256.0 / 945, 161051.0 / 393120, 89.0 / 1080}},
            B21 = {
            {0},
            {1.0 / 10, 1.0 / 10},
            {1.0 / 18, -5.0 / 54, 10.0 / 27},
            {49.0 / 576, 5.0 / 128, 55.0 / 384, -5.0 / 288},
            {329.0 / 2178, 20.0 / 1331, -40240.0 / 115797, 39520.0 / 51183, 4095.0 / 29986},
            {-1067.0 / 6408, -5.0 / 178, 11060.0 / 7743, -34360.0 / 37647, 658845.0 / 970456, 0},
            {23.0 / 288, 0, 125.0 / 1392, 1000.0 / 2961, 161051.0 / 392544, 83.0 / 1008, 0}};

    //Structural groups of the ODE rhs
    std::function<Scalar(Scalar, Eigen::VectorX<Scalar>)> _rhsF;
    std::function<Scalar(Scalar, Eigen::VectorX<Scalar>)> _rhsG;

    bool _firstIter = true;
    Scalar _K1Last, _K2Last;

    void makeStep() override {
        //Calculate the correction functions for the step
        std::vector<Scalar> K1(7), K2(7);
        K1[0] = _firstIter ? this->_rhsF(this->_x, this->_y) : _K1Last;
        K2[0] = _firstIter ? this->_rhsG(this->_x, this->_y) : _K2Last;
        _firstIter = false;
        for (int i = 1; i < 7; i++) {
            double s2x_k1 = 0, s2y_k1 = 0, s2x_k2 = 0, s2y_k2 = 0, s2xs = 0, s2ys = 0;

            for (int j = 0; j < i; j++) {
                s2x_k1 += A11[i][j] * K1[j];
                s2y_k1 += A12[i][j] * K2[j];
                s2ys += B12[i][j] * K2[j];
            }
            Eigen::VectorX<Scalar> yK1(this->_y);
            yK1[0] += C1[i] * this->_h * this->_y[1] + this->_h * this->_h * s2x_k1;
            yK1[2] += C1[i] * this->_h * this->_y[3] + this->_h * this->_h * s2y_k1;
            yK1[3] += this->_h * s2ys;
            K1[i] = _rhsF(this->_x + C1[i] * this->_h, yK1);

            for (int j = 0; j < i; j++) {
                s2x_k2 += A21[i][j] * K1[j];
                s2y_k2 += A22[i][j] * K2[j];
                s2xs += B21[i][j] * K1[j];
            }
            s2x_k2 += A21[i][i] * K1[i];
            s2xs += B21[i][i] * K1[i];
            Eigen::VectorX<Scalar> yK2(this->_y);
            yK2[0] += C2[i] * this->_h * this->_y[1] + this->_h * this->_h * s2x_k2;
            yK2[1] += this->_h * s2xs;
            yK2[2] += C2[i] * this->_h * this->_y[3] + this->_h * this->_h * s2y_k2;

            K2[i] = _rhsG(this->_x + C2[i] * this->_h, yK2);
            if (this->_logging) this->_logger->addCalls(1);
        }
        Scalar tempK1(_K1Last), tempK2(_K2Last);
        _K1Last = K1[6];
        _K2Last = K2[6];
        Eigen::VectorX<Scalar> Y(this->_y), Yi;
        Y[0] += this->_h * this->_y[1];
        Y[2] += this->_h * this->_y[3];
        Yi = Y;
        for (int i = 0; i < 6; i++) {
            Y[0] += this->_h * this->_h * D1[i] * K1[i];
            Y[1] += this->_h * Q1[i] * K1[i];
            Y[2] += this->_h * this->_h * D2[i] * K2[i];
            Y[3] += this->_h * Q2[i] * K2[i];
        }
        for (int i = 0; i < 7; i++) {
            Yi[0] += this->_h * this->_h * Di1[i] * K1[i];
            Yi[1] += this->_h * Qi1[i] * K1[i];
            Yi[2] += this->_h * this->_h * Di2[i] * K2[i];
            Yi[3] += this->_h * Qi2[i] * K2[i];
        }
        Scalar truncationError = euclidNorm((Yi - Y).eval());

        if (truncationError < this->_eps) {
            //It is important to perform solveOnStep before modifying y so we compare the same Cauchy problems
            Eigen::VectorX<Scalar> y_test;
            if (this->_logging) {
                this->_logger->setStepTaken(true);
                this->_logging = false;
                _firstIter = true;
                y_test = this->solveOnStep();
                this->_logging = true;
            }

            this->_y = Y;
            this->_x += this->_h;

            //This logging step is only done to accepted steps
            if (this->_logging) {
                this->_logger->logErrorEstimation(this->_x, this->_y, y_test, truncationError);
            }
        } else {
            _K1Last = tempK1;
            _K2Last = tempK2;
        }
        //This logging step is done even for steps that were not accepted due to precision requirements
        if (this->_logging) {
            this->_logger->logStep(this->_x, this->_y, this->_h);
        }

        this->_h = std::min(0.9 * this->_h * std::pow((this->_eps / truncationError), 0.2), this->_xf - this->_x);
    }
public:
    SRK64Solver(const Scalar &x, const Scalar &xf, const Scalar &eps, Eigen::VectorX<Scalar> y,
                std::function<Scalar(Scalar, Eigen::VectorX<Scalar>)> rhsF,
                std::function<Scalar(Scalar, Eigen::VectorX<Scalar>)> rhsG,
                bool logging = false,
                const std::string &steps = "", const std::string &stepSizes = "",
                const std::string &stepErrorEstimation = ""){
        this->_x = x;
        this->_xf = xf;
        this->_eps = eps;
        this->_y = y;
        this->_rhsF = rhsF;
        this->_rhsG = rhsG;
        this->_rhs = [&](Scalar x, Eigen::VectorX<Scalar> y){
            Eigen::VectorXd result(y.size());
            result[0] = y[1];
            result[1] = rhsF(x, y);
            result[2] = y[3];
            result[3] = rhsG(x, y);
            return result;
        };
        this->findStartingStep();
        if (logging) {
            this->_logging = true;
            this->_logger = new(RKLogger<Scalar>)(steps, stepSizes, stepErrorEstimation);
            this->_logger->addCalls(1);
        }
    }
};

#endif //EMBEDDED_METHODS_CPP_SRK64SOLVER_H
