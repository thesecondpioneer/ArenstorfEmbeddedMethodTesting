#ifndef EMBEDDED_METHODS_CPP_RK45SOLVER_H
#define EMBEDDED_METHODS_CPP_RK45SOLVER_H

#include "EmbeddedRungeKutta.h"

template<typename Scalar>
class RK45Solver : public EmbeddedRungeKutta<Scalar> {
protected:
    using EmbeddedRungeKutta<Scalar>::EmbeddedRungeKutta;
    //method coefficients
    const Scalar a2 = 1.0 / 4, a3 = 3.0 / 8, a4 = 12.0 / 13, a5 = 1.0, a6 = 1.0 / 2,
            b21 = 1.0 / 4, b31 = 3.0 / 32, b41 = 1932.0 / 2197, b51 = 439.0 / 216, b61 = -8.0 / 27,
            b32 = 9.0 / 32, b42 = -7200.0 / 2197, b52 = -8.0, b62 = 2.0,
            b43 = 7296.0 / 2197, b53 = 3680.0 / 513, b63 = -3544.0 / 2565,
            b54 = -845.0 / 4104, b64 = 1859.0 / 4104,
            b65 = -11.0 / 40,
            c1 = 25.0 / 216, c3 = 1408.0 / 2565, c4 = 2197.0 / 4104, c5 = -1.0 / 5,
            ch1 = 16.0 / 135, ch3 = 6656.0 / 12825, ch4 = 28561.0 / 56430, ch5 = -9.0 / 50, ch6 = 2.0 / 55,
            ct1 = -1.0 / 360, ct3 = 128.0 / 4275, ct4 = 2197.0 / 75240, ct5 = -1.0 / 50, ct6 = -2.0 / 55;

    void makeStep() override {
        //Calculate the correction functions for the step
        Eigen::VectorXd k1 = this->_h * this->_rhs(this->_x, this->_y),
                k2 = this->_h * this->_rhs(this->_x + a2 * this->_h, this->_y + b21 * k1),
                k3 = this->_h * this->_rhs(this->_x + a3 * this->_h, this->_y + b31 * k1 + b32 * k2),
                k4 = this->_h * this->_rhs(this->_x + a4 * this->_h, this->_y + b41 * k1 + b42 * k2 + b43 * k3),
                k5 = this->_h * this->_rhs(this->_x + a5 * this->_h, this->_y + b51 * k1 + b52 * k2 + b53 * k3 + b54 * k4),
                k6 = this->_h *
                     this->_rhs(this->_x + a6 * this->_h, this->_y + b61 * k1 + b62 * k2 + b63 * k3 + b64 * k4 + b65 * k5);
        if (this->_logging) this->_logger->addCalls(6);

        Scalar truncationError = euclidNorm((ct1 * k1 + ct3 * k3 + ct4 * k4 + ct5 * k5 + ct6 * k6).eval());
        if (truncationError < this->_eps) {
            //It is important to perform solveOnStep before modifying y so we compare the same Cauchy problems
            Eigen::VectorX<Scalar> y_test;
            if (this->_logging) {
                this->_logger->setStepTaken(true);
                this->_logging = false;
                y_test = this->solveOnStep();
                this->_logging = true;
            }

            this->_y += c1 * k1 + c3 * k3 + c4 * k4 + c5 * k5;
            this->_x += this->_h;

            //This logging step is only done to accepted steps
            if (this -> _logging) {
                this->_logger->logErrorEstimation(this->_x, this->_y, y_test, truncationError);
            }
        }
        //This logging step is done even for steps that were not accepted due to precision requirements
        if (this->_logging){
            this->_logger->logStep(this->_x, this->_y, this->_h);
        }

        this->_h = std::min(0.9 * this->_h * std::pow((this->_eps / truncationError), 0.2), this->_xf - this->_x);
    }
};
#endif //EMBEDDED_METHODS_CPP_RK45SOLVER_H
