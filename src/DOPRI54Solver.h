#ifndef EMBEDDED_METHODS_CPP_DOPRI54SOLVER_H
#define EMBEDDED_METHODS_CPP_DOPRI54SOLVER_H

#include "EmbeddedRungeKutta.h"

template<typename Scalar>
class DOPRI54Solver : public EmbeddedRungeKutta<Scalar> {
protected:
    bool _firstIter = true;
    Eigen::VectorX<Scalar> _rhsLast;
    //method coefficients
    const Scalar a2 = 1.0 / 5, a3 = 3.0 / 10, a4 = 4.0 / 5, a5 = 8.0 / 9, a6 = 1.0, a7 = 1.0,
            b21 = 1.0 / 5, b31 = 3.0 / 40, b41 = 44.0 / 45, b51 = 19372.0 / 6561, b61 = 9017.0 / 3168, b71 = 35.0 / 384,
            b32 = 9.0 / 40, b42 = -56.0 / 15, b52 = -25360.0 / 2187, b62 = -355.0 / 33, b72 = 0,
            b43 = 32.0 / 9, b53 = 64448.0 / 6561, b63 = 46732.0 / 5247, b73 = 500.0 / 1113,
            b54 = -212.0 / 729, b64 = 49.0 / 176, b74 = 125.0 / 192,
            b65 = -5103.0 / 18656, b75 = -2187.0 / 6784, b76 = 11.0 / 84,
            c1 = 35.0 / 384, c3 = 500.0 / 1113, c4 = 125.0 / 192, c5 = -2187.0 / 6784, c6 = 11.0 / 84,
            ch1 = 5179.0 / 57600, ch3 = 7571.0 / 16695, ch4 = 393.0 / 640, ch5 = -92097.0 / 339200,
            ch6 = 187.0 / 2100, ch7 = 1.0 / 40,
            ct1 = 71.0 / 57600, ct3 = -71.0 / 16695, ct4 = 71.0 / 1920, ct5 = -17253.0 / 339200,
            ct6 = 22.0 / 525, ct7 = -1.0 / 40;

    void makeStep() override {
        //Calculate the correction functions for the step
        Eigen::VectorXd k1 = _firstIter ? this->_h * this->_rhs(this->_x, this->_y) : this->_h * _rhsLast,       //FSAL
        k2 = this->_h * this->_rhs(this->_x + a2 * this->_h, this->_y + b21 * k1),
                k3 = this->_h * this->_rhs(this->_x + a3 * this->_h, this->_y + b31 * k1 + b32 * k2),
                k4 = this->_h * this->_rhs(this->_x + a4 * this->_h, this->_y + b41 * k1 + b42 * k2 + b43 * k3),
                k5 = this->_h *
                     this->_rhs(this->_x + a5 * this->_h, this->_y + b51 * k1 + b52 * k2 + b53 * k3 + b54 * k4),
                k6 = this->_h *
                     this->_rhs(this->_x + a6 * this->_h,
                                this->_y + b61 * k1 + b62 * k2 + b63 * k3 + b64 * k4 + b65 * k5);
        Eigen::VectorXd temp(_rhsLast); //In case the step isn't accepted, we roll back to previous _rhsLast
        _rhsLast = this->_rhs(this->_x + a7 * this->_h,
                              this->_y + b71 * k1 + b72 * k2 + b73 * k3 + b74 * k4 + b75 * k5 + b76 * k6);
        Eigen::VectorXd k7 = this->_h * _rhsLast;
        if (this->_logging) this->_logger->addCalls(6);
        _firstIter = false;

        Scalar truncationError = euclidNorm((ct1 * k1 + ct3 * k3 + ct4 * k4 + ct5 * k5 + ct6 * k6 + ct7 * k7).eval());
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

            this->_y += c1 * k1 + c3 * k3 + c4 * k4 + c5 * k5 + c6 * k6;
            this->_x += this->_h;

            //This logging step is only done to accepted steps
            if (this->_logging) {
                this->_logger->logErrorEstimation(this->_x, this->_y, y_test, truncationError);
            }
        } else {
            _rhsLast = temp;
        }
        //This logging step is done even for steps that were not accepted due to precision requirements
        if (this->_logging) {
            this->_logger->logStep(this->_x, this->_y, this->_h);
        }

        this->_h = std::min(0.9 * this->_h * std::pow((this->_eps / truncationError), 0.2), this->_xf - this->_x);
    }
public:
    DOPRI54Solver(const Scalar &x, const Scalar &xf, const Scalar &eps, Eigen::VectorX<Scalar> y,
                  std::function<Eigen::VectorX<Scalar>(Scalar, Eigen::VectorX<Scalar>)> rhs,
                  bool logging = false,
                  const std::string &steps = "", const std::string &stepSizes = "",
                  const std::string &stepErrorEstimation = "") :
            EmbeddedRungeKutta<Scalar>(x, xf, eps, y, rhs, logging, steps, stepSizes, stepErrorEstimation) {
        this->_logger->addCalls(1);
    }
};

#endif //EMBEDDED_METHODS_CPP_DOPRI54SOLVER_H
