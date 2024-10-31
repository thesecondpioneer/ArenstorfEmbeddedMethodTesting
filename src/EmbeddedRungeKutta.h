#ifndef EMBEDDED_METHODS_CPP_EMBEDDEDRUNGEKUTTA_H
#define EMBEDDED_METHODS_CPP_EMBEDDEDRUNGEKUTTA_H

#include "RKUtils.h"
#include "RKLogger.h"
#include <cmath>
#include <iostream>
#include <iomanip>

template<typename Scalar>
//Base class for all embedded Runge-Kutta methods
class EmbeddedRungeKutta {
    static_assert(std::is_floating_point_v<Scalar>);

protected:
    RKLogger<Scalar> *_logger;

    //current value of the independent variable of the ODE system
    Scalar _x;

    //final value of x, to which we solve the ODE system
    Scalar _xf;

    //precision
    Scalar _eps;

    //maximal precision
    Scalar _epsm = 1e-15;

    //current step size
    Scalar _h;

    //current approximation of the solution
    Eigen::VectorX<Scalar> _y;

    //ODE system right side functor
    std::function<Eigen::VectorX<Scalar>(Scalar, Eigen::VectorX<Scalar>)> _rhs;

    bool _logging = false;

    void findStartingStep() {
        Eigen::VectorX<Scalar> y_h1;
        Scalar delta = std::pow(1.0 / std::max(std::abs(_x), std::abs(_xf)), 5) + std::pow(euclidNorm(_rhs(_x, _y)), 5);
        _h = std::pow(_eps / delta, 0.2);
        y_h1 = eulerStep(_rhs, _x, _y, _h);
        delta = std::pow(1.0 / std::max(std::abs(_x), std::abs(_xf)), 5) + std::pow(euclidNorm(_rhs(_x + _h, _y)), 5);
        _h = std::pow(_eps / delta, 0.2);
    }

    virtual void makeStep() = 0;
public:
    EmbeddedRungeKutta() = default;

    EmbeddedRungeKutta(const Scalar &x, const Scalar &xf, const Scalar &eps, Eigen::VectorX<Scalar> y,
                       std::function<Eigen::VectorX<Scalar>(Scalar, Eigen::VectorX<Scalar>)> rhs,
                       bool logging = false,
                       const std::string &steps = "", const std::string &stepSizes = "",
                       const std::string &stepErrorEstimation = "")
            : _x(x), _xf(xf), _y(y), _eps(eps), _rhs(rhs) {
        this->findStartingStep();
        if (logging) {
            _logging = true;
            _logger = new(RKLogger<Scalar>)(steps, stepSizes, stepErrorEstimation);
        }
    }

    void setMaxPrecision(const Scalar &epsm) {
        _epsm = epsm;
    }

    std::vector<Eigen::VectorX<Scalar>> solve() {
        std::vector<Eigen::VectorX<Scalar>> result;
        while (_x < _xf) {
            makeStep();
        }
        if(_logging) std::cout << "Solution finished with " << _logger->getCalls() << " calls to the right-hand side\n";
        return result;
    }

    //Method for solving the ODE with x as the start and x + h as the end at max precision.
    //This is needed for the purposes of evaluating the quality of truncation error estimation
    Eigen::VectorX<Scalar> solveOnStep(){
        //save the solver state before the method call
        Scalar x = _x, eps = _eps, xf = _xf, h = _h;
        Eigen::VectorX<Scalar> y = _y, result;

        //change the solver parameters
        _eps = _epsm;
        _xf = _x + _h;
        while (_x < _xf) {
            makeStep();
        }
        result = _y;

        _x = x;
        _eps = eps;
        _xf = xf;
        _h = h;
        _y = y;
        return result;
    }
};

#endif //EMBEDDED_METHODS_CPP_EMBEDDEDRUNGEKUTTA_H
