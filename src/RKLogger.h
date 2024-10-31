#ifndef EMBEDDED_METHODS_CPP_RKLOGGER_H
#define EMBEDDED_METHODS_CPP_RKLOGGER_H

#include <fstream>
#include <iomanip>
#include "Eigen/Core"
#include "RKUtils.h"
#include <functional>

template<typename Scalar>
class RKLogger {
protected:
    //The amount of calls to the right hand side of the ODE that the method makes
    int _calls = 0;

    //Check if the step size was appropriate this step or if it has to be changed due to precision restrictions
    bool _stepTaken = false;

    std::ofstream _steps;

    std::ofstream _stepSizes;

    std::ofstream _stepErrorEstimation;

public:
    RKLogger() = default;

    RKLogger(const std::string &steps, const std::string &stepSizes, const std::string &stepErrorEstimation) :
            _steps(steps),
            _stepSizes(stepSizes),
            _stepErrorEstimation(stepErrorEstimation) {};

    void setStepLoggingOutput(const std::string &steps) {
        _steps = std::ofstream(steps);
    }

    void setStepSizeLoggingOutput(const std::string &stepSizes) {
        _stepSizes = std::ofstream(stepSizes);
    }

    void setStepErrorEstimationLoggingOutput(const std::string &steps) {
        _steps = std::ofstream(steps);
    }

    void setStepTaken(const bool &b) {
        _stepTaken = b;
    }

    void logStep(const Scalar &x, const Eigen::VectorX<Scalar> &y, const Scalar &h) {
        if (_stepTaken) {
            for (auto i: y) {
                _steps << std::setprecision(16) << std::fixed << i << ',';
            }
            _steps << std::setprecision(16) << std::fixed << x << std::endl;
        }
        _stepSizes << std::setprecision(16) << std::fixed << h << ',' << x << std::endl;
        _stepTaken = false;
    }

    void logErrorEstimation(const Scalar &x, const Eigen::VectorX<Scalar> &y, const Eigen::VectorX<Scalar> &y_test,
                            const Scalar &truncationError){
        _stepErrorEstimation << std::setprecision(16) << std::fixed << x << ','
                             << std::abs(truncationError - euclidNorm<Scalar>(y - y_test)) << std::endl;
    };

    void addCalls(int c){
        _calls += c;
    }

    int getCalls(){
        return _calls;
    }
};

#endif //EMBEDDED_METHODS_CPP_RKLOGGER_H
