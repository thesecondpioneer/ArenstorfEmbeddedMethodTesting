#ifndef EMBEDDED_METHODS_CPP_RKUTILS_H
#define EMBEDDED_METHODS_CPP_RKUTILS_H

#include "Eigen/Core"
#include <functional>
#include <fstream>

//Returns he Euclidean norm of an Eigen vector
template<typename Scalar, typename std::enable_if<std::is_arithmetic_v<Scalar>, bool>::type = true>
inline Scalar euclidNorm(const Eigen::VectorX<Scalar> &v){
    return std::sqrt(v.array().square().sum());
}

//Returns the infinity norm of an Eigen vector
template<typename Scalar, typename std::enable_if<std::is_arithmetic_v<Scalar>, bool>::type = true>
inline Scalar infNorm(const Eigen::VectorX<Scalar> &v){
    return v.array().abs().maxCoeff();
}

//Returns one Euler step approximation of the solution y with step size h
template<typename Scalar, typename std::enable_if<std::is_arithmetic_v<Scalar>, bool>::type = true>
inline Eigen::VectorX<Scalar> eulerStep(std::function<Eigen::VectorX<Scalar>(Scalar, Eigen::VectorX<Scalar>)> f, const Scalar &x0,
                        const Eigen::VectorX<Scalar> &y0, Scalar h){
    return y0 + h * f(x0, y0);
};

#endif //EMBEDDED_METHODS_CPP_RKUTILS_H
