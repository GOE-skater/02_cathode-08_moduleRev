#pragma once
#include <complex>
#include <cmath>
#include <stdexcept>

using namespace std;

struct CmpMat2 {

    complex<double> A, B, C, D;

    complex<double> det() const {
        return A * D - B * C;
    }

    CmpMat2 inv() const {
        complex<double> d = det();
        return {
             D / d,
            -B / d,
            -C / d,
             A / d
        };
    }
};

CmpMat2 operator*(const CmpMat2& X, const CmpMat2& Y) {
    return {
        X.A * Y.A + X.B * Y.C,
        X.A * Y.B + X.B * Y.D,
        X.C * Y.A + X.D * Y.C,
        X.C * Y.B + X.D * Y.D
    };
}
