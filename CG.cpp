//
//  CG.cpp
//  CongujateGradient
//
//  Created by Li Xinrui on 12/13/14.
//  Copyright (c) 2014 Li Xinrui. All rights reserved.
//

#include "CG.h"
int CGSolver(Matrix& A, mVector& b, mVector& res, mVector& guess, double precision, int kMax){
    mVector x = guess;
    int k = 0;
    mVector r = b - A * x;
    double rho = InnerProduct(r, r);
    double rhoSlash = 0;
    mVector p(b.dim());
    while (std::sqrt(rho) > precision * b.Norm_2() && k < kMax) {
        k++;
        if (k == 1) {
            p = r;
        }
        else {
            double beta = rho / rhoSlash;
            p *= beta;
            p += r;
        }
        mVector w = A * p;
        double alpha = rho / (InnerProduct(p, w));
        x += p * alpha;
        r -= w * alpha;
        rhoSlash = rho;
        rho = InnerProduct(r, r);
    }
    res = x;
    return k;
}
