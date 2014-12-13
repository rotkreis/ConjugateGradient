//
//  main.cpp
//  CongujateGradient
//
//  Created by Li Xinrui on 12/13/14.
//  Copyright (c) 2014 Li Xinrui. All rights reserved.
//

#include <iostream>
#include "Matrix.h"
#include "CG.h"

static Matrix CreateHilbertMatrix(int dim) {
    Matrix temp(dim,dim);
    for (int i = 0; i != dim; i++) {
        for (int j = 0 ; j != dim; j++) {
            temp(i, j) = 1.0 / (i + j + 3);
        }
    }
    return temp;
}

static mVector CreateHilbertRHS(int dim) {
    mVector temp(dim);
    for (int i = 0; i != dim; i++) {
        for (int j = 0; j != dim; j++) {
            temp[i] += 1.0 / ((i + j + 3) * 3.0);
        }
    }
    return temp;
}

int main(int argc, const char * argv[]) {
//    mVector v(2);
//    v[0] = 1;
//    v[1] = 2;
//    Matrix A(2,2);
//    A = {1,2,3,4};
//    mVector p = A * v;
//    
//    mVector x(4);
//    x = {.33,.33,.33,.33};
//    Matrix h = CreateHilbertMatrix(4);
//    mVector r = CreateHilbertRHS(4);
//    std::cout << h << r << std::endl;
//    std::cout << r - h * x;
    int dim = 10;
    mVector res(dim);
    mVector guess(dim);
    for (int i = 0; i != dim; i++) {
        guess[i] = 1;
    }
    Matrix A = CreateHilbertMatrix(dim);
    mVector b = CreateHilbertRHS(dim);
    int kMax = 10000;
    std::cout << CGSolver(A, b, res, guess, 10E-14, kMax) << endl;
    std::cout << res;
    return 0;
}
