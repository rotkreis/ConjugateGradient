//
//  main.cpp
//  CongujateGradient
//
//  Created by Li Xinrui on 12/13/14.
//  Copyright (c) 2014 Li Xinrui. All rights reserved.
//

#include <iostream>
#include <cmath>
#include "Matrix.h"
#include "CG.h"
#include "IterativeMethods.h"

double f(double x, double y){
    return sin(x * y);
}
double boundary(double x, double y){
    return x*x + y*y;
}

void problem1(int N){
    int nNodes = N - 1;
    double h = 1.0 / N;
    Matrix A(nNodes * nNodes, nNodes * nNodes);
    mVector rhs(nNodes * nNodes);
    for (int i = 0; i != nNodes; i++) {
        for (int j = 0; j != nNodes; j++) {
            rhs[i * nNodes + j] = h * h / 4.0 * f((i + 1) * h, (j + 1) * h);
            if (i == 0) {
                rhs[i * nNodes + j] += boundary(0, (j + 1) * h) / 4;
            }
            if (i == nNodes - 1) {
                rhs[i * nNodes + j] += boundary(1, (j + 1) * h) / 4;
            }
            if (j == 0) {
                rhs[i * nNodes + 1] += boundary((i + 1) * h, 0) / 4;
            }
            if (j == nNodes - 1) {
                rhs[i * nNodes + 1] += boundary((i + 1) * h, 1) / 4;
            }
        }
    }
    for (int i = 0; i != nNodes; i++) {
        for (int j = 0; j != nNodes; j++) {
            A(i * nNodes + j, i * nNodes + j) = 1 + h * h / 4;
            if (i != 0) {
                A(i * nNodes + j, (i - 1) * nNodes + j) = -1.0 / 4;
            }
            if (i != nNodes - 1) {
                A(i * nNodes + j, (i + 1) * nNodes + j) = -1.0 / 4;
            }
            if (j != 0) {
                A(i * nNodes + j, i * nNodes + j - 1) = -1.0 / 4;
            }
            if (j != nNodes - 1) {
                A(i * nNodes + j, i * nNodes + j + 1) = -1.0 / 4;
            }
        }
    }
    mVector guess(nNodes * nNodes);
    for (int i = 0; i != guess.dim(); i++) {
        guess[i] = 5;
    }
    mVector res(nNodes * nNodes);
    std::cout << CGSolver(A, rhs, res, guess, 10E-10, 100) << std::endl;
    std::cout << res;
}

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
    //-------------- Problem 1
//    problem1(20);
  
    
    // -------------- Problem 2
//    int dim = 20;
//    mVector res(dim);
//    mVector guess(dim);
//    for (int i = 0; i != dim; i++) {
//        guess[i] = 1;
//    }
//    Matrix A = CreateHilbertMatrix(dim);
//    mVector b = CreateHilbertRHS(dim);
//    std::cout << CGSolver(A, b, res, guess, 10E-14, 100) << endl;
//    std::cout << res;
    

    // --------------Problem 3
    Matrix prob3(5,5);
    prob3 = {10,1,2,3,4,1,9,-1,2,-3,2,-1,7,3,-5,3,2,3,12,-1,4,-3,-5,-1,15};
    mVector rhs3(5);
    rhs3 = {12,-27,14,-17,12};
    mVector res3(5);
    mVector guess3(5);
    mVector Jacobi(5);
    mVector GS(5);
    double precision = 10E-5;
    std::cout <<"CG: " <<CGSolver(prob3, rhs3, res3, guess3, precision, 100) << std::endl;
    std::cout << res3;
    std::cout <<"Jacobi: "<< IterativeMethods::JacobiSolve(prob3, rhs3, guess3, Jacobi, precision) << std::endl;
    std::cout << Jacobi;
    std::cout << "GS: " << IterativeMethods::GaussSolve(prob3, rhs3, guess3, GS, precision) << std::endl;
    std::cout << GS;
    return 0;
}
