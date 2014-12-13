//
//  IterativeMethods.cpp
//  Numerical_C4
//
//  Created by Li Xinrui on 11/26/14.
//  Copyright (c) 2014 Li Xinrui. All rights reserved.
//

#include "IterativeMethods.h"
int IterativeMethods::JacobiSolve(Matrix& A, mVector& b, mVector& guess, mVector& res, double precision){
    assert(A.GetNumCols() == A.GetNumRows());
    int n = A.GetNumCols();
    Matrix inverseD = GenerateInverseD(A);
    Matrix L = GenerateL(A);
    Matrix U = GenerateU(A);
    Matrix B = inverseD * (L + U);
    mVector g(n);
    g = inverseD * b;
    mVector xPre(n);
    mVector xPost(n);
    xPost = guess;
    mVector diff(n);
    int count = 0;
    do {
        xPre = xPost;
        xPost = B * xPre + g;
        diff = xPre - xPost;
        count++;
    } while ((diff).NormInf() >= precision);
    res = xPost;
    return count;
}

int IterativeMethods::GaussSolve(Matrix &A, mVector &b, mVector &guess, mVector &res, double precision){
    return SORSolve(A, b, guess, res, 1, precision);
}

int IterativeMethods::GaussSolveNorm2(Matrix &A, mVector &b, mVector &guess, mVector &res, double precision){
    return SORSolveNorm2(A, b, guess, res, 1, precision);
}

int IterativeMethods::SORSolve(Matrix &A, mVector &b, mVector &guess, mVector &res, double w, double precision){
    assert(A.GetNumCols() == A.GetNumRows());
    int n = A.GetNumRows();
    Matrix inverseD = GenerateInverseD(A);
    Matrix L = GenerateL(A);
    Matrix U = GenerateU(A);
    Matrix B = inverseD * (L + U);
    mVector g(n);
    g = inverseD * b;
    mVector xPre(n), xPost(n);
    xPost = guess;
    mVector diff(n);
    
    int count = 0;
    do {
        xPre = xPost;
        for (int i = 0; i != n; i++) {
            double s1 = 0;
            for (int j = 0; j <= i - 1; j++) {
                s1 += xPost[j] * B(i,j);
            }
            double s2 = 0;
            for (int j = i + 1; j != n; j++) {
                s2 += xPre[j] * B(i,j);
            }
            xPost[i] = (1 - w) * xPre[i] + w * (s1 + s2 + g[i]);
        }
        diff = xPost - xPre;
        count++;
    } while (diff.NormInf() >= precision);
    res = xPost;
    return count;
}

int IterativeMethods::DirectGaussSolve(Matrix &A, mVector &b, mVector &guess, mVector &res, double precision){
    assert(A.GetNumCols() == A.GetNumRows());
    int n = A.GetNumRows();
    mVector xPre(n), xPost(n);
    xPost = guess;
    mVector diff(n);
    int count = 0;
    
    
    return count;
}

int IterativeMethods::SORSolveNorm2(Matrix &A, mVector &b, mVector &guess, mVector &res, double w, double precision){
    assert(A.GetNumCols() == A.GetNumRows());
    int n = A.GetNumRows();
    mVector xPre(n), xPost(n);
    xPost = guess;
    mVector diff(n);
    mVector g(n);
    for (int i = 0; i != n; i++) {
        g[i] = b[i] / A(i,i);
    }
    int count = 0;
    do {
        xPre = xPost;
        for (int i = 0; i != n; i++) {
            double s1 = 0;
            for (int j = 0; j <= i - 1; j++) {
                double temp = - A(i,j)/ A(i,i);
                s1 += xPost[j] * temp;
            }
            double s2 = 0;
            for (int j = i + 1; j != n; j++) {
                double temp = - A(i,j) / A(i,i);
                s2 += xPre[j] * temp;
            }
            xPost[i] = (1 - w) * xPre[i] + w * (s1 + s2 + g[i]);
        }
        diff = xPost - xPre;
        count++;
        if (count % 10 ==0) {
            std::cout << count << std::endl;
        }
    } while (diff.Norm_2() >= precision);
    res = xPost;
    return count;
}

Matrix IterativeMethods::GenerateInverseD(Matrix &A) {
    Matrix temp(A.GetNumCols());
    for (int i = 0; i != A.GetNumCols(); i++) {
        assert(A(i,i) != 0);
        temp(i,i) = 1.0 / A(i,i);
    }
    return temp;
}
Matrix IterativeMethods::GenerateD(Matrix &A) {
    Matrix temp(A.GetNumCols());
    for (int i = 0; i != A.GetNumCols(); i++) {
        temp(i,i) = A(i,i);
    }
    return temp;
}

Matrix IterativeMethods::GenerateL(Matrix &A){
    int n = A.GetNumCols();
    Matrix temp(n,n);
    for (int i = 0; i != n ; i++) {
        for (int j = 0; j != i; j++) {
            temp(i,j) = - A(i,j);
        }
    }
    return temp;
}

Matrix IterativeMethods::GenerateU(Matrix &A){
    int n = A.GetNumCols();
    Matrix temp(n,n);
    for (int i = 0; i != n ; i++) {
        for (int j = 0; j != i; j++) {
            temp(j,i) = - A(j,i);
        }
    }
    return temp;
}