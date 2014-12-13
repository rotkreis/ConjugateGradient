//
//  IterativeMethods.h
//  Numerical_C4
//
//  Created by Li Xinrui on 11/26/14.
//  Copyright (c) 2014 Li Xinrui. All rights reserved.
//

#ifndef __Numerical_C4__IterativeMethods__
#define __Numerical_C4__IterativeMethods__

#include <stdio.h>
#include "Matrix.h"

namespace IterativeMethods {
    int JacobiSolve(Matrix& A, mVector& b, mVector& guess, mVector& res, double precision);
    int GaussSolve(Matrix& A, mVector& b, mVector& guess, mVector& res, double precision);
    int GaussSolveNorm2(Matrix& A, mVector& b, mVector& guess, mVector& res, double precision);
    int SORSolve(Matrix& A, mVector& b, mVector& guess, mVector& res, double w, double precision);
    int SORSolveNorm2(Matrix& A, mVector& b, mVector& guess, mVector& res, double w, double precision);
    int DirectGaussSolve(Matrix& A, mVector& b, mVector& guess, mVector& res, double precision);
    Matrix GenerateInverseD(Matrix& A);
    Matrix GenerateD(Matrix& A);
    Matrix GenerateL(Matrix& A);
    Matrix GenerateU(Matrix& A);
}

#endif /* defined(__Numerical_C4__IterativeMethods__) */
