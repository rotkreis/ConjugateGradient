//
//  CG.h
//  CongujateGradient
//
//  Created by Li Xinrui on 12/13/14.
//  Copyright (c) 2014 Li Xinrui. All rights reserved.
//

#ifndef __CongujateGradient__CG__
#define __CongujateGradient__CG__

#include <stdio.h>
#include "Matrix.h"
int CGSolver(Matrix& A, mVector& b, mVector& res, mVector& guess, double precision, int max);


#endif /* defined(__CongujateGradient__CG__) */
