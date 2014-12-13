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

int main(int argc, const char * argv[]) {
    mVector v(2);
    v[0] = 1;
    v[1] = 2;
    Matrix A(2,2);
    A = {1,2,3,4};
    mVector p = A * v;
    std::cout << A << v << std::endl;
    std::cout << p << std::endl;

    return 0;
}
