//
//  Matrix.h
//  MatrixClass
//
//  Created by Li Xinrui on 10/1/14.
//  Copyright (c) 2014 Li Xinrui. All rights reserved.
//

#ifndef __MatrixClass__Matrix__
#define __MatrixClass__Matrix__
#define TYPE double
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <cstdio>
#include <ctime>
#include <cassert>

using std::vector;
using std::cout;
using std::endl;


class Matrix{
//  Members
protected:
    int mNumRows = 0;
    int mNumCols = 0;
    TYPE **data;
public:
// Constructors
    Matrix();
    Matrix(int row, int col);
    Matrix(int dim);
    Matrix(const Matrix& m);
    
// Deconstructor
    ~Matrix(){
        for (int i = 0; i != this -> mNumCols; i++) {
            delete [] data[i];
        }
        delete [] data;
    }
    
public:
    // Operators
    TYPE* operator[](int index);
    const TYPE* operator[](int index) const;
    TYPE& operator()(int rowIndex, int colIndex);
    const TYPE& operator() (int rowIndex, int colIndex) const;
    // Members
    int GetNumRows(){
        return mNumRows;
    }
    int GetNumRows() const{
        return mNumRows;
    }
    int GetNumCols() const{
        return mNumCols;
    }
    int GetNumCols(){
        return mNumCols;
    }
    void ResizeRow(int row){
        mNumRows = row;
    }
    void ResizeCol(int col) {
        mNumCols = col;
    }

    // Operations
    Matrix& operator+=(const Matrix& rhs);
    Matrix& operator-=(const Matrix& rhs);
    Matrix& operator*=(const TYPE& rhs);
    Matrix& operator/=(const TYPE& rhs);
    Matrix& operator=(const TYPE* rhs);
    Matrix& operator=(const vector<TYPE>& rhs);
    Matrix& operator=(const Matrix& rhs);


public:
    Matrix& rowSwap(int i, int j);
    Matrix& colSwap(int i, int j);
    Matrix& transpose();
    double Norm2Vec();
    int transpose(Matrix& res);
    bool QSquare();
    void SetIdentity(){
        if( !QSquare() ){
            return;
        }
        for (int i = 0; i != this -> mNumCols; i++) {
            data[i][i] = 1;
        }
    }

public:
    friend std::ostream& operator<<(std::ostream& out, const Matrix& m);
    friend std::ofstream& operator<<(std::ofstream& out, Matrix& m);
    friend std::ostream& latex(std::ostream& out, const Matrix& m, int line);

};
// Operations
Matrix operator+(const Matrix& m1, const Matrix& m2);
Matrix operator-(const Matrix& m1, const Matrix& m2);
Matrix operator*(const Matrix& m1, const Matrix& m2);
Matrix operator*(const Matrix& m1, const TYPE& num);
Matrix operator*(const TYPE& num, const Matrix& m1);
Matrix operator/(const Matrix& m1, const TYPE& num);



// Algorithms
int LUDecomp(const Matrix& m, Matrix& l, Matrix& u);
int LUPivotDecomp(const Matrix& m, Matrix& p, Matrix& l, Matrix& u);
int CholeskyDecomp(const Matrix& m, Matrix& l);
int LDLDecomp(const Matrix& m, Matrix &l, Matrix& d);

int LUSolve(const Matrix& m,  const Matrix& rhs, Matrix& res);
int LUPivotSolve(const Matrix& m, const Matrix& rhs, Matrix& res);
int CholeskySolve(const Matrix&m, const Matrix& rhs, Matrix& res);
int LDLSolve(const Matrix& m, const Matrix& rhs, Matrix& res);

int lowerSolve(const Matrix& l,const Matrix& rhs, Matrix& res);
int upperSolve(const Matrix& u,const Matrix& rhs, Matrix& res);


class mVector : public Matrix {
public:
    using Matrix::operator=;
    using Matrix::operator/=;
    using Matrix::operator*=;
    using Matrix::operator+=;
    using Matrix::operator-=;
public:
    int dim(){
        return this -> mNumRows;
    }
    int dim() const{
        return this -> mNumRows;
    }
public:
    TYPE& operator[](int index){
        return data[0][index];
    }
    const TYPE& operator[](int index) const{
        return data[0][index];
    }
public:
    mVector() {}
    mVector(int dim){
        mNumRows = dim;
        mNumCols = 1;
        data = new TYPE*[1];
        data[0] = new TYPE[dim];
        for (int i = 0; i != dim; i++) {
            data[0][i] = 0;
        }
    }
    mVector(const mVector& rhs){
        mNumRows = rhs.dim();
        mNumCols = 1;
        data = new TYPE*[1];
        data[0] = new TYPE[rhs.dim()];
        for (int i = 0; i != rhs.dim(); i++) {
            data[0][i] = rhs[i];
        }
    }
    double NormInf();
    double NormInf(int& position);
    double Norm_2();
    double Norm_1();
};
mVector operator+(const mVector& m1, const mVector& m2);
mVector operator-(const mVector& m1, const mVector& m2);
mVector operator*(const Matrix& m1, const mVector& m2);
mVector operator*(const mVector& m1, const TYPE& num);
mVector operator*(const TYPE& num, const mVector& m1);
mVector operator/(const mVector& m1, const TYPE& num);



double InnerProduct(mVector& v1, mVector& v2);








#endif /* defined(__MatrixClass__Matrix__) */
