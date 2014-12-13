//
//  Matrix.cpp
//  MatrixClass
//
//  Created by Li Xinrui on 10/1/14.
//  Copyright (c) 2014 Li Xinrui. All rights reserved.
//

#include "Matrix.h"
#include <cmath>
#include <cassert>
using std::cout;
using std::endl;
using std::cin;
// Algorithms
int LUDecomp(const Matrix& m, Matrix& l, Matrix& u){ // l, u blank Matrices
    int n = m.GetNumRows();
    if(m.GetNumRows() != m.GetNumCols()){
        std::cout << "Not Square" << std::endl;
        return 0;
    }
    int i, j, k; // Copy
    for(int i = 0; i != n; ++i)
        for(int j = 0; j != n; ++j){
            u(i,j) = m(i,j);
        }
    
    for(i = 0; i != n; ++i) l(i,i) = 1;
    
    for(k = 0; k != n - 1; ++k){ // L
        for(i = k + 1; i != n; ++i){
            l(i,k) = u(i,k) / u(k,k);
        }
        for(i = k + 1; i != n; ++i){ // U
            u(i,k) = 0;
            for(j = k + 1; j != n; ++j){
                u(i,j) = u(i,j) - l(i,k) * u(k,j);
            }
        }
    }
    return 0;
}

int LUPivotDecomp(const Matrix& m, Matrix& p, Matrix& l, Matrix& u){ // Mathematical Subscripts
    
    int n = m.GetNumRows();
    int i, j, k;
    int *per = new int[n - 1]; // at k, exchange GetNumRows() k with GetNumRows() per[k - 1]
    TYPE max = 0;
    TYPE temp = 0;
    //Test
    if(m.GetNumRows() != m.GetNumCols()){
        std::cout << "Not Square" << std::endl;
        return 0;
    }
    //Copy M into U
    for(i = 0; i != n; ++i)
        for(j = 0; j != n; ++j){
            u(i,j) = m(i,j);
        }
    //Initiate L
    for(i = 0; i != n; ++i) l(i,i) = 1;
    
    //Computation
    for(k = 1; k <= n; k++){
        // Find Permutations
        max = std::abs(u(k-1,k-1));
        per[k-1] = k;
        for(i = k; i <= n; i++){
            if(std::abs(u(i-1,k-1)) > max){
                max = std::abs(u(i-1,k-1));
                per[k-1] = i;
            }
        }
        u.rowSwap(k,per[k - 1]);
        // Compute L
        for(i = k + 1; i <= n; i++){
            l(i-1,k-1) = u(i-1,k-1) / u(k-1,k-1);
        }
        // Update U
        for(i = k + 1; i <= n; i++){
            u(i - 1,k - 1) = 0;
            for(j = k + 1; j <= n; ++j){
                u(i - 1,j - 1) = u(i - 1,j - 1) - l(i - 1,k - 1) * u(k - 1,j - 1);
            }
        }
    }
    // Compute P
    for(i = 1; i <= n; i++) p(i-1,i-1) = 1;
    for(k = 1; k <= n-1; k++) p.rowSwap(k, per[k-1]);
    
    // Compute L
    for(k = 1; k != n - 1; k++){
        if(per[k] != k + 1 ){
            for(j = 1; j != k + 1; j++){  // Exchange GetNumRows() k+1 and GetNumRows() per[k] in the Lower Left Matrix
                temp = l(k,j-1);
                l(k,j-1) = l( (per[k] - 1) ,j-1);
                l( (per[k] - 1) ,j-1) = temp;
            }
        }
    }
    
    return 0;
}

int CholeskyDecomp(const Matrix& m, Matrix& l){
    int i, p, k;
    int n = m.GetNumRows();
    TYPE temp;
    
    l(0,0) = sqrt(m(0,0));
    for(i = 1; i <= n; i++){
        l(i - 1,0) = m(i - 1,0) / l(0,0);
    }
    for(k = 2; k <= n; k++){
        temp = 0;
        for(p = 1; p <= k - 1; p++){
            temp += l(k - 1,p - 1) * l(k - 1,p - 1);
        }
        if(m(k - 1,k - 1) - temp < 0){ // Negative Diagonal
            cout << "Cholesky: Diagonal cannot Sqrt at " << k << endl;
            return 1;
        }
        else{
            l(k - 1,k - 1) = sqrt(m(k - 1,k - 1) - temp);
        }
        for(i = k + 1; i <= n; i++){
            temp = 0;
            for(p = 1; p <= k - 1; p++){
                temp += l(i - 1,p - 1) * l(k - 1,p - 1);
            }
            l(i - 1,k - 1) = (m(i - 1,k - 1) - temp) / l(k-1,k-1);
        }
    }
    
    return 0;
}

int LDLDecomp(const Matrix& m, Matrix& l, Matrix& d){
    int i, j, k;
    int n = m.GetNumRows();
    TYPE temp;
    for(i = 1; i <= n; i++){
        l(i - 1,i - 1) = 1;
    }
    for(j = 1; j <= n; j++){
        TYPE* v = new TYPE[j - 1];
        temp = 0;
        for(k = 1; k <= j - 1; k++){
            v[k - 1] = d(k - 1,k - 1) * l(j - 1,k - 1);
            temp += l(j - 1,k - 1) * v[k - 1];
        }
        d(j - 1,j - 1) = m(j - 1,j - 1) - temp;
        for(i = j + 1; i <= n; i++){
            temp = 0;
            for(k = 1; k <= j - 1; k++){
                temp += l(i -1,k - 1) * v[k - 1];
                
            }
            l(i - 1,j - 1) = (m(i-1,j-1) - temp) / d(j-1,j-1);
        }
        delete v;
    }
    return 0;
}


int LUSolve(const Matrix& m, const Matrix& rhs, Matrix& res){ // NEED RECYCLE
    if(rhs.GetNumCols() == 1){
        int n = m.GetNumRows();
        Matrix l(n,n), u(n,n);
        LUDecomp(m, l, u);
        Matrix temp(n,1);
        lowerSolve(l, rhs, temp);
        upperSolve(u, temp, res);
    }
    else{
        int n = m.GetNumRows();
        Matrix l(n,n), u(n,n);
        LUDecomp(m, l, u);
        for(int i = 1; i <= rhs.GetNumCols(); i++){
            Matrix temp1(n,1);
            Matrix temp2(n,1);
            Matrix right(n,1);
            for(int j = 1; j <= rhs.GetNumRows(); j++){
                right(j - 1,0) = rhs(j - 1,i - 1);
            }
            lowerSolve(l, right, temp1);
            upperSolve(u, temp1, temp2);
            for(int j = 1; j <= rhs.GetNumRows(); j++){
                res(j - 1,i - 1) = temp2(j - 1,0);
            }
        }
    }
    return 0;
}

int LUPivotSolve(const Matrix& m, const Matrix& rhs, Matrix& res){ // NEED
    if(rhs.GetNumCols() == 1){
        int n = m.GetNumRows();
        Matrix p(n,n), l(n,n),u(n,n);
        LUPivotDecomp(m, p, l, u);
        Matrix newb(n,1);
        newb = p * rhs;
        Matrix temp(n,1);
        lowerSolve(l, newb, temp);
        upperSolve(u, temp, res);
    }
    else{
        int n = m.GetNumRows();
        Matrix p(n,n), l(n,n),u(n,n);
        LUPivotDecomp(m, p, l, u);
        for(int i = 1; i <= rhs.GetNumCols(); i++){
            Matrix right(n,1);
            for(int j = 1; j <= rhs.GetNumRows(); j++){
                right(j-1,0) = rhs(j - 1,i - 1);
            }
            right = p * right;
            Matrix temp1(n,1);
            Matrix temp2(n,1);
            lowerSolve(l, right, temp1);
            upperSolve(u, temp1, temp2);
            for(int j = 1; j <= rhs.GetNumRows(); j++){
                res(j - 1,i - 1) = temp2(j - 1,0);
            }
        }
    }
    return 0;
}

int CholeskySolve(const Matrix&m, const Matrix& rhs, Matrix& res){
    int n = m.GetNumRows();
    Matrix y(n,1);
    Matrix l(n,n);
    Matrix trans(n,n);
    if(CholeskyDecomp(m, l)){
        cout << "Cholesky Solve: Cannot Decomp" << endl;
        return 1;
    }
    lowerSolve(l, rhs, y);
    l.transpose(trans);
    upperSolve(trans, y, res);
    return 0;
}
int LDLSolve(const Matrix& m, const Matrix& rhs, Matrix& res){
    int n = m.GetNumRows();
    Matrix y(n,1);
    Matrix l(n,n), d(n,n);
    Matrix trans(n,n);
    LDLDecomp(m, l, d);
    l.transpose(trans);
    lowerSolve(l, rhs, y);
    upperSolve(d * trans, y, res);
    return 0;
}

int lowerSolve(const Matrix& l,const Matrix& rhs, Matrix& res) // Array or Vector? I love vector.. Without Check on Dim
{
    int n = l.GetNumRows();
    res(0,0) = rhs(0,0) / l(0,0);
    for(int i = 1; i != n; ++i){
        TYPE temp = 0;
        for(int j = 0; j <= i - 1; ++j){
            temp += l(i,j) * res(j,0);
        }
        res(i,0) = (rhs(i,0) - temp) / l(i,i);
    }
    return 0;
}

int upperSolve(const Matrix& u, const Matrix& rhs, Matrix& res){
    assert(u.GetNumRows() == u.GetNumCols());
    int n = u.GetNumRows();
    res(n - 1,0) = rhs(n - 1,0) / u(n - 1,n - 1);
    for(int i = n - 2; i >= 0; --i){
        TYPE temp = 0;
        for(int j = i + 1; j != n; ++j){
            temp += u(i,j) * res(j,0);
        }
        res(i,0) = (rhs(i,0) - temp) / u(i,i);
    }
    return 0;
}


              


// Constructors
Matrix::Matrix():
mNumRows(0), mNumCols(0), data(nullptr) {
}

Matrix::Matrix(int inputRow, int inputCol):
mNumRows(inputRow), mNumCols(inputCol), data(new TYPE*[mNumCols]) {
    int i,j;
    for(i = 0; i != mNumCols; ++i){
        data[i] = new TYPE[mNumRows];
    }
    for (i = 0; i != mNumCols; ++i) {
        for (j = 0; j != mNumRows; ++j) {
            data[i][j] = 0;
        }
    }
}

Matrix::Matrix(int dim):
mNumRows(dim), mNumCols(dim), data(new TYPE*[dim]) {
    int i,j;
    for(i = 0; i != dim; ++i){
        data[i] = new TYPE[dim];
    }
    for (i = 0; i != mNumCols; ++i) {
        for (j = 0; j != mNumRows; ++j) {
            data[i][j] = 0;
        }
    }
}

Matrix::Matrix(const Matrix& m):
mNumRows(m.mNumRows),mNumCols(m.mNumCols),data(new TYPE*[mNumCols]){
    for(int i = 0; i != mNumCols; ++i){
        data[i] = new TYPE[mNumRows];
    }
    for(int i = 0; i != mNumCols; ++i)
        for(int j = 0; j!= mNumRows; ++j){
            data[i][j] = m[i][j];
        }
}

// Destructor

// Operators

TYPE* Matrix::operator[](int index){
    return data[index];
}

const TYPE* Matrix::operator[](int index) const{
    return data[index];
}

TYPE& Matrix::operator()(int rowIndex, int colIndex){
//    assert(rowIndex < mNumRows && colIndex < mNumCols);
    return data[colIndex][rowIndex];
}

const TYPE& Matrix::operator()(int rowIndex, int colIndex) const{
    return data[colIndex][rowIndex];
}



Matrix& Matrix::operator+=(const Matrix& rhs){
    if(mNumCols != rhs.mNumCols || mNumRows != rhs.mNumRows){
        std::cout << "Unequal Dimensions" << std::endl;
    }
    for(size_t i = 0; i != rhs.mNumCols; ++i){
        for(size_t j = 0; j != rhs.mNumRows; ++j){
            data[i][j] += rhs.data[i][j];
        }
    }
    return *this;
}
Matrix& Matrix::operator-=(const Matrix& rhs){
    return *this += rhs * -1;
}

Matrix& Matrix::operator*=(const TYPE& rhs){
    for(int i = 0; i != mNumCols; i++)
        for(int j = 0; j != mNumRows; j++) {
            data[i][j] *= rhs;
        }
    return *this;
}

Matrix& Matrix::operator/=(const TYPE& rhs){
    return *this *= 1/rhs;
}

Matrix& Matrix::operator=(const TYPE* rhs){
    int k = 0;
    for(int i = 0; i != mNumCols; i++)
        for(int j = 0; j != mNumRows; j++) {
            data[i][j] = rhs[k];
            k++;
        }
    return *this;
}

Matrix& Matrix::operator=(const vector<TYPE>& rhs){
    int k = 0;
    for(int i = 0; i != mNumCols; i++)
        for(int j = 0; j != mNumRows; j++) {
            data[i][j] = rhs[k];
            k++;
        }
    return *this;
}
Matrix& Matrix::operator=(const Matrix& rhs){
    int i,j;
    if(data == nullptr){
        mNumRows = rhs.mNumRows;
        mNumCols = rhs.mNumCols;
        data = new TYPE*[rhs.mNumCols];
        for(int i = 0; i != rhs.mNumCols; ++i){
            data[i] = new TYPE[rhs.mNumRows];
        }
    }
    if(mNumCols != rhs.mNumCols || mNumRows != rhs.mNumRows){
        cout << "unequal dimensions cannot assign matrix" << endl;
    }
    else{
        for(i = 0; i != rhs.mNumCols; i++)
            for(j = 0; j != rhs.mNumRows; j++){
                data[i][j] = rhs[i][j];
            }
    }
    return *this;
}


Matrix& Matrix::colSwap(int i, int j){
    if(i != j){
        TYPE *p = data[i-1];
        data[i-1] = data[j-1];
        data[j-1] = p;
    }
    return *this;
}

Matrix& Matrix::rowSwap(int i, int j) { // Mathematical Subscripts
    TYPE temp;
    for(int k = 0; k != mNumCols; k++){
        temp = data[k][i - 1];
        data[k][i - 1] = data[k][j - 1];
        data[k][j - 1] = temp;
    }
    return *this;
}
bool Matrix::QSquare(){
    assert(mNumCols == mNumRows);
    if(mNumRows == mNumCols){
        return 1;
    }
    else return 0;
}

double Matrix::Norm2Vec(){
    int i,j;
    double sum = 0;
    for (i = 0; i != mNumCols; i++) {
        for (j = 0; j != mNumRows; j++) {
            sum += data[i][j] * data[i][j];
        }
    }
    return sqrt(sum);
}

double mVector::NormInf(){
    double max = 0;
    for (int i = 0; i != this->mNumRows; i++) {
        if (std::abs((double) this->data[0][i]) > max ) {
            max = std::abs((double) this->data[0][i]);
        }
    }
    return max;
};

double mVector::NormInf(int& position){
    double max = 0;
    for (int i = 0; i != this->mNumRows; i++) {
        if (std::abs((double) this->data[0][i]) > max ) {
            max = std::abs((double) this->data[0][i]);
            position = i;
        }
    }
    return max;
};

double InnerProduct(mVector &v1, mVector &v2){
    double product = 0;
    assert(v1.dim() == v2.dim());
    for (int i = 0; i != v1.dim(); i++) {
        product += v1[i] * v2[i];
    }
    return product;
}

double mVector::Norm_2(){
    double sum = 0;
    for (int i = 0; i != this -> mNumRows; i++) {
        sum += this -> data[0][i] * data[0][i];
    }
    return sqrt(sum);
}

double mVector::Norm_1(){
    double sum = 0;
    for (int i = 0; i != this -> mNumRows; i++) {
        sum += std::abs(data[0][i]);
    }
    return sum;
}

Matrix& Matrix::transpose(){
    int i, j;
    TYPE temp;
    // This changes the Matrix itself
    if(mNumRows != mNumCols){
        Matrix copy(*this);
        for (i = 0; i != mNumCols; i++) {
            delete [] data[i];
        }
        delete [] data;
        mNumCols = copy.mNumRows;
        mNumRows = copy.mNumCols;
        data = new TYPE*[mNumCols];
        for(i = 0; i != mNumCols; i++){
            data[i] = new TYPE[mNumRows];
        }
        for(i = 1; i <= mNumCols; i++)
            for(j = 1; j <= mNumRows; j++){
                data[i - 1][j - 1] = copy[j - 1][i - 1];
            }

        return *this;
    }
    for(i = 1; i <= mNumCols; i++)
        for(j = 1; j < i; j++){
            if(data[i - 1][j - 1] != data[j - 1][i - 1]){
                temp = data[i - 1][j - 1];
                data[i - 1][j - 1] = data[j - 1][i - 1];
                data[j - 1][i - 1] = temp;
            }
        }
    return *this;
}

int Matrix::transpose(Matrix& res){
    int i, j;
    for(i = 1; i <= mNumCols; i++)
        for(j = 1; j <= mNumRows; j++){
            res[j - 1][i - 1] = data[i - 1][j - 1];
        }
    return 0;
}


// Public
Matrix operator+(const Matrix& m1, const Matrix& m2){
    Matrix sum(m1.GetNumRows(), m1.GetNumCols());
    sum = m1;
    sum += m2;
    return sum;
}

Matrix operator-(const Matrix& m1, const Matrix& m2){
    Matrix sum(m1.GetNumRows(), m1.GetNumCols());
    sum = m1;
    sum -= m2;
    return sum;
}
Matrix operator*(const Matrix& m1, const Matrix& m2){
    if(m1.GetNumCols() != m2.GetNumRows()) {
        std::cout << "Illegal Multiplication" << std::endl;
    }
    Matrix res(m1.GetNumRows(), m2.GetNumCols());
    for(int i = 0; i != m1.GetNumRows(); ++i)
        for(int j = 0; j != m2.GetNumCols(); ++j)
            for(int k = 0; k != m1.GetNumCols(); ++k){
                res(i,j) += (m1(i,k) * m2(k,j));
            }
    return res;
}
Matrix operator*(const Matrix& m1, const TYPE& num){
    Matrix res(m1.GetNumRows(),m1.GetNumCols());
    res = m1;
    res *= num;
    return res;
}
Matrix operator*(const TYPE& num, const Matrix& m1){
    return m1 * num;
}

Matrix operator/(const Matrix& m1, const TYPE& num){
    Matrix res(m1.GetNumRows(),m1.GetNumCols());
    res = m1;
    res /= num;
    return res;
}
mVector operator+(const mVector& m1, const mVector& m2){
    mVector res(m1.dim());
    for (int i = 0; i != m1.dim(); i++) {
        res[i] = m1[i] + m2[i];
    }
    return res;
}
mVector operator-(const mVector& m1, const mVector& m2){
    mVector res(m1.dim());
    for (int i = 0; i != m1.dim(); i++) {
        res[i] = m1[i] - m2[i];
    }
    return res;
}
mVector operator*(const Matrix& m1, const mVector& m2){
    assert(m1.GetNumCols() == m2.GetNumRows());
    mVector res(m1.GetNumRows());
    for(int i = 0; i != m1.GetNumRows(); ++i)
            for(int k = 0; k != m1.GetNumCols(); ++k){
                res[i] += (m1(i,k) * m2[k]);
            }
    return res;
}
mVector operator*(const mVector& m1, const TYPE& num){
    mVector res(m1.GetNumRows());
    for (int i = 0; i != res.dim(); i++) {
        res[i] = m1[i] * num;
    }
    return res;
}
mVector operator*(const TYPE& num, const mVector& m1){
    return m1 * num;
}
mVector operator/(const mVector& m1, const TYPE& num){
    return m1 * 1/ num;
}

// Output & Style
std::ostream& operator<<(std::ostream& out, const Matrix& m){
    int i, j;
    if(m.data == nullptr) cout << "Empty Matrix";
    for (i = 0; i != m.mNumRows; i++){
        for(j = 0; j!= m.mNumCols; j++){
            out << m[j][i] << " ";
        }
        out << std::endl;
    }
    return out;
}

std::ofstream& operator<<(std::ofstream& out, Matrix& m){
    int i, j;
    if(m.data == nullptr) cout << "Empty Matrix";
    for (i = 0; i != m.mNumRows; i++){
        for(j = 0; j!= m.mNumCols; j++){
            out << m[j][i] << " " ;
        }
        out << std::endl;
    }
    return out;
}

std::ostream& latex(std::ostream& out, const Matrix& m, int line){
    int i, j;
    int t = 0;
    if(m.data == nullptr) cout << "Empty Matrix";
    for (i = 0; i != m.mNumRows; i++){
        for(j = 0; j!= m.mNumCols; j++){
            if(i == m.mNumRows - 1 && j == m.mNumCols - 1){
                out << m[j][i] << " " << "\\\\";
            }
            else{
                t++;
                if(t == line){
                    out << m[j][i] << " ";
                    out << "\\\\";
                    out << endl;
                    out << "\\hline";
                    out << endl;
                    t = 0;
                }
                else{
                    out << m[j][i] << " " << "&" << " ";
                }
            }
        }
    }
    return out;
}


