#ifndef PATERNREC_MATRICOP_H
#define PATERNREC_MATRICOP_H

#endif //PATERNREC_MATRICOP_H

#define rows 2
#define cols 2

#include "bits/stdc++.h"

using namespace std;
// Matrix class is implement to simulate Matrix and its operations
// + , - ,* ,== ,<  operator are overloaded for direct use in main class
// transpose, determinant and inverse function: trans(), deter(),inv()
// are abstracted in Matrix class

// other functions
//getSigmaSqr() to return value of sigma in case covariance matrix is diagonal and have same variance

//getDouble() - it returns the value at 0,0 . it is used in case matrix is 1X1

class Matrix {

    double arr[rows][cols];
    pair<int, int> dimen;

public:
    Matrix() {}

    /*~Matrix() {
        delete arr;
    }*/


    Matrix(double **ar, pair<int, int> dimen) {
        this->dimen = dimen;
        for (int i = 0; i < dimen.first; i++) {
            for (int j = 0; j < dimen.second; j++) {
                arr[i][j] = ar[i][j];
            }
        }
    }
    Matrix(vector<vector<double> > &A, pair<int, int> dimen) {
        this->dimen = dimen;
        for (int i = 0; i < dimen.first; i++) {
            for (int j = 0; j < dimen.second; j++) {
                arr[i][j] = A[i][j];
            }
        }
    }
    void input(vector<vector<double> > &A, pair<int, int> dimen);
    void display();

    Matrix operator+(Matrix x);
    Matrix operator-(Matrix x);
    Matrix operator*(Matrix x);
    Matrix operator*(double x);
    bool operator==(Matrix x);
    bool operator<(const Matrix &X) const;
    Matrix trans();
    Matrix inv();
    double deter();


    double getSigmaSqr() {
        double temp = arr[0][0];
        for (int i = 0; i < dimen.first; i++) {
            for (int j = 0; j < dimen.second; j++) {
                if (!((arr[i][j] == temp && i == j) || (arr[i][j] == 0 && i != j)))
                    return 0;
            }
        }
        return temp;
    }

    double getDouble() { return arr[0][0]; }
    vector<double> getVector() { return {arr[0][0], arr[1][0]}; }

    double deterR(double a[rows][cols], int k) {
        double s = 1, det = 0;
        double b[rows][cols];
        int i, j, m, n, c;
        if (k == 1) {
            return (a[0][0]);
        } else {
            det = 0;
            for (c = 0; c < k; c++) {
                m = 0;
                n = 0;
                for (i = 0; i < k; i++) {
                    for (j = 0; j < k; j++) {
                        b[i][j] = 0;
                        if (i != 0 && j != c) {
                            b[m][n] = a[i][j];
                            if (n < (k - 2))
                                n++;
                            else {
                                n = 0;
                                m++;
                            }
                        }
                    }
                }
                det = det + s * (a[0][c] * deterR(b, k - 1));
                s = -1 * s;
            }
        }
        return (det);
    }

    double **cofactor(double num[rows][cols], int f) {
        double b[rows][cols];
        double **fac = new double *[f];
        int m, n;
        for (int q = 0; q < f; q++) {
            fac[q] = new double[f];

            for (int p = 0; p < f; p++) {
                m = 0;
                n = 0;
                for (int i = 0; i < f; i++) {
                    for (int j = 0; j < f; j++) {
                        if (i != q && j != p) {
                            b[m][n] = num[i][j];
                            if (n < (f - 2))
                                n++;
                            else {
                                n = 0;
                                m++;
                            }
                        }
                    }
                }
                fac[q][p] = pow(-1, q + p) * deterR(b, f - 1);
            }
        }
        return fac;
    }
};

void Matrix::input(vector<vector<double> > &A, pair<int, int> dimen) {
    this->dimen = dimen;
    for (int i = 0; i < dimen.first; i++) {
        for (int j = 0; j < dimen.second; j++) {
            arr[i][j] = A[i][j];
        }
    }
}

void Matrix::display() {
    cout << "dimentions " << this->dimen.first << "x" << this->dimen.second << '\n';
    for (int i = 0; i < dimen.first; i++) {
        for (int j = 0; j < dimen.second; j++) {
            cout << arr[i][j] << ' ';
        }
        cout << endl;
    }
}

Matrix Matrix::operator+(Matrix x) {

    vector<vector<double> > mat;
    mat.resize(x.dimen.first);
    for (int i = 0; i < dimen.first; i++) {
        mat[i].resize(x.dimen.second);
        for (int j = 0; j < x.dimen.second; j++) {
            mat[i][j] = arr[i][j]+ x.arr[i][j];
        }
    }
    return Matrix(mat, {dimen.first, x.dimen.second});
}

Matrix Matrix::operator-(Matrix x) {

    vector<vector<double> > mat;
    mat.resize(x.dimen.first);
    for (int i = 0; i < dimen.first; i++) {
        mat[i].resize(x.dimen.second);
        for (int j = 0; j < x.dimen.second; j++) {
            mat[i][j] = arr[i][j]  - x.arr[i][j];
        }
    }
    return Matrix(mat, {dimen.first, x.dimen.second});
}

Matrix Matrix::operator*(Matrix x) {

    vector<vector<double> > mat;
    mat.resize(dimen.first);
    for (int i = 0; i < dimen.first; i++) {
        mat[i].resize(x.dimen.second);
        for (int j = 0; j < x.dimen.second; j++) {
            mat[i][j] = 0;
            for (int k = 0; k < x.dimen.first; k++) {
                mat[i][j] += arr[i][k]
                             * (x.arr[k][j]);
            }
        }
    }
    return Matrix(mat, {dimen.first, x.dimen.second});
}

Matrix Matrix::operator*(double x) {

    vector<vector<double> > mat;
    mat.resize(dimen.first);
    for (int i = 0; i < dimen.first; i++) {
        mat[i].resize(dimen.second);
        for (int j = 0; j < dimen.second; j++) {
            mat[i][j] = x * arr[i][j];
        }
    }
    return Matrix(mat, dimen);
}


bool Matrix::operator<(const Matrix &X) const {
    bool less = false;
    if (arr[0][0] < X.arr[0][0]) { less = true; }
    else if (arr[0][0] == X.arr[0][0] and arr[1][0] < X.arr[1][0]) { less = true; }
    return less;
}

bool Matrix::operator==(Matrix x) {
    bool check = true;
    if (dimen != x.dimen) {
        check = false;
    } else {
        for (int i = 0; i < dimen.first; i++) {
            for (int j = 0; j < x.dimen.second; j++) {
                if (arr[i][j] != x.arr[i][j]) {
                    check = false;
                    break;
                }
            }
        }
    }
    return check;
}


Matrix Matrix::trans() {
    vector<vector<double> > b;
    b.resize(dimen.second);
    for (int i = 0; i < dimen.second; i++) {

        b[i].resize(dimen.first);
        for (int j = 0; j < dimen.first; j++) {
            b[i][j] = arr[j][i];
        }
    }
    return Matrix(b, {dimen.second, dimen.first});
}

double Matrix::deter() {return deterR(arr, dimen.second);}

Matrix Matrix::inv() {

    double d;
    int i, j;
    int k = dimen.second;
    double **fac;
    double **inverse = new double *[k];
    d = deterR(arr, k);
    if (d == 0)
        printf("\nInverse of Entered Matrix is not possible\n");
    else {
        fac = cofactor(arr, k);
        for (int i = 0; i < k; i++) {
            inverse[i] = new double[k];
            for (int j = 0; j < k; j++) {
                inverse[i][j] = fac[j][i] / d;
            }
        }
    }
    return Matrix(inverse, dimen);
}
