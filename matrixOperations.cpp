
#include<sstream>
#include "lamdafunctions.h"
#include "bits/stdc++.h"

using namespace std;

template<typename T>
struct Normalize {
    T operator()(const T &Left, const T &Right) const {
        return (Left / pow(Right, 0.5));
    }
};

vector<double> getMean(vector<vector<double>> images, int sample_size) {

    vector<double> mean;
    mean.resize(images[0].size());

    for (vector<double> samples :  images) {
        std::transform(mean.begin(), mean.end(), samples.begin(), mean.begin(), std::plus<double>());
    }

    std::transform(mean.begin(), mean.end(), mean.begin(),
                   [sample_size](double c) { return (double) c / (double) sample_size; });
    return mean;
};

vector<vector<double> > getCovairiance(vector<vector<double>> A, int sample_size, int image_size) {
    vector<vector<double> > AT(A[0].size(), vector<double>());

    for (int i = 0; i < A.size(); i++) {
        for (int j = 0; j < A[i].size(); j++) {
            AT[j].push_back(A[i][j]);
        }
    }
    //displayMatrix(AT);
    pair<int, int> dimen1 = {sample_size, image_size};
    pair<int, int> dimen2 = {image_size, sample_size};
    vector<vector<double> > mul;
    mul.resize(dimen1.first);
    for (int i = 0; i < dimen1.first; i++) {
        mul[i].resize(dimen2.second);
        for (int j = 0; j < dimen2.second; j++) {
            mul[i][j] = 0;
            for (int k = 0; k < dimen2.first; k++) {
                mul[i][j] += A[i][k]
                             * (AT[k][j]);
            }
            mul[i][j] = mul[i][j] / sample_size;
        }
    }
    return mul;
}


void normalize(vector<vector<double>> &A) {
    vector<double> length;
    length.resize(A[0].size());

    for (vector<double> samples :  A) {
        std::transform(length.begin(), length.end(), samples.begin(), length.begin(), square<double>());
    }
    //displayMatrix(length);
    std::for_each(std::begin(A), std::end(A),
                  [length](vector<double> (&row)) {
                      std::transform(std::begin(row), std::end(row), std::begin(length), std::begin(row),
                                     Normalize<double>());

                  });
}

void vitualization(vector<vector<double>> &A) {

    std::for_each(std::begin(A), std::end(A),
                  [](vector<double> (&row)) {
                      pair<vector<double>::iterator, vector<double>::iterator> min_max = minmax_element(row.begin(),
                                                                                                        row.end());
                      double min = row[min_max.first - row.begin()];
                      double max = row[min_max.second - row.begin()];
                      std::transform(std::begin(row), std::end(row), std::begin(row),
                                     [min, max](double x) { return 255 * (x - min) / (max - min); });
                  });
}


vector<vector<double>>
getrealEigenVectors(vector<vector<double>> A, vector<vector<double>> V, int sample_size, int image_size) {
    vector<vector<double> > AT(A[0].size(), vector<double>());

    for (int i = 0; i < A.size(); i++) {
        for (int j = 0; j < A[i].size(); j++) {
            AT[j].push_back(A[i][j]);
        }
    }
    pair<int, int> dimen1 = {image_size, sample_size};
    pair<int, int> dimen2 = {sample_size, sample_size};
    vector<vector<double> > mul;
    mul.resize(dimen1.first);
    for (int i = 0; i < dimen1.first; i++) {
        mul[i].resize(dimen2.second);
        for (int j = 0; j < dimen2.second; j++) {
            mul[i][j] = 0;
            for (int k = 0; k < dimen2.first; k++) {
                mul[i][j] += AT[i][k]
                             * (V[k + 1][j + 1]);
            }
        }
    }
    //displayMatrix(mul);
    normalize(mul);
    return mul;
}

vector<vector<double> >
getEigenspace(vector<vector<double>> U, vector<vector<double>> X, int sample_size, int image_size, int eigensize) {

    pair<int, int> dimen1 = {sample_size, image_size};
    pair<int, int> dimen2 = {image_size, eigensize};
    vector<vector<double> > mul;
    mul.resize(dimen1.first);
    for (int i = 0; i < dimen1.first; i++) {
        mul[i].resize(dimen2.second);
        for (int j = 0; j < dimen2.second; j++) {
            mul[i][j] = 0;
            for (int k = 0; k < dimen2.first; k++) {
                mul[i][j] += X[i][k]
                             * (U[k][j]);
            }
        }
    }
    return mul;
}

vector<vector<double>>
transpose(vector<vector<double>> A) {
    vector<vector<double> > AT(A[0].size(), vector<double>());

    for (int i = 0; i < A.size(); i++) {
        for (int j = 0; j < A[i].size(); j++) {
            AT[j].push_back(A[i][j]);
        }
    }
    return AT;
}

vector<vector<double> >
matrixMultiply(vector<vector<double>> U, vector<vector<double>> X, pair<int, int> dimen1,pair<int, int> dimen2) {

    vector<vector<double> > mat;
    mat.resize(dimen1.first);
    for (int i = 0; i < dimen1.first; i++) {
        mat[i].resize(dimen2.second);
        for (int j = 0; j < dimen2.second; j++) {
            mat[i][j] = 0;
            for (int k = 0; k < dimen2.first; k++) {
                mat[i][j] += U[i][k]
                             * (X[k][j]);
            }
        }
    }
    return mat;
}