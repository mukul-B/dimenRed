//
// Created by dell on 4/23/2021.
//

#ifndef DIMENRED_MATRIXOPERATIONS_H
#define DIMENRED_MATRIXOPERATIONS_H

vector<double> getMean(vector<vector<double>> images, int sample_size);

vector<vector<double>> getCovairiance(vector<vector<double>> A, int sample_size, int image_size);

vector<vector<double>>
getrealEigenVectors(vector<vector<double>> A, vector<vector<double>> V, int sample_size, int image_size);

void normalize(vector<vector<double>> &A);

vector<vector<double>>
getEigenspace(vector<vector<double>> U, vector<vector<double>> X, int sample_size, int image_size, int eigensize);

void vitualization(vector<vector<double>> &A);

#endif //DIMENRED_MATRIXOPERATIONS_H
