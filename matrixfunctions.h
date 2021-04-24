//
// Created by dell on 4/23/2021.
//

#ifndef DIMENRED_MATRIXFUNCTIONS_H
#define DIMENRED_MATRIXFUNCTIONS_H

vector<vector<double>> getX(string direc);
vector<double> getReferenceImageMatrix(string sample_file);
vector<string> listFile(string);

void writeImages(vector<vector<double>> eigenfaces);
void writeImages(vector<double> meanfaces);
void displayMatrix(vector<double> vector);
void displayMatrix(vector<vector<double>> matrix);

int computeKVal(int sampleSize, vector<double> lambda, float threshold);

void createTrainingDir();
bool exists_test(const string& name);

void writeMean(vector<double> mean);
void writeOmega(vector<vector<double>> omega);
void writeLambda(vector<double> lambda);
void writeEigenface(vector<vector<double>> eigenface);

vector<double> readMean(string filename);
vector<vector<double>> readOmega(string filename);
vector<double> readLambda(string filename);
void readEigenface(vector<vector<double>> eigenface);

void convertToDouble(vector<string> convert);

#endif //DIMENRED_MATRIXFUNCTIONS_H
