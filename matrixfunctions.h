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


#endif //DIMENRED_MATRIXFUNCTIONS_H
