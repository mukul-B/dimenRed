//
// Created by dell on 4/23/2021.
//

#ifndef DIMENRED_MATRIXFUNCTIONS_H
#define DIMENRED_MATRIXFUNCTIONS_H

vector<vector<double>> getX(string direc);
vector<double> getReferenceImageMatrix(string sample_file);
vector<string> listFile(string);

void writeImages(vector<vector<double>> eigenfaces, string resolution);
void writeImages(vector<double> meanfaces, string resolution);
void displayMatrix(vector<double> vector);
void displayMatrix(vector<vector<double>> matrix);
void writeResults(vector<vector<double>> vector,string file_name);
void writeResults(vector<double> eigenvalue,string file_name);


vector<vector<double>> readResults(string file_name);

#endif //DIMENRED_MATRIXFUNCTIONS_H
