#include <iostream>
#include "image.h"
#include<string.h>
#include<sstream>
#include<dirent.h>
#include "bits/stdc++.h"

using namespace std;

void readImageHeader(string, int &, int &, int &, bool &);
void readImage(string, ImageType &);
void writeImage(string, ImageType &);

void writeImages(vector<double> meanfaces) {
    int N = 20, M = 16;
    ImageType image(N, M, 255);
    for (int i = 0; i < N; i++)
        for (int j = 0; j < M; j++) {
            image.setPixelVal(i, j, meanfaces[i * M + j]);
        }
    string filename = string("meanface/") + string("meanface.pgm");
    writeImage(filename, image);

}

vector<double> getReferenceImageMatrix(string sample_file) {
    int sval;
    int SM, SN, SQ;
    bool Stype;
    readImageHeader(sample_file, SN, SM, SQ, Stype);
    // allocate memory for the image array
    ImageType Simage(SN, SM, SQ);
    // read image
    readImage(sample_file, Simage);
    /*SM = 2;
    SN = 2;*/
    // cout << SN<<","<< SM<<endl;
    vector<double> ref;
    for (int i = 0; i < SN; i++)
        for (int j = 0; j < SM; j++) {
            Simage.getPixelVal(i, j, sval);
            ref.push_back(sval);
        }
    return ref;
}
void writeImages(vector<vector<double>> eigenfaces) {
    int N = 20, M = 16;
    int fcno = 0;
    for (vector<double> eg : eigenfaces) {
        fcno++;
        ImageType image(N, M, 255);
        for (int i = 0; i < N; i++) {

            for (int j = 0; j < M; j++) {
                image.setPixelVal(i, j, eg[i * M + j]);
            }
        }
        string filename = string("eigenfaces/") + to_string(fcno) + string("eigenfacesC.pgm");
        writeImage(filename, image);
    }
}
vector<string> listFile(string direc) {
    ifstream inn;
    DIR *pDIR;
    vector<string> file_list;
    struct dirent *entry;
    string path = "./" + direc;
    char *c = strcpy(new char[path.length() + 1], path.c_str());

    if (pDIR = opendir(c)) {
        while (entry = readdir(pDIR)) {
            if (strcmp(entry->d_name, ".") != 0 && strcmp(entry->d_name, "..") != 0) {
                string str(entry->d_name);
                file_list.emplace_back(direc + "/" + str);
            }

        }
        closedir(pDIR);
    }

    return file_list;

}


vector<vector<double>> getX(string direc) {
    vector<vector<double>> images;
    vector<string> file_list = listFile(direc);
    int k = 399;
    for (string fn :  file_list) {
        vector<double> ref_image = getReferenceImageMatrix(fn);
        images.push_back(ref_image);
        k--;
        if (k < 0)
            break;
    }
    return images;
}

void write_sol_file(vector<vector<double>> matrix) {
    ofstream myfile("advertise_pb_out.txt");
    if (myfile.is_open()) {
        for (vector<double> vector :  matrix) {
            for (double j :  vector) {
                cout << j << ",";
            }
            cout << " \n";
        }
        myfile.close();
    }
}

void write_sol_file(vector<double> vector) {
    ofstream myfile("mean.csv");
    if (myfile.is_open()) {
        for (double j :  vector) {
            cout << j << ",";
        }
        cout << " \n";
    }
}

void displayMatrix(vector<vector<double>> matrix) {
    cout << endl << "size: " << matrix.size() << endl;

    for (vector<double> vector :  matrix) {


        for (double j :  vector) {

            cout << j << " ";

        }
        cout << " \n";
    }
}

void displayMatrix(vector<double> vector) {
    cout << endl << "size: " << vector.size() << endl;
    for (double j :  vector) {

        cout << j << " ";

    }
    cout << " \n";
}

int computeKVal(int sampleSize, vector<double> lambda, float threshold) {
    double eigenSum = 0;
    double totalEigenSum = 0;
    for (int i = 1; i < sampleSize; i++) {
        totalEigenSum += lambda[i];
    }

    for (int i = 1; i < sampleSize; i++) {
        eigenSum += lambda[i];
        //cout << i << " " << lambda[i] << " " << eigenSum << " " << (eigenSum/totalEigenSum) << endl;
        if ((eigenSum/totalEigenSum) > threshold) {
            return i;
        }
    }
}

void createTrainingDir() {
    string str2 = "mkdir training";
    const char *command2 = str2.c_str();
    system(command2);
}

bool exists_test(const string& name) {
    ifstream f(name.c_str());
    return f.good();
}

void writeMean(vector<double> mean) {
    string filename = "training/mean.csv";
    ofstream outFile;
    outFile.open(filename, ios::trunc);

    for (int i = 0; i < mean.size(); i++) {
        outFile << mean[i] << endl;
    }
}

void writeOmega(vector<vector<double>> omega) {
    string filename = "training/omega.csv";
    ofstream outFile;
    outFile.open(filename, ios::trunc);

    for (int i = 0; i < omega.size(); i++) {
        for(int j = 0; j < omega[i].size(); j++) {
            outFile << omega[i][j] << ",";
        }
        outFile << endl;
    }
}

void writeLambda(vector<double> lambda) {
    string filename = "training/lambda.csv";
    ofstream outFile;
    outFile.open(filename, ios::trunc);

    for (int i = 0; i < lambda.size(); i++) {
        outFile << lambda[i] << endl;
    }
}

void writeEigenface(vector<vector<double>> eigenface) {
    string filename = "training/eigenface.csv";
    ofstream outFile;
    outFile.open(filename, ios::trunc);

    for (int i = 0; i < eigenface.size(); i++) {
        for(int j = 0; j < eigenface[i].size(); j++) {
            outFile << eigenface[i][j] << ",";
        }
        outFile << endl;
    }
}

vector<double> readMean(string filename) {
    vector<string> convert;
    ifstream inFile(filename);
    string data;
    if (inFile.is_open()) {
        string line;
        while (!inFile.eof()) {
            getline(inFile, data, '\n');
            convert.push_back(data);
        }
        inFile.close();
    }
    vector<double> mean;
    for (int i = 0; i < convert.size()-1; i++) {
        mean.push_back(stod(convert[i]));
    }
    return mean;
    //convertToDouble(convert);
}

vector<vector<double>> readOmega(string filename) {
    vector<vector<double>> omega;

    ifstream inFile(filename);
    string data;
    int count = 0;
    if (inFile.is_open()) {
        string line;
        while (!inFile.eof()) {
            vector<string> convert;

            while (inFile.peek() != '\n') {
                getline(inFile, data, ',');
                convert.push_back(data);
            }
            getline(inFile, data, '\n');

            vector<double> singleLine;
            for (int i = 0; i < convert.size()-1; i++) {
                singleLine.push_back(stod(convert[i]));
            }
            omega.push_back(singleLine);
        }
        inFile.close();
    }

    return omega;
}

vector<double> readLambda(string filename) {
    vector<string> convert;
    ifstream inFile(filename);
    string data;
    if (inFile.is_open()) {
        string line;
        while (!inFile.eof()) {
            getline(inFile, data, '\n');
            convert.push_back(data);
        }
        inFile.close();
    }
    vector<double> lambda;
    for (int i = 0; i < convert.size()-1; i++) {
        lambda.push_back(stod(convert[i]));
    }
    return lambda;
}

void readEigenface(vector<vector<double>> eigenface) {

}

void convertToDouble(vector<string> convert) {

}