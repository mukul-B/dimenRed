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
       /* k--;
        if (k < 0)
            break;*/
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



vector<vector<double>> readResults(string file_name) {
    fstream fin;
    fin.open(file_name, ios::in);
    vector<double> row;
    vector<vector<double>> result;
    string line,temp,word;

    while (fin >> temp) {
        row.clear();
        getline(fin, line);
        stringstream s(temp);
        while (getline(s, word, ',')) {
            row.push_back(stod(word));
        }
        result.push_back(row);
    }
    return result;
}

void writeResults(vector<vector<double>> eigenvectors,string file_name) {
    ofstream myfile(file_name);
    if (myfile.is_open()) {
        for (vector<double> eg : eigenvectors) {
            for (double e : eg) {
                myfile <<e<<",";
            }
            myfile<<endl;
        }
        myfile.close();
    }
}
void writeResults(vector<double> eigenvalue,string file_name) {
    ofstream myfile(file_name);
    if (myfile.is_open()) {

        for (double e : eigenvalue) {
            myfile <<e<<",";
        }

        myfile.close();
    }
}


