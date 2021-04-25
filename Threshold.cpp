#include <iostream>
#include <math.h>
#include<sstream>
#include "jacobi.h"
#include "lamdafunctions.h"
#include "bits/stdc++.h"
#include "matrixOperations.h"
#include "matrixfunctions.h"

using namespace std;
using namespace std::chrono;

int training();

int testing(vector<double> mean, vector<vector<double> > eigespace, vector<double> w,
            vector<vector<double> > U,int rank);

int testing(string train_dir, string test_dir);
int reconstructionCheck(string train_dir, string test_dir);

int main(int argc, char *argv[]) {
    string training_folder = "fa_H";
    string testing_folder = "fb_H";
    testing(training_folder,testing_folder);
    //training();
    //reconstructionCheck(training_folder, testing_folder);
    return 0;
}


int testing(vector<double> mean, vector<vector<double> > eigespace, vector<double> w,
            vector<vector<double> > U,string train_dir, string test_dir ,int rank= 50) {

    vector<string> train_file_list = listFile(train_dir);
    vector<string> test_file_list = listFile(test_dir);
    std::transform(test_file_list.begin(), test_file_list.end(), test_file_list.begin(), findimageId<string>());
    std::transform(train_file_list.begin(), train_file_list.end(), train_file_list.begin(), findimageId<string>());

    vector<vector<double>> testImages = getX(test_dir);

    int test_size = testImages.size();
    int image_size = mean.size();
    int galery_size = eigespace.size();
    int iop = 0;
    int correctCount = 0;
    int kValue=0;
    double eigenSum = 0;
    double totalEigenSum = 0;
    for (int i = 0; i < image_size; i++) {
        totalEigenSum += w[i];
    }
    for (int i = 0; i < image_size; i++) {
        eigenSum += w[i];
        if ((eigenSum / totalEigenSum) > 0.8) {
            kValue = i - 1;
            break;
        }
    }
    // kValue=1000;
    /*cout << kValue << endl;
    cout << galery_size << endl;
    cout << image_size << endl;
    cout << test_size << endl;*/
    //displayMatrix(w);
    for (vector<double> X : testImages) {
//        vector<double> X = getReferenceImageMatrix("fb_l/00019_940422_fb.pgm");
        vector<vector<double> > test;
        vector<double> testvalue;
        testvalue.resize(X.size());
        std::transform(X.begin(), X.end(), mean.begin(), testvalue.begin(), std::minus<double>());
        test.push_back(testvalue);

        vector<vector<double> > omega = getEigenspace(U, test, 1, image_size, kValue);

        set<pair<double, int >> topN;
        set<pair<double, int> >::iterator it;
        pair<double, int> min_dis = {999999, 0};
       // int rank = 50;

        for (int i = 0; i < galery_size; i++) {
            double Mahalanobis_distance = 0;
            for (int j = 0; j < kValue; j++) {
                Mahalanobis_distance += pow(eigespace[i][j] - omega[0][j], 2) / w[j + 1];
            }

            if (topN.size() >= rank) {
                it = prev(topN.end());
                if ((*it).first > Mahalanobis_distance) {
                    topN.erase(it);
                    min_dis = {Mahalanobis_distance, i};
                    topN.insert(min_dis);
                } else
                    min_dis = (*it);
            } else {
                min_dis = {Mahalanobis_distance, i};
                topN.insert(min_dis);
            }
        }

        for (auto itm = topN.begin(); itm != topN.end(); ++itm) {
            if (train_file_list[(*itm).second] == test_file_list[iop]) {
                correctCount++;
                break;
            }
        }
        iop++;
    }
   // cout << "correctCount : " << correctCount << " accuracy ratio : " << (double) correctCount / test_size << endl;
   cout << rank <<","<<(double) correctCount / test_size <<endl;
   return 0;
}

int testing(string train_dir, string test_dir) {
    vector<double> mean = readResults("mean2.csv")[0];
    vector<vector<double> > eigespace = readResults("eigespace2.csv");
    vector<double> w = readResults("eigenvalues2.csv")[0];
    vector<vector<double> > U = readResults("eigenvectors2.csv");
    for(int rank=1;rank<=51;rank+=10)
    testing(mean, eigespace, w,U,train_dir,  test_dir,rank);

    /*1,0.538462
    11,0.777592
    21,0.833612
    31,0.867057
    41,0.882107
    51,0.894649*/
    return 0;
}

int training() {

    vector<vector<double>> images = getX("fa_H");
    int sample_size = images.size();
    cout << sample_size << endl;

    vector<double> mean = getMean(images, sample_size);
    int image_size = mean.size();
    cout << image_size << endl;

    vector<vector<double>> A;
    for (vector<double> samples :  images) {
        std::transform(samples.begin(), samples.end(), mean.begin(), samples.begin(), std::minus<double>());
        A.push_back(samples);
    }
    vector<vector<double> > mul = getCovairiance(A, sample_size, image_size);
    vector<double> w;
    vector<vector<double>> V;
    int adj_size = sample_size + 1;
    w.resize(adj_size);
    V.resize(adj_size);
    for (int i = 1; i <= sample_size; i++)
        V[i].resize(adj_size);

    auto start = high_resolution_clock::now();
    jacobi(mul, sample_size, w, V);
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<seconds>(stop - start);
    cout << "Time taken by function: "
         << duration.count() << " seconds" << endl;

    vector<vector<double>> U = getrealEigenVectors(A, V, sample_size, image_size);
    vector<vector<double> > UT(U[0].size(), vector<double>());
    for (int i = 0; i < U.size(); i++) {
        for (int j = 0; j < U[i].size(); j++) {
            UT[j].push_back(U[i][j]);
        }
    }
    vector<vector<double> > eigespace = getEigenspace(U, A, sample_size, image_size, sample_size);
    // testing(mean, eigespace, w, U);

   /* writeResults(U, "eigenvectors2.csv");
    writeResults(eigespace, "eigespace2.csv");
    writeResults(mean, "mean2.csv");
    writeResults(w, "eigenvalues2.csv");*/
    /*writeImages(mean);
    vitualization(UT);
    writeImages(UT);*/

    return 0;
}
int reconstructionCheck(string train_dir, string test_dir) {
    vector<double> mean = readResults("mean2.csv")[0];
    vector<vector<double> > U = readResults("eigenvectors2.csv");
    int image_size = mean.size();
    int sample_size = U[0].size();

    vector<double> X = getReferenceImageMatrix("fa_H/00001_930831_fa_a.pgm");
    vector<vector<double> > test;
    vector<double> testvalue;
    testvalue.resize(X.size());
    std::transform(X.begin(), X.end(), mean.begin(), testvalue.begin(), std::minus<double>());
    test.push_back(testvalue);

    vector<vector<double>> omega = getEigenspace(U, test, 1, image_size, sample_size);
    pair<int, int> dimen1 = {1, sample_size};
    pair<int, int> dimen2 = { sample_size,image_size};

    vector<double> recon = matrixMultiply(omega,transpose(U),dimen1,dimen2)[0];
    std::transform(recon.begin(), recon.end(), mean.begin(), recon.begin(), std::plus<double>());
    double Mahalanobis_distance = 0;
    for (int j = 0; j < image_size; j++) {
        Mahalanobis_distance += pow(X[j] - recon[j], 2) ;
    }
    cout << Mahalanobis_distance;
    return 0;
}