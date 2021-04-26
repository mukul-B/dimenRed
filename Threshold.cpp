#include <iostream>
#include <math.h>
#include<sstream>
#include "jacobi.h"
#include "lamdafunctions.h"
#include "bits/stdc++.h"
#include "matrixOperations.h"
#include "matrixfunctions.h"

#define lowReso 0
#define highReso 1

#define testPhase 0
#define trainPhase 1

using namespace std;
using namespace std::chrono;

int training(string train_dir, string resolution);

int testing(vector<double> mean, vector<vector<double> > eigespace, vector<double> w,
            vector<vector<double> > U, double preservationRatio, int rank, bool detail);

int testing(string train_dir, string test_dir, string resolution);

int reconstructionCheck(string train_dir, string test_dir);

int main(int argc, char *argv[]) {
   /* int prob;

    printf("Programing Assignment 2:  question 3\n");
    printf("------------------------------------------\n");
    printf("Main Menu\n");
    printf("1.  Training process with High resolution Images .\n");
    printf("2.  Testing process with High resolution Images \n");
    printf("3.  Training process with low resolution Images.\n");
    printf("4.  Testing process with low resolution Images\n");
    printf(" Please enter an option from the main menu: ");

    fflush(stdin);
    cin >> prob;*/
    string training_folder,testing_folder,resolution;

    int high_low_resolution = highReso;
    int training_testing = testPhase;

    if(high_low_resolution==lowReso){
        training_folder = "fa_L";
        testing_folder = "fb_L";
        resolution="L";
    }
    else if(high_low_resolution==highReso){
        training_folder = "fa_H";
        testing_folder = "fb_H";
        resolution="H";
    }
    if(training_testing==trainPhase)
        training(training_folder,resolution);
    else if(training_testing==testPhase)
        testing(training_folder, testing_folder, resolution);

    //reconstructionCheck(training_folder, testing_folder);

    return 0;
}


int testing(vector<double> mean, vector<vector<double> > eigespace, vector<double> w,
            vector<vector<double> > U, string train_dir, string test_dir, double preservationRatio, int rank,
            bool detail = false) {

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
    int kValue = 0;
    double eigenSum = 0;
    double totalEigenSum = 0;

    for (int i = 0; i < galery_size; i++) {
        totalEigenSum += w[i + 1];
    }
    for (int i = 0; i < galery_size; i++) {
        eigenSum += w[i + 1];
        if ((eigenSum / totalEigenSum) > preservationRatio) {
            kValue = i - 1;
            break;
        }
    }
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
        int found = 0;

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
                found = 1;
                break;
            }
        }
        if (detail) {
            string matchNomatch;
            if (found == 1)
                 matchNomatch = " correctly";
            else
            matchNomatch = " incorrectly";

                cout << test_file_list[iop] << " at " << iop << matchNomatch <<" matched with " <<train_file_list[min_dis.second] << " at " <<  min_dis.second << endl;
        }


        iop++;
    }
    cout << rank << "," << (double) correctCount / test_size << endl;
    return 0;
}

int testing(string train_dir, string test_dir, string resolution) {
   vector<double> mean = readResults(resolution+"mean.csv")[0];
    vector<vector<double> > eigespace = readResults(resolution+"eigespace.csv");
    vector<double> w = readResults(resolution+"eigenvalues.csv")[0];
    vector<vector<double> > U = readResults(resolution+"eigenvectors.csv");
    /*vector<double> mean = readResults("mean2.csv")[0];
    vector<vector<double> > eigespace = readResults("eigespace2.csv");
    vector<double> w = readResults("eigenvalues2.csv")[0];
    vector<vector<double> > U = readResults("eigenvectors2.csv");*/
   // testing(mean, eigespace, w, U, train_dir, test_dir, 0.8, 1, true);

   cout << 0.8 <<endl;
    for(int rank=1;rank<=51;rank+=5)
    testing(mean, eigespace, w,U,train_dir,  test_dir,0.8,rank);
    cout << 0.9 <<endl;
    for(int rank=1;rank<=51;rank+=5)
        testing(mean, eigespace, w,U,train_dir,  test_dir,0.9,rank);
    cout << 0.95 <<endl;
    for(int rank=1;rank<=51;rank+=5)
        testing(mean, eigespace, w,U,train_dir,  test_dir,0.95,rank);

    return 0;
}

int training(string train_dir, string resolution) {

    vector<vector<double>> images = getX(train_dir);
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

    writeResults(U, resolution + "eigenvectors.csv");
    writeResults(eigespace, resolution + "eigespace.csv");
    writeResults(mean, resolution + "mean.csv");
    writeResults(w, resolution + "eigenvalues.csv");
    /* writeImages(mean,resolution);
     vitualization(UT);
     writeImages(UT,resolution);*/

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
    pair<int, int> dimen2 = {sample_size, image_size};

    vector<double> recon = matrixMultiply(omega, transpose(U), dimen1, dimen2)[0];
    std::transform(recon.begin(), recon.end(), mean.begin(), recon.begin(), std::plus<double>());
    double Mahalanobis_distance = 0;
    for (int j = 0; j < image_size; j++) {
        Mahalanobis_distance += pow(X[j] - recon[j], 2);
    }
    cout << Mahalanobis_distance;
    return 0;
}