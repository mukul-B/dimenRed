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
            vector<vector<double> > U);

int main(int argc, char *argv[]) {
    training();
    return 0;
}


int testing(vector<double> mean, vector<vector<double> > eigespace, vector<double> w,
            vector<vector<double> > U) {

    vector<string> file_list = listFile("fa_L");
    vector<string> file_list2 = listFile("fb_L");
    std::transform(file_list2.begin(), file_list2.end(), file_list2.begin(), findimageId<string>());
    std::transform(file_list.begin(), file_list.end(), file_list.begin(), findimageId<string>());

    vector<vector<double>> testImages = getX("fb_L");
    int image_size = eigespace[0].size();
    int galery_size = eigespace.size();
    cout << galery_size << endl;

    int iop = 0;

    int TP = 0, TN = 0;

    for (vector<double> X : testImages) {
//        vector<double> X = getReferenceImageMatrix("fb_l/00019_940422_fb.pgm");
        vector<vector<double> > test;
        vector<double> testvalue;
        testvalue.resize(X.size());
        std::transform(X.begin(), X.end(), mean.begin(), testvalue.begin(), std::minus<double>());
        test.push_back(testvalue);
        vector<vector<double> > omega = getEigenspace(U, test, 1, 320, 190);

        set<pair<double, int >> topN;
        set<pair<double, int> >::iterator it;
        int N = 50;
        pair<double, int> min_dis = {999999, 0};

        //pair< double,int> min_dis = { 999999,0};

        for (int i = 0; i < galery_size; i++) {
            double Mahalanobis_distance = 0;
            for (int j = 0; j < image_size; j++) {
                Mahalanobis_distance = pow(eigespace[i][j] - omega[0][j], 2) / w[i + 1];
            }

            if (topN.size() >= N) {
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

            /*if (min_dis.first > Mahalanobis_distance)
            { min_dis = { Mahalanobis_distance,i};
            }*/

        }

        for (auto itm = topN.begin(); itm != topN.end(); ++itm) {
            //cout << (*itm).first << " ";

            if (file_list[(*itm).second] == file_list2[iop]) {
                TP++;
                continue;
            }

        }
        //cout << min_dis.first + 1 << " with " << min_dis.second << endl;

        //cout << file_list[min_dis.first] << " <> " << file_list2[iop] << endl;
        iop++;
    }
    cout << "TP : " << TP << " TN : " << TN << endl;
    return 0;
}

int training() {
    vector<vector<double>> images = getX("fa_L");
    int sample_size = images.size();
    cout << sample_size;

    //displayMatrix(images);
    //get Mean
    vector<double> mean = getMean(images, sample_size);
    int image_size = mean.size();
    cout << image_size;
    //displayMatrix(mean);
    // get A
    vector<vector<double>> A;
    for (vector<double> samples :  images) {
        std::transform(samples.begin(), samples.end(), mean.begin(), samples.begin(), std::minus<double>());
        A.push_back(samples);
    }
    //displayMatrix(A);
// get coviance
    vector<vector<double> > mul = getCovairiance(A, sample_size, image_size);
    //displayMatrix(mul);
    //cout << endl << sample_size << endl;
    vector<double> w;
    vector<vector<double>> V;
    int adj_size = sample_size + 1;
    w.resize(adj_size);
    V.resize(adj_size);
    for (int i = 1; i <= sample_size; i++)
        V[i].resize(adj_size);
    auto start = high_resolution_clock::now();
    // Call the function, here sort()sort(values.begin(), values.end());
    jacobi2(mul, sample_size, w, V);
    // Get ending timepoint
    auto stop = high_resolution_clock::now();
    // Get duration. Substart timepoints to
    // get durarion. To cast it to proper unit
    // use duration cast method
    auto duration = duration_cast<microseconds>(stop - start);
    cout << "Time taken by function: "
         << duration.count() << " microseconds" << endl;
    //displayMatrix(w);
    //displayMatrix(V);
    vector<vector<double>> U = getrealEigenVectors(A, V, sample_size, image_size);
    //displayMatrix(realEigenVectors);
    // int Kvalue = image_size;
    vector<vector<double> > UT(U[0].size(), vector<double>());
    for (int i = 0; i < U.size(); i++) {
        for (int j = 0; j < U[i].size(); j++) {
            UT[j].push_back(U[i][j]);
        }
    }
    // //displayMatrix(U);
    cout << UT.size() << " " << UT[0].size() << endl;
    //displayMatrix(eigenfaces);
    /*writeImages(mean);
    vitualization(UT);
    writeImages(UT);*/
    int kValue;
    double eigenSum = 0;
    double totalEigenSum = 0;
    for (int i = 0; i < sample_size; i++) {
        totalEigenSum += U[i][i];
    }

    for (int i = 0; i < sample_size; i++) {
        if ((eigenSum/totalEigenSum) > 0.9) {
            kValue = i-1;
            break;
        } else {
            eigenSum += U[i][i];
        }
    }

    vector<vector<double> > eigespace = getEigenspace(U, A, sample_size, image_size, kValue);
    testing(mean, eigespace, w, U);
    return 0;
}