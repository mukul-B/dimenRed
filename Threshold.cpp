#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include "image.h"
#include "MatricOP.h"
#include "DiscriminantCases.h"

#include<stdio.h>
#include<string.h>
#include<sstream>
#include<dirent.h>


typedef unsigned dimension;
typedef unsigned iterations;
#define ROTATE(S, i, j, k, l) g=S[i][j];h=S[k][l];S[i][j]=g-s*(h+g*tau); \
              S[k][l]=h+s*(g-h*tau)


/* Maximum number of iterations allowed in jacobi() */
static unsigned long jacobi_max_iterations = 3;

using namespace std;

vector<string> listFile() {
    ifstream inn;
    DIR *pDIR;
    vector<string> file_list;
    struct dirent *entry;
    if (pDIR = opendir("./fa_L")) {
        while (entry = readdir(pDIR)) {
            if (strcmp(entry->d_name, ".") != 0 && strcmp(entry->d_name, "..") != 0) {
                string str(entry->d_name);
                file_list.emplace_back("fa_L/" + str);
            }

        }
        closedir(pDIR);
    }
    return file_list;

}

int jacobi2(vector<vector<double>> zas, dimension n, vector<double> &w, vector<vector<double>> &V) {
    iterations i, j, k, iq, ip;
    double tresh, theta, tau, t, sm, s, h, g, c;
    double p;
    int nrot;
    vector<double> b, z;
    int m = n + 1;
    b.resize(m);
    z.resize(m);
    vector<vector<double>> S;
    S.resize(m);
    for (i = 1; i <= n; i++) {
        S[i].resize(m);
        for (j = 1; j <= n; j++) {
            S[i][j] = zas[i - 1][j - 1];
        }
    }
    for (ip = 1; ip <= n; ip++) {
        for (iq = 1; iq <= n; iq++)
            V[ip][iq] = 0.0;
        V[ip][ip] = 1.0;
    }

    for (ip = 1; ip <= n; ip++) {
        b[ip] = w[ip] = S[ip][ip];
        z[ip] = 0.0;
    }

    nrot = 0;

    for (i = 1; i <= jacobi_max_iterations; i++) {
        sm = 0.0;
        for (ip = 1; ip <= n - 1; ip++) {
            for (iq = ip + 1; iq <= n; iq++) {
                sm += fabs(S[ip][iq]);
            }
        }
        if (sm == 0.0) {
            for (i = 1; i < n; i++) {
                p = w[k = i];
                for (j = i + 1; j <= n; j++) if (w[j] >= p) p = w[k = j];
                if (k != i) {
                    w[k] = w[i];
                    w[i] = p;
                    for (j = 1; j <= n; j++) {
                        p = V[j][i];
                        V[j][i] = V[j][k];
                        V[j][k] = p;
                    }
                }
            }


            for (i = 2; i <= n; i++) {
                for (j = 1; j < i; j++) S[j][i] = S[i][j];
            }
            return (nrot);
        }
        if (i < 4) tresh = 0.2 * sm / (n * n); else tresh = 0.0;
        for (ip = 1; ip <= n - 1; ip++) {
            for (iq = ip + 1; iq <= n; iq++) {
                g = 100.0 * fabs(S[ip][iq]);
                if (i > 4 && fabs(w[ip]) + g == fabs(w[ip]) && fabs(w[iq]) + g == fabs(w[iq]))
                    S[ip][iq] = 0.0;
                else if (fabs(S[ip][iq]) > tresh) {
                    h = w[iq] - w[ip];
                    if (fabs(h) + g == fabs(h))
                        t = (S[ip][iq]) / h;
                    else {
                        theta = 0.5 * h / (S[ip][iq]);
                        t = 1.0 / (fabs(theta) + sqrt(1.0 + theta * theta));
                        if (theta < 0.0) t = -t;
                    }
                    c = 1.0 / sqrt(1 + t * t);
                    s = t * c;
                    tau = s / (1.0 + c);
                    h = t * S[ip][iq];
                    z[ip] -= h;
                    z[iq] += h;
                    w[ip] -= h;
                    w[iq] += h;
                    S[ip][iq] = 0.0;
                    for (j = 1; j <= ip - 1; j++) {
                        ROTATE(S, j, ip, j, iq);
                    }
                    for (j = ip + 1; j <= iq - 1; j++) {
                        ROTATE(S, ip, j, j, iq);
                    }
                    for (j = iq + 1; j <= n; j++) {
                        ROTATE(S, ip, j, iq, j);
                    }
                    for (j = 1; j <= n; j++) {
                        ROTATE(V, j, ip, j, iq);
                    }
                    ++nrot;
                }
            }
        }
        for (ip = 1; ip <= n; ip++) {
            b[ip] += z[ip];
            w[ip] = b[ip];
            z[ip] = 0.0;
        }
    }
    return (-1);
}

void readImageHeader(string, int &, int &, int &, bool &);

void readImage(string, ImageType &);

void writeImage(string, ImageType &);

vector<double> getReferenceImageMatrix(string);


int main(int argc, char *argv[]) {

    int M, N, Q;
    bool type;
    vector<string> file_list = listFile();
    vector<vector<double>> images;
    vector<double> mean;
    int k = 2;
    for (string fn :  file_list) {
        vector<double> ref_image = getReferenceImageMatrix(fn);
        images.push_back(ref_image);
       /* k--;
        if (k < 0)
            break;*/
    }

    for (vector<double> samples :  images) {

        if (mean.empty())
            mean = samples;
        else
            std::transform(mean.begin(), mean.end(), samples.begin(), mean.begin(), std::plus<double>());
    }

    int sample_size = images.size();
    int image_size = mean.size();
    cout << image_size;
    std::transform(mean.begin(), mean.end(), mean.begin(),
                   [sample_size](double c) { return (double) c / (double) sample_size; });

    vector<vector<double>> A;
    for (vector<double> samples :  images) {
        std::transform(samples.begin(), samples.end(), mean.begin(), samples.begin(), std::minus<double>());
        A.push_back(samples);
    }

    vector<vector<double> > AT(A[0].size(), vector<double>());

    for (int i = 0; i < A.size(); i++) {
        for (int j = 0; j < A[i].size(); j++) {
            AT[j].push_back(A[i][j]);
        }
    }

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
            mul[i][j]=mul[i][j]/sample_size;
        }
    }

//    std::transform(mul.begin(), mul.end(), mul.begin(),[sample_size](double c) { return (double) c / (double) sample_size; });

    cout << endl << sample_size << endl;
    vector<double> w;
    vector<vector<double>> V;
    int adj_size = sample_size + 1;
    w.resize(adj_size);
    V.resize(adj_size);
    for (int i = 1; i <= sample_size; i++)
        V[i].resize(adj_size);

     jacobi2(mul, sample_size, w, V);

    /* cout << endl << V.size() << endl;*/

   /* for (vector<double> samples :  images) {


        for (double j :  samples) {

            cout << j << " ";

        }
        cout << "S \n";
    }
    for (double j :  mean) { cout << j << " "; }
    cout << "Meam \n";

    for (vector<double> samples :  A) {

        for (double j :  samples) {

            cout << j << " ";

        }
        cout << "A \n";
    }

    for (vector<double> samples :  AT) {


        for (double j :  samples) {

            cout << j << " ";

        }
        cout << "AT \n";


    }

    for (vector<double> sampl :  mul) {
        for (double j :  sampl) {
            cout << j << " ";
        }
        cout << "covariance \n";
    }
    cout << mul.size() << endl;
    for (vector<double> sampl :  V) {
        for (double j :  sampl) {
            cout << j << " ";
        }
        cout << "V \n";
    }
*/
    for (double j :  w) {

        cout << j << " ";

    }
    cout << "lamda \n";
    return 0;
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
    vector<double> ref;
    for (int i = 0; i < SN; i++)
        for (int j = 0; j < SM; j++) {
            Simage.getPixelVal(i, j, sval);
            ref.push_back(sval);
        }
    return ref;
}

