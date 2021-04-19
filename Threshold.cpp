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
static unsigned long jacobi_max_iterations = 500;

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
vector<vector<double>> getX() {
    vector<vector<double>> images;
    vector<string> file_list = listFile();
    int k = 300;
    for (string fn :  file_list) {
        vector<double> ref_image = getReferenceImageMatrix(fn);
        images.push_back(ref_image);
        k--;
        if (k < 0)
            break;
    }
    return images;
}

vector<double> getMean(vector<vector<double>> images, int sample_size) {

    vector<double> mean;

    for (vector<double> samples :  images) {

        if (mean.empty())
            mean = samples;
        else
            std::transform(mean.begin(), mean.end(), samples.begin(), mean.begin(), std::plus<double>());
    }


    std::transform(mean.begin(), mean.end(), mean.begin(),
                   [sample_size](double c) { return (double) c / (double) sample_size; });
    return mean;
}

vector<vector<double> > getCovairiance(vector<vector<double>> A, int sample_size, int image_size) {
    vector<vector<double> > AT(A[0].size(), vector<double>());

    for (int i = 0; i < A.size(); i++) {
        for (int j = 0; j < A[i].size(); j++) {
            AT[j].push_back(A[i][j]);
        }
    }
    //displayMatrix(AT);
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
            mul[i][j] = mul[i][j] / sample_size;
        }
    }
    return mul;
}
template<typename T>
struct square
{
    T operator()(const T& Left, const T& Right) const
    {
        return (Left+ Right*Right  );
        //fabs(
    }
};
void normalize(vector<vector<double>> &A) {
    std::for_each(std::begin(A), std::end(A),
                  [](vector<double> (&row)) {
                      double sum=std::accumulate(std::begin(row), std::end(row), 0.0,square<double>());
                      std::transform(std::begin(row), std::end(row), std::begin(row),
                                     [sum](double x) { return  x/pow(sum,0.5); });
                  });
}

vector<vector<double>> getrealEigenVectors(vector<vector<double>> A,vector<vector<double>> V, int sample_size, int image_size)
{
    vector<vector<double> > AT(A[0].size(), vector<double>());

    for (int i = 0; i < A.size(); i++) {
        for (int j = 0; j < A[i].size(); j++) {
            AT[j].push_back(A[i][j]);
        }
    }
    pair<int, int> dimen1 = {image_size,sample_size  };
    pair<int, int> dimen2 = {sample_size, sample_size};
    vector<vector<double> > mul;
    /*cout << AT.size() << " " << AT[0].size() << endl;
    cout <<sample_size << " " <<image_size << endl;
    cout << V.size() << " " << V[1].size() << endl;*/
    mul.resize(dimen1.first);
    for (int i = 0; i < dimen1.first; i++) {
        mul[i].resize(dimen2.second);
        for (int j = 0; j < dimen2.second; j++) {
            mul[i][j] = 0;
            //cout << AT[i][j] << " ";
            for (int k = 0; k < dimen2.first; k++) {
                mul[i][j] += AT[i][k]
                             * (V[k+1][j+1]);
            }
            //mul[i][j] = mul[i][j] / sample_size;
        }
    }
   // displayMatrix(mul);
    //cout << "why?";
    normalize(mul);
    return mul;
}
vector<vector<double> > getEigenfaces(vector<vector<double>> U,vector<vector<double>> X, int sample_size, int image_size,int eigensize) {
    vector<vector<double> > UT(U[0].size(), vector<double>());

    for (int i = 0; i < U.size(); i++) {
        for (int j = 0; j < U[i].size(); j++) {
            UT[j].push_back(U[i][j]);
        }
    }
   // displayMatrix(U);
    pair<int, int> dimen1 = {eigensize,sample_size};
    pair<int, int> dimen2 = {sample_size, image_size};
    vector<vector<double> > mul;
    mul.resize(dimen1.first);
    for (int i = 0; i < dimen1.first; i++) {
        mul[i].resize(dimen2.second);
        for (int j = 0; j < dimen2.second; j++) {
            mul[i][j] = 0;
            for (int k = 0; k < dimen2.first; k++) {
                mul[i][j] += U[i][k]
                             * (X[k][j]);
            }

        }
    }
    return mul;
}
int main(int argc, char *argv[]) {

    // get x
    vector<vector<double>> images = getX();
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
   // displayMatrix(mul);

    //cout << endl << sample_size << endl;
     vector<double> w;
     vector<vector<double>> V;
     int adj_size = sample_size + 1;
     w.resize(adj_size);
     V.resize(adj_size);
     for (int i = 1; i <= sample_size; i++)
         V[i].resize(adj_size);


    jacobi2(mul, sample_size, w, V);

    //displayMatrix(V);

    vector<vector<double>> realEigenVectors = getrealEigenVectors(A,V,sample_size,image_size);

    //displayMatrix(realEigenVectors);

    vector<vector<double> > eigenfaces = getEigenfaces(realEigenVectors,images,sample_size,image_size,10);

    //displayMatrix(eigenfaces);
    int N = 20, M = 16;
    int fcno=0;
    for( vector<double> eg : eigenfaces){
        fcno++;
    ImageType image(N, M, 255);
    for (int i = 0; i < N; i++)
        for (int j = 0; j < M; j++) {
            image.setPixelVal(i, j,eg[i * M + j]);
        }
    string filename=to_string(fcno)+string("eigenfaces.pgm");
    writeImage( filename , image);
    }
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
    /*SM = 10;
    SN = 1;*/
    // cout << SN<<","<< SM<<endl;
    vector<double> ref;
    for (int i = 0; i < SN; i++)
        for (int j = 0; j < SM; j++) {
            Simage.getPixelVal(i, j, sval);
            ref.push_back(sval);
        }
    return ref;
}

