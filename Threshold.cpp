#include <iostream>
#include <math.h>
#include "image.h"
#include "MatricOP.h"
#include "DiscriminantCases.h"
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
using namespace std::chrono;

void readImageHeader(string, int &, int &, int &, bool &);

void readImage(string, ImageType &);

void writeImage(string, ImageType &);


void displayMatrix(vector<double> vector);

void displayMatrix(vector<vector<double>> matrix);

vector<double> getReferenceImageMatrix(string sample_file);

vector<string> listFile();

vector<vector<double>> getX();

vector<double> getMean(vector<vector<double>> images, int sample_size);

vector<vector<double> > getCovairiance(vector<vector<double>> A, int sample_size, int image_size);

int jacobi2(vector<vector<double>> zas, dimension n, vector<double> &w, vector<vector<double>> &V);

vector<vector<double>>
getrealEigenVectors(vector<vector<double>> A, vector<vector<double>> V, int sample_size, int image_size);

void normalize(vector<vector<double>> &A);

vector<vector<double> >
getEigenspace(vector<vector<double>> U, vector<vector<double>> X, int sample_size, int image_size, int eigensize);

void writeImages(vector<vector<double>> eigenfaces);

void writeImages(vector<double> meanfaces);

void vitualization(vector<vector<double>> &A);

int training();
int testing(vector<double> mean, vector<vector<double> > eigespace, vector<double> X, vector<double> w,
            vector<vector<double> > U);

int main(int argc, char *argv[]) {
    // get x
    training();

    return 0;
}

template<typename T>
struct Mahalanobis {
    T operator()(const T &Left, const T &Right) const {
        return (Left / pow(Right, 0.5));
    }
};

int testing(vector<double> mean, vector<vector<double> > eigespace, vector<double> Y, vector<double> w,
            vector<vector<double> > U) {

    vector<double> X = getReferenceImageMatrix("fb_l/00019_940422_fb.pgm");
    vector<vector<double> > test;
    // test.resize(1);
    vector<double> testvalue;
    testvalue.resize(X.size());
    // displayMatrix(test);
    std::transform(X.begin(), X.end(), mean.begin(), testvalue.begin(), std::minus<double>());
    //displayMatrix(testvalue);
    test.push_back(testvalue);
   // displayMatrix(test);
    vector<double> testvalueP = test[0];
   // vector<double> testvalueP = eigespace[60];

    vector<vector<double> > omega = getEigenspace(U, test, 1, 320, 10);

    int image_size=eigespace[0].size();

    pair<int,double> min_dis={0,999999};
   // cout <<image_size <<endl;
    for (int i=0 ;i < eigespace.size() ; i ++) {

        double Mahalanobis_distance =0;
        for(int j=0 ;j < image_size ; j ++){

            Mahalanobis_distance = pow(eigespace[i][j]-omega[0][j],2) /*/ w[i+1]*/;
        }
        //cout << i <<","<< Mahalanobis_distance << endl;
        if(min_dis.second > Mahalanobis_distance)
            min_dis={i,Mahalanobis_distance};



    }

cout<< min_dis.first+1<< " with "<< min_dis.second<<endl;
    vector<string> file_list = listFile();
    cout << file_list[min_dis.first] ;
    return 0;
}

int training() {

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
    int Kvalue = 10;
    vector<vector<double> > eigespace = getEigenspace(U, A, sample_size, image_size, Kvalue);

    vector<vector<double> > UT(U[0].size(), vector<double>());
    for (int i = 0; i < U.size(); i++) {
        for (int j = 0; j < U[i].size(); j++) {
            UT[j].push_back(U[i][j]);
        }
    }
    // //displayMatrix(U);


    //displayMatrix(eigenfaces);
    /*writeImages(mean);
    vitualization(UT);
    writeImages(UT);*/

    testing(mean, eigespace, images[120], w, U);

    return 0;
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


vector<vector<double>> getX() {
    vector<vector<double>> images;
    vector<string> file_list = listFile();
    int k = 700;
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
    mean.resize(images[0].size());

    for (vector<double> samples :  images) {
        std::transform(mean.begin(), mean.end(), samples.begin(), mean.begin(), std::plus<double>());
    }

    std::transform(mean.begin(), mean.end(), mean.begin(),
                   [sample_size](double c) { return (double) c / (double) sample_size; });
    return mean;
};

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
struct square {
    T operator()(const T &Left, const T &Right) const {
        return (Left + Right * Right);
    }
};

template<typename T>
struct Normalize {
    T operator()(const T &Left, const T &Right) const {
        return (Left / pow(Right, 0.5));
    }
};

void normalize(vector<vector<double>> &A) {
    vector<double> length;
    length.resize(A[0].size());

    for (vector<double> samples :  A) {
        std::transform(length.begin(), length.end(), samples.begin(), length.begin(), square<double>());
    }
    //displayMatrix(length);
    std::for_each(std::begin(A), std::end(A),
                  [length](vector<double> (&row)) {
                      std::transform(std::begin(row), std::end(row), std::begin(length), std::begin(row),
                                     Normalize<double>());
                  });
}

void vitualization(vector<vector<double>> &A) {

    std::for_each(std::begin(A), std::end(A),
                  [](vector<double> (&row)) {
                      pair<vector<double>::iterator, vector<double>::iterator> min_max = minmax_element(row.begin(),
                                                                                                        row.end());
                      double min = row[min_max.first - row.begin()];
                      double max = row[min_max.second - row.begin()];
                      std::transform(std::begin(row), std::end(row), std::begin(row),
                                     [min, max](double x) { return 255 * (x - min) / (max - min); });
                  });
}


vector<vector<double>>
getrealEigenVectors(vector<vector<double>> A, vector<vector<double>> V, int sample_size, int image_size) {
    vector<vector<double> > AT(A[0].size(), vector<double>());

    for (int i = 0; i < A.size(); i++) {
        for (int j = 0; j < A[i].size(); j++) {
            AT[j].push_back(A[i][j]);
        }
    }
    pair<int, int> dimen1 = {image_size, sample_size};
    pair<int, int> dimen2 = {sample_size, sample_size};
    vector<vector<double> > mul;
    mul.resize(dimen1.first);
    for (int i = 0; i < dimen1.first; i++) {
        mul[i].resize(dimen2.second);
        for (int j = 0; j < dimen2.second; j++) {
            mul[i][j] = 0;
            for (int k = 0; k < dimen2.first; k++) {
                mul[i][j] += AT[i][k]
                             * (V[k + 1][j + 1]);
            }
        }
    }
    //displayMatrix(mul);
    normalize(mul);
    return mul;
}

vector<vector<double> >
getEigenspace(vector<vector<double>> U, vector<vector<double>> X, int sample_size, int image_size, int eigensize) {

    pair<int, int> dimen1 = {sample_size, image_size};
    pair<int, int> dimen2 = {image_size, eigensize};
    vector<vector<double> > mul;
    mul.resize(dimen1.first);
    for (int i = 0; i < dimen1.first; i++) {
        mul[i].resize(dimen2.second);
        for (int j = 0; j < dimen2.second; j++) {
            mul[i][j] = 0;
            for (int k = 0; k < dimen2.first; k++) {
                mul[i][j] += X[i][k]
                             * (U[k][j]);
            }
        }
    }
    return mul;
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