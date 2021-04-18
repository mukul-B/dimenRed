#include <iostream>




/* JACOBI.C  module

   NOTES: The jacobi() function is a modified version of the one found in
          'Numerical Recipes in C'.



*/



#include <math.h>
#include <stdio.h>
#include <vector>


typedef unsigned dimension;
typedef unsigned iterations;
#define ROTATE(S, i, j, k, l) g=S[i][j];h=S[k][l];S[i][j]=g-s*(h+g*tau); \
              S[k][l]=h+s*(g-h*tau)


/* Maximum number of iterations allowed in jacobi() */
static unsigned long jacobi_max_iterations = 500;

using namespace std;

int printsol(dimension n, vector<double> w, vector<vector<double>> V) {

    for (int l = 1; l <= n; l++) {
        for (int m = 1; m <= n; m++)
            cout << V[l][m] << " ";
        cout << endl;
    }

    for (int m = 1; m <= n; m++)
        cout << w[m] << endl;
    return 0;
}

int jacobi(vector<vector<double>> zas, dimension n, vector<double> &w, vector<vector<double>> &V) {
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

int mainjl() {

   /* vector<vector<double>> zas = {{420.222,  694.556,  -1114.78},
                                  {694.556,  1154.89,  -1849.44},
                                  {-1114.78, -1849.44, 2964.22}};

    vector<double> w;
    vector<vector<double>> V;
    int n = 3;
    int  ad=n+1;
    iterations i, j;

    w.resize(ad);
    V.resize(ad);
    for (i = 1; i <= n; i++)
        V[i].resize(ad);

    for (int l = 0; l < n; l++) {
        for (int m = 0; m < n; m++)
            cout << zas[l][m] << " ";
        cout << endl;
    }
    printsol(n, w, V);
    jacobi(zas, n, w, V);*/

   int sample_size=3;
    vector<double> w;
    vector<vector<double>> V;
    int adj_size = sample_size + 1;
    w.resize(adj_size);
    V.resize(adj_size);
    for (int i = 1; i <= sample_size; i++)
        V[i].resize(adj_size);
    vector<vector<double>> zas = {{420.222,  694.556,  -1114.78},
                                  {694.556,  1154.89,  -1849.44},
                                  {-1114.78, -1849.44, 2964.22}};

    jacobi(zas, sample_size, w, V);
    printsol(sample_size, w, V);
    return 0;
}