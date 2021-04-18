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
int printsol2(double **S,dimension n,double *w,double **V) {

    for (int l = 1; l <= n; l++)
    {for (int m = 1; m <= n; m++)
            cout<< S[l][m] << ",";
        cout<<endl;}

    for (int l = 1; l <= n; l++)
    {for (int m = 1; m <= n; m++)
            cout<< V[l][m] << ",";
        cout<<endl;}

    for (int m = 1; m <= n; m++)
        cout<< w[m] << endl;
    return 0;
}
int jacobi2(double **S,dimension n,double *w,double **V) {

    iterations i, j, k, iq, ip;
    double tresh, theta, tau, t, sm, s, h, g, c;
    double p;
    double *b;
    double *z;
    int nrot;




    b = (double *) malloc(n * sizeof(double));
    b--;
    z = (double *) malloc(n * sizeof(double));
    z--;

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
        for (ip = 1; ip <= n-1 ; ip++) {
            for (iq = ip + 1; iq <= n; iq++) {sm += fabs(S[ip][iq]);
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
            z++;
            free(z);
            b++;
            free(b);
            //break;

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


    free(z++);
    free(b++);
    return (-1);
}

int main2() {

    double zas[2][2] = {{3, -2},
                        {-2, 6}};

    dimension n=2;
    double **S ;
    double *w;
    double **V;
    iterations i, j, k, iq, ip;
    S = (double **)malloc(n * sizeof(double *));
    for (i=0; i<=n; i++)
        S[i] = (double *)malloc(n * sizeof(double));

    for (i = 1; i <= n; i++)
        for (j = 1; j <= n; j++)
            S[i][j] = zas[i-1][j-1];

    V = (double **)malloc(n * sizeof(double *));
    for (i=0; i<=n; i++)
        V[i] = (double *)malloc(n * sizeof(double));

    w = (double *) malloc(n * sizeof(double));

    jacobi2(S,n,w,V);
    printsol2(S,n,w,V);
    return 0;
}