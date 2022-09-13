#include "spkmeans.h"

/*
Musn't to compile it! it's only fot the organization of the code!
*/

double **jacobi(double **A, int N){
    /*Calculates and prints the eigenvalues and eigenvectors*/
    int curr_iter, max_iter;
    double diff;
    double **V;
    curr_iter = 0, max_iter = 100;
    V = create_I_mat(N);
    diff = 1;
    while (diff > eps && curr_iter < max_iter){ 
        double off_A = off(A, N);
        curr_iter++;
        calc_A_tag(A, V, N); 
        diff = off_A - off(A, N);
    }
    return V;
}

void calc_A_tag(double **A, double **V, int N){
    /*calc A tag*/
    int i, j, r;
    int *arr_max;
    double theta, s, t, c, var_1, var_2, var_3, var_4;

    arr_max = max_indices_off_diag(A, N); 
    i = arr_max[0], j = arr_max[1];
    theta = (A[j][j] - A[i][i]) / (2 * A[i][j]);
    t = sign(theta) / (fabs(theta) + sqrt((pow(theta, 2)) + 1));
    c = 1 / sqrt((pow(t, 2)) + 1);
    s = t * c;
    V_kaful_P(V, s, c, N, i, j); /*do the mult*/
    for (r = 0; r < N; r++){
        if ((r != j) && (r != i)){
            var_1 = c * A[r][i] - s * A[r][j];
            var_2 = c * A[r][j] + s * A[r][i];
            A[r][i] = var_1;
            A[r][j] = var_2;
            A[j][r] = var_2;
            A[i][r] = var_1;
        }
    }
    var_3 = pow(c, 2) * A[i][i] + pow(s, 2) * A[j][j] - 2 * s * c * A[i][j];
    var_4 = pow(c, 2) * A[j][j] + pow(s, 2) * A[i][i] + 2 * s * c * A[i][j];
    A[j][j] = var_4;
    A[i][i] = var_3;
    A[i][j] = 0.0;
    A[j][i] = 0.0;

    free(arr_max);
}
