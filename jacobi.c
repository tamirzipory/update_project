#include "spkmeans.h"

double **jacobi(double **A, int N)
{
    /*Calculates and prints the eigenvalues and eigenvectors as described in 1.2.1*/
    int curr_iter, max_iter;
    double diff;
    double **V;

    curr_iter = 0;
    max_iter = 100;
    V = gen_id_mat(N);
    diff = 1;

    while (diff > eps && curr_iter < max_iter)
    { /*Stops when off(A)^2- - off(A')^2 <= epsilon OR 100 iterations*/
        double off_A = off(A, N);
        curr_iter++;
        A_to_A_tag(A, V, N); /*Also change V during A->A' iteration so that eventually V = P1*P2*P3*...*/
        diff = off_A - off(A, N);
    }
    return V;
}

void A_to_A_tag(double **A, double **V, int N)
{
    /*Calculates A' from A using the relation between them as explained in 1.2.6*/
    int i, j, r;
    int *arr_max;
    double theta, s, t, c, a_ri, a_rj, a_ii, a_jj;

    arr_max = max_indices_off_diag(A, N); /*Finding pivot - 1.2.3*/
    i = arr_max[0];
    j = arr_max[1];
    theta = (A[j][j] - A[i][i]) / (2 * A[i][j]);
    t = sign(theta) / (fabs(theta) + sqrt((pow(theta, 2)) + 1));
    c = 1 / sqrt((pow(t, 2)) + 1);
    s = t * c;
    V_multi_P(V, s, c, N, i, j); /*V = VP*/

    for (r = 0; r < N; r++)
    {
        if ((r != j) && (r != i))
        {
            a_ri = c * A[r][i] - s * A[r][j];
            a_rj = c * A[r][j] + s * A[r][i];
            A[r][i] = a_ri;
            A[r][j] = a_rj;
            A[j][r] = a_rj;
            A[i][r] = a_ri;
        }
    }
    /*1.2.3*/
    a_ii = pow(c, 2) * A[i][i] + pow(s, 2) * A[j][j] - 2 * s * c * A[i][j];
    a_jj = pow(c, 2) * A[j][j] + pow(s, 2) * A[i][i] + 2 * s * c * A[i][j];
    A[j][j] = a_jj;
    A[i][i] = a_ii;
    A[i][j] = 0.0;
    A[j][i] = 0.0;

    free(arr_max);
}