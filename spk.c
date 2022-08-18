#define eps pow(10, -5)
#include "spkmeans.h"

int compare_eigen(const void *a, const void *b)
{
    struct pair *A, *B;
    A = (pair *)a;
    B = (pair *)b;
    return A->eigenValue == B->eigenValue ? A->index-B->index : A->eigenValue > B-> eigenValue ? 1 : -1;
}



double ** init_double_mat(int rows){
    int i;
    double** ret = (double**) calloc(rows, rows*sizeof(double));
    print_err_with_assertg(ret != NULL);
    for(i = 0; i < rows; i++){
        ret[i] = (double *) calloc(rows, sizeof(double));
        print_err_with_assertg(ret[i] != NULL);
    }
    return ret;
}

void sortpairs()
{
    int i;
    pairs = (pair *)calloc(vector_num, vector_num * sizeof(pair));
    print_err_with_assertg(pairs != NULL);
    for (i = 0; i < vector_num; i++)
    {
        pairs[i].index = i;
        pairs[i].eigenValue = eigenValues[i];
    }

    qsort(pairs, vector_num, sizeof(pair), compare_eigen);
    for (i = 0; i < vector_num; i++)
    {
        eigenValues[i] = pairs[i].eigenValue;
    }
}



int eigengapHeuristic()
{
    int i, maxIndex, k = 0;
    double maxGap = -1.0;
    double *eigenGaps;
    calcJacobi(Lnorm());
    sortpairs();
    eigenGaps = (double *)calloc(vector_num - 1, sizeof(double));
    for (i = 0; i < vector_num - 1; i++)
    {
        eigenGaps[i] = fabs(eigenValues[i] - eigenValues[i + 1]);
    }
    maxIndex = (int)floor(vector_num / 2);
    for (i = 0; i < maxIndex; i++)
    {
        if (eigenGaps[i] > maxGap)
        {
            maxGap = eigenGaps[i];
            k = i;
        }
    }
    free(eigenGaps);
    return k + 1;
}



void createTMat()
{ /*using vectors_mat*/
    int i, j;
    double sum;
    /* Transpose vectors_mat */
    get_t_mat(vectors_mat);
    UMat = (double **)calloc(vector_num, k * sizeof(double));
    print_err_with_assertg(UMat != NULL);
    for (i = 0; i < vector_num; i++)
    {
        UMat[i] = (double *)calloc(k, sizeof(double));
        print_err_with_assertg(UMat[i] != NULL);
        for (j = 0; j < k; j++)
        {
            UMat[i][j] = vectors_mat[(pairs[j]).index][i];
        }
    }
    /*normalizing UMat*/
    for (i = 0; i < vector_num; i++)
    {
        sum = 0;
        for (j = 0; j < k; j++)
        {
            sum += (UMat[i][j]) * (UMat[i][j]);
        }
        sum = sqrt(sum);
        if (sum != 0)
        {
            for (j = 0; j < k; j++)
            {
                UMat[i][j] = (UMat[i][j]) / sum;
            }
        }
    }
    TMat = UMat; /* TMat is a normalized Umat*/
    /*
    printf("here is Tmat\n");
   print_mat(TMat, vector_num, k);
    printf("\nhere is Tmat\n");
    */
}

void init_centroids(int k)
{
    int i, j, z;
    print_err_with_assertg(k < vector_num);
    centroids = (double **)calloc(k, vector_len * sizeof(double));
    print_err_with_assertg(centroids != NULL);
    for (i = 0; i < k; i++)
    {
        centroids[i] = (double *)calloc(vector_len, sizeof(double));
        print_err_with_assertg(centroids[i] != NULL);
    }
    for (z = 0; z < k; z++)
    {
        for (j = 0; j < vector_len; j++)
        {
            centroids[z][j] = vector_list[z][j];
        }
    }
    clusters = (double **)calloc(k, sizeof(double *));
    assert_double_mat(clusters);
}

double get_distance(double *vector_1, double *vector_2)
{
    double ret;
    int i;
    ret = 0;
    for (i = 0; i < vector_len; i++)
    {
        ret+=pow(vector_1[i]-vector_2[i], 2);
    }
    return ret;
}

int min_dist_centroid(double *v)
{
    double min;
    int ind, i;
    double temp;
    min= get_distance(v, centroids[0]);
    ind = 0;
    for (i = 0; i < k; i++)
    {
        temp = get_distance(v, centroids[i]);
        if (temp < min)
        {
            min = temp;
            ind = i;
        }
    }
    return ind;
}

void vector_to_cluster(int k)
{
    int *clusterssizes; /*for realloc*/
    int ind;
    int i;
    free(clustersindexes);
    clustersindexes = (int *)calloc(k, sizeof(int));
    print_err_with_assertg(clustersindexes != NULL);
    clusterssizes = (int *)calloc(k, sizeof(int));
    print_err_with_assertg(clusterssizes != NULL);
    for (i = 0; i < k; i++)
    { /*initialize each cluster's size to 100*/
        clusterssizes[i] = 100;
    }
    for (i = 0; i < k; i++)
    {
        free(clusters[i]);
        clusters[i] = (double *)calloc(100, sizeof(double));
        print_err_with_assertg(clusters[i] != NULL);
    }
    for (i = 0; i < vector_num; i++)
    {
        ind = min_dist_centroid(vector_list[i]);
        if (clustersindexes[ind] > ((clusterssizes[ind]) / 2))
        { /*Increase if necessary*/
            clusters[ind] = (double *)realloc(clusters[ind], 2 * clusterssizes[ind] * sizeof(double));
            clusterssizes[ind] *= 2;
        }
        clusters[ind][clustersindexes[ind]] = i;
        clustersindexes[ind]++; /*increase number of vectors in specified cluster*/
    }
    free(clusterssizes);
}

double *cluster_to_centroid(int index)
{
    int i;
    int j;
    int vector_index;
    int num = clustersindexes[index]; /* number of vectors in given cluster */
    double *res = (double *)calloc(vector_len, sizeof(double));
    print_err_with_assertg(res != NULL);
    if (num != 0)
    {
        for (i = 0; i < vector_len; i++)
        {
            for (j = 0; j < num; j++)
            {
                vector_index = (int)clusters[index][j]; /* not actual vector but index in vector_list */
                res[i] += vector_list[vector_index][i]; /*relevant cluster*/
            }
        }

        for (i = 0; i < vector_len; i++)
        {
            res[i] = res[i] / (num);
        }
    }
    else
    {
        for (i = 0; i < vector_len; i++)
        {
            res[i] = centroids[index][i];
        }
    }
    return res;
}

int areequal(double *arr1, double *arr2)
{
    int i;
    for (i = 0; i < vector_len; i++)
    {
        if (arr1[i] != arr2[i])
        {
            return 0;
        }
    }
    return 1;
}

int update_centroids()
{
    int is_changed = 0;
    int x, i, j;
    for (i = 0; i < k; i++)
    {
        double *newcentroid, *ret;
        newcentroid = (double *)calloc(vector_len, sizeof(double));
        print_err_with_assertg(newcentroid != NULL);
        ret = cluster_to_centroid(i);
        for (j = 0; j < vector_len; j++)
        {
            newcentroid[j] = ret[j];
        }
        if (areequal(centroids[i], newcentroid) == 0) {is_changed = 1;}
        for (x = 0; x < vector_len; x++)
            centroids[i][x] = newcentroid[x];
        free(newcentroid);
        free(ret);
    }
    return (is_changed != 0);
}

double **calccentroids(int max_iter)
{
    int counter, isequal;
    print_err_with_assertg(clusters_num == 0 || clusters_num > 0);
    counter = 0, isequal = 1;
    clusters = (double **)calloc(clusters_num, sizeof(double *)); /* originally in init_centroids */
    while (counter < max_iter && isequal == 1)
    {
        vector_to_cluster(clusters_num);
        isequal = update_centroids();
        counter++;
    }
    free_full_array_by_size(clusters, clusters_num);
    return centroids;
}
void fullSpectral()
{
    int counter, isequal, max_iter = 300; 
    counter = 0, isequal = 1;
    k = k == 0 ? eigengapHeuristic() : k;
    vector_len = k;
    createTMat();
    free_full_array_by_size(vector_list, vector_num);
    vector_list = UMat;
    init_centroids(k);
    while (counter < max_iter && isequal == 1)
    {
        vector_to_cluster(k);
        isequal = update_centroids();
        counter++;
    }
    
}