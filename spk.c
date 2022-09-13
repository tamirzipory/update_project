#define eps pow(10, -5)
#include "spkmeans.h"

/*
Musn't to compile it! it's only fot the organization of the code!
*/


/*the comparator for the qsort*/
int compare_eigen(const void *a, const void *b){
    struct pair *A, *B;
    A = (pair *)a;
    B = (pair *)b;
    return A->eigenValue == B->eigenValue ? A->index-B->index : A->eigenValue > B-> eigenValue ? 1 : -1;
}

double **calloc_mat(int rows){
    double** ret = (double**) calloc(rows, rows*sizeof(double));
    catch_err_of_int(ret != NULL);
    return ret;
}

double ** init_double_mat(int rows){
    int i;
    double **ret;
    ret = calloc_mat(rows);
    for(i = 0; i < rows; i++){
        ret[i] = (double *) calloc(rows, sizeof(double));
        catch_err_of_int(ret[i] != NULL);
    }
    return ret;
}

void sortpairs(){
    int i;
    pairs = (pair *)calloc(vector_num, vector_num * sizeof(pair));
    catch_err_of_int(pairs != NULL);
    for (i = 0; i < vector_num; i++){
        pairs[i].index = i;
        pairs[i].eigenValue = eigenValues[i];
    }

    qsort(pairs, vector_num, sizeof(pair), compare_eigen);
    for (i = 0; i < vector_num; i++)
        eigenValues[i] = pairs[i].eigenValue;  
}



int eigengapHeuristic(){
    int i, maxIndex, k = 0;
    double maxGap = -1.0;
    double *eigenGaps;
    calc_jacobi_patterm(Lnorm());
    sortpairs();
    eigenGaps = (double *)calloc(vector_num - 1, sizeof(double));
    for (i = 0; i < vector_num - 1; i++)
        eigenGaps[i] = fabs(eigenValues[i] - eigenValues[i + 1]);
    maxIndex = (int)floor(vector_num / 2);
    for (i = 0; i < maxIndex; i++){
        if (eigenGaps[i] > maxGap){
            maxGap = eigenGaps[i];
            k = i;
        }
    }
    free(eigenGaps);
    return k + 1;
}

void create_T_mat(){ /*using vectors_mat*/
    int i, j;
    double sum;
    /* Transpose vectors_mat */
    get_t_mat(vectors_mat);
    the_U_mat = (double **)calloc(vector_num, k * sizeof(double));
    catch_err_of_int(the_U_mat != NULL);
    for (i = 0; i < vector_num; i++){
        the_U_mat[i] = (double *)calloc(k, sizeof(double));
        catch_err_of_int(the_U_mat[i] != NULL);
        for (j = 0; j < k; j++)
            the_U_mat[i][j] = vectors_mat[(pairs[j]).index][i];
        
    }
    /*normalizing the_U_mat*/
    for (i = 0; i < vector_num; i++){
        sum = 0;
        for (j = 0; j < k; j++)
            sum += (the_U_mat[i][j]) * (the_U_mat[i][j]);
        
        sum = sqrt(sum);
        if (sum != 0){
            for (j = 0; j < k; j++)
                the_U_mat[i][j] = (the_U_mat[i][j]) / sum; 
        }
    }
    the_T_mat = the_U_mat; /* the_T_mat is a normalized the_U_mat*/

    /*
    to check the T-
    printf("here is the_T_mat\n");
   print_mat(the_T_mat, vector_num, k);
    printf("\nhere is the_T_mat\n");
    */
}

void init_centroids(int k){
    int i, j, z;
    catch_err_of_int(k < vector_num);
    centroids = (double **)calloc(k, vector_len * sizeof(double));
    catch_err_of_int(centroids != NULL);
    for (i = 0; i < k; i++){
        centroids[i] = (double *)calloc(vector_len, sizeof(double));
        catch_err_of_int(centroids[i] != NULL);
    }
    for (z = 0; z < k; z++){
        for (j = 0; j < vector_len; j++)
            centroids[z][j] = vector_list[z][j];
        
    }
    clusters = (double **)calloc(k, sizeof(double *));
    assert_double_mat(clusters);
}

double get_distance(double *vector_1, double *vector_2){
    double ret;
    int i;
    ret = 0;
    for (i = 0; i < vector_len; i++)
        ret+=pow(vector_1[i]-vector_2[i], 2);
    
    return ret;
}

int min_dist_centroid(double *v){
    double min, temp;
    int ind, i;   
    min= get_distance(v, centroids[0]);
    ind = 0;
    for (i = 0; i < k; i++){
        temp = get_distance(v, centroids[i]);
        if (temp < min){
            min = temp;
            ind = i;
        }
    }
    return ind;
}

void vector_to_cluster(int k){
    int *clusterssizes; 
    int ind, i;
    free(clustersindexes);
    clustersindexes = (int *)calloc(k, sizeof(int));
    catch_err_of_int(clustersindexes != NULL);
    clusterssizes = (int *)calloc(k, sizeof(int));
    catch_err_of_int(clusterssizes != NULL);
    for (i = 0; i < k; i++)
        clusterssizes[i] = 100;
    for (i = 0; i < k; i++){
        free(clusters[i]);
        clusters[i] = (double *)calloc(100, sizeof(double));
        catch_err_of_int(clusters[i] != NULL);
    }
    for (i = 0; i < vector_num; i++){
        ind = min_dist_centroid(vector_list[i]);
        if (clustersindexes[ind] > ((clusterssizes[ind]) / 2)){
            clusters[ind] = (double *)realloc(clusters[ind], 2 * clusterssizes[ind] * sizeof(double));
            clusterssizes[ind] *= 2;
        }
        clusters[ind][clustersindexes[ind]] = i;
        clustersindexes[ind]++; /*increase number of vectors in specified cluster*/
    }
    free(clusterssizes);
}

double *cluster_to_centroid(int index){
    int i, j, vector_index, num;
    num = clustersindexes[index]; /* number of vectors in given cluster */
    double *res = (double *)calloc(vector_len, sizeof(double));
    catch_err_of_int(res != NULL);
    if (num != 0){
        for (i = 0; i < vector_len; i++){
            for (j = 0; j < num; j++){
                vector_index = (int)clusters[index][j]; /* not actual vector but index in vector_list */
                res[i] += vector_list[vector_index][i]; /*relevant cluster*/
            }
        }

        for (i = 0; i < vector_len; i++)
            res[i] = res[i] / (num);
        
    }
    else{
        for (i = 0; i < vector_len; i++)
            res[i] = centroids[index][i];   
    }
    return res;
}

int areequal(double *arr1, double *arr2){
    int i, ret;
    for (i = 0; i < vector_len; i++){
        if (arr1[i] != arr2[i])
            return 0;
    }
    ret = 1;
    return ret;
}

int update_centroids(){
    int x, i, j, is_changed;
    is_changed = 0;
    for (i = 0; i < k; i++){
        double *newcentroid, *ret;
        newcentroid = (double *)calloc(vector_len, sizeof(double));
        catch_err_of_int(newcentroid != NULL);
        ret = cluster_to_centroid(i);
        for (j = 0; j < vector_len; j++)
            newcentroid[j] = ret[j];
        if (areequal(centroids[i], newcentroid) == 0) {is_changed = 1;}
        for (x = 0; x < vector_len; x++)
            centroids[i][x] = newcentroid[x];
        free(newcentroid);
        free(ret);
    }
    return (is_changed != 0);
}

double **calccentroids(int max_iter){
    int count_var, isequal;
    catch_err_of_int(clusters_num == 0 || clusters_num > 0);
    count_var = 0, isequal = 1;
    clusters = (double **)calloc(clusters_num, sizeof(double *)); /* originally in init_centroids */
    while (count_var < max_iter && isequal == 1){
        vector_to_cluster(clusters_num);
        isequal = update_centroids();
        count_var++;
    }
    free_full_array_by_size(clusters, clusters_num);
    return centroids;
}
void fullSpectral(){
    int count_var, isequal, max_iter = 300; 
    count_var = 0, isequal = 1;
    k = k == 0 ? eigengapHeuristic() : k;
    vector_len = k;
    createthe_T_mat();
    free_full_array_by_size(vector_list, vector_num);
    vector_list = the_U_mat;
    init_centroids(k);
    while (count_var < max_iter && isequal == 1){
        vector_to_cluster(k);
        isequal = update_centroids();
        count_var++;
    }
    
}
