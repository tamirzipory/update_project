#define PY_SSIZE_T_CLEAN

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include "spkmeans.h"
#define eps pow(10, -5)



/* GOLOBAL ELEMENTS: */
typedef struct pair
{
    double eigenValue;
    int index;
} pair;

double **clusters;
double **centroids;
double **vector_list;
int *clustersindexes;
int vector_len;
int vector_num;
int k;
int clusters_num;
double **weighted_mat;
double **diag_degree_mat;
double **norm_mat;
pair *pairs;
double *eigenValues;
double **vectors_mat;
double **the_U_mat;
double **the_T_mat;
int max_iter;
float inputK;
char *goal;
double **python_centroids;

int eigenComperator(const void *a, const void *b){
    struct pair *A = (struct pair *)a;
    struct pair *B = (struct pair *)b;
    A = (pair *)a;
    B = (pair *)b;
    if (A->eigenValue == B->eigenValue)
        return A->index - B->index;
    else
       return A->eigenValue > B-> eigenValue ? 1 : -1;
    
}

double ** init_double_mat(int rows){
    int i;
    double** ret = (double**) calloc(rows, rows*sizeof(double));
    catch_err_of_int(ret != NULL);
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
    qsort(pairs, vector_num, sizeof(pair), eigenComperator);
    for (i = 0; i < vector_num; i++)
        eigenValues[i] = pairs[i].eigenValue;
    
}

int eigengapHeuristic(){
    int i, maxIndex, k;
    double temp_maximum;
    double *e_arr_gaps;
    temp_maximum = -1.0;
    k = 0;
    calc_jacobi_patterm(Lnorm());
    sortpairs();
    e_arr_gaps = (double *)calloc(vector_num - 1, sizeof(double));
    for (i = 0; i < vector_num - 1; i++)
        e_arr_gaps[i] = fabs(eigenValues[i] - eigenValues[i + 1]);
    maxIndex = (int)floor(vector_num / 2);
    for (i = 0; i < maxIndex; i++){
        if (e_arr_gaps[i] > temp_maximum){
            temp_maximum = e_arr_gaps[i];
            k = i;
        }
    }
    free(e_arr_gaps);
    return k + 1;
}



void create_T_mat(){ /*using vectors_mat*/
    int i, j;
    double sum;
    /* Transpose vectors_mat */
    transpMat(vectors_mat);
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
    printf("here is the_T_mat\n");
   print_mat(the_T_mat, vector_num, k);
    printf("\nhere is the_T_mat\n");
    */
}
/*Adar - to do until today please*/
void init_centroids(int k){
    int i, j, x;
    catch_err_of_int(k < vector_num);
    centroids = (double **)calloc(k, vector_len * sizeof(double));
    catch_err_of_int(centroids != NULL);
    for (i = 0; i < k; i++){
        centroids[i] = (double *)calloc(vector_len, sizeof(double));
        catch_err_of_int(centroids[i] != NULL);
    }
    for (x = 0; x < k; x++){
        for (j = 0; j < vector_len; j++)
            centroids[x][j] = vector_list[x][j];
    }
    clusters = (double **)calloc(k, sizeof(double *));
    assert_double_mat(clusters);
}

double get_distance(double *v1, double *v2){
    double ret;
    int i;
    ret = 0;
    for (i = 0; i < vector_len; i++)
        ret+=pow(v1[i]-v2[i], 2);
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
/*Tamir section*/
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
        clustersindexes[ind]++; 
    }
    free(clusterssizes);
}

double *cluster_to_centroid(int index){
    int i, j, num, vector_index;
    double *ret = (double *)calloc(vector_len, sizeof(double));
    num = clustersindexes[index]; /* number of vectors in given cluster */
    catch_err_of_int(ret != NULL);
    if (num != 0){
        for (i = 0; i < vector_len; i++){
            for (j = 0; j < num; j++){
                vector_index = (int)clusters[index][j]; 
                ret[i] += vector_list[vector_index][i]; 
            }
        }
        for (i = 0; i < vector_len; i++)
            ret[i] = ret[i] / (num);
    }
    else{
        for (i = 0; i < vector_len; i++)
            ret[i] = centroids[index][i];
    }
    return ret;
}

int areequal(double *arr1, double *arr2){
    int i;
    for (i = 0; i < vector_len; i++){
        if (arr1[i] != arr2[i])
            return 0; 
    }
    return 1;
}

int update_centroids(){
    int x, i, j, changed;
    changed = 0;
    for (i = 0; i < k; i++){
        double *newcentroid;
        double *res;
        newcentroid = (double *)calloc(vector_len, sizeof(double));
        catch_err_of_int(newcentroid != NULL);
        res = cluster_to_centroid(i);
        for (j = 0; j < vector_len; j++)
            newcentroid[j] = res[j];
        if (areequal(centroids[i], newcentroid) == 0)
            changed++;
        for (x = 0; x < vector_len; x++)
            centroids[i][x] = newcentroid[x];
        free(newcentroid);
        free(res);
    }
    return (changed != 0);
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
    freearray(clusters, clusters_num);
    return centroids;
}
void fullSpectral(){
    int count_var, isequal, max_iter;
    max_iter = 300; 
    count_var = 0;
    isequal = 1;
    k = k == 0 ? eigengapHeuristic() : k;
    vector_len = k;
    create_T_mat();
    freearray(vector_list, vector_num);
    vector_list = the_U_mat;
    init_centroids(k);
    while (count_var < max_iter && isequal == 1){
        vector_to_cluster(k);
        isequal = update_centroids();
        count_var++;
    }  
}

/*calculate the len of the vector*/

/*Adar section*/
int find_vector_num(FILE *fp){
    int ret;
    char ch, before;
    before = 0;
    ret = 0;
    while((ch = fgetc(fp))!=EOF){
        if(ch == '\n')
            ret++;
        before = ch;
    }
    if(before!='\n')
        ret++;
    vector_num = ret;
    return ret;  
}



/* calculate the number of lines*/
int find_vector_len(FILE *file){
    char ch;
    int the_len;
    the_len = 1;
    while ((ch = fgetc(file)) != '\n'){
        if (ch == ',')
            the_len++;
    }
    vector_len = the_len;
    return the_len;
}

void init_vectors_mat(int rows, int cols){
    int i;
    vector_list = (double **) calloc(rows, sizeof(double));
    catch_err_of_int(vector_list != NULL);
    assert_double_mat(vector_list);
    for(i = 0; i < rows; i++){
        vector_list[i] = (double *) calloc(cols, sizeof(double));
        catch_err_of_int(vector_list[i] != NULL);
    }
}

void catch_err_of_int(int x){
    if (!x){
        printf("An Error Has Occured\n");
        assert(x);
    }
}

/*read the file*/
void read_file(FILE *file){
    char *set_split;
    int i, j;
    char line[1000];
    find_vector_len(file); 
    rewind(file);
    find_vector_num(file); 
    rewind(file);
    init_vectors_mat(vector_num, vector_len);
    j = 0;
    while (fgets(line, 1000, file) != NULL){ 
        set_split = strtok(line, ","); 
        for (i = 0; i < vector_len; i++){
            vector_list[j][i] = atof(set_split); 
            set_split = strtok(NULL, ",");
        }
        j++;
    }
    fclose(file);
}


/*kefel matritsot*/
double **matMult(double **A, double **B){
    int j, k, l;
    double **res = init_double_mat(vector_num);
    for (j = 0; j < vector_num; j++){
        for (k = 0; k < vector_num; k++){
            res[j][k] = 0;
            for (l = 0; l < vector_num; l++)
                res[j][k] += A[j][l] * B[l][k];
        }
    }
    return res;
}

/* ret wam mat*/
double **weightedAdjMat(){
    int j, k;
    double *v1;
    weighted_mat = init_double_mat(vector_num);
    for (j = 0; j < vector_num; j++){
        v1 = vector_list[j];
        for (k = j + 1; k < vector_num; k++){
            double *v2 = vector_list[k];
            weighted_mat[j][k] = weightedAdjMat_calc(v1, v2);
            weighted_mat[k][j] = weighted_mat[j][k]; /*symmetry*/
        }
    }
    return weighted_mat;
}
/*cat wam_agj*/
double weightedAdjMat_calc(double *v1, double *v2){
    double ret, dist;
    int i;
    dist = 0;
    for (i = 0; i < vector_len; i++)
        dist += (v1[i] - v2[i]) * (v1[i] - v2[i]);
    ret = sqrt(dist);
    ret = -((ret) / 2);
    ret = exp(ret);
    return ret;
}

double **Lnorm(){
    int i, j, k;
    double **res;
    weightedAdjMat();
    diagDegMat();
    for (i = 0; i < vector_num; i++)
     /* D^(-1/2) */
        diag_degree_mat[i][i] = 1 / (sqrt(diag_degree_mat[i][i]));
    
    res = matMult(diag_degree_mat, weighted_mat);
    norm_mat = matMult(res, diag_degree_mat);
    for (j = 0; j < vector_num; j++){
        for (k = 0; k < vector_num; k++)
            norm_mat[j][k] = j == k ? 1 - norm_mat[j][k]: - (norm_mat[j][k]);
    }
    freearray(res, vector_num);
    return norm_mat;
}

/* print the matix by num of lines and rows */
void print_mat(double **mat_to_print, int num_of_row, int num_of_col){
    int i, j;
    for (i = 0; i < num_of_row; i++){
        for (j = 0; j < num_of_col; j++){

            if ((mat_to_print[i][j] < 0) && (mat_to_print[i][j] > -0.00005))
                mat_to_print[i][j] = 0;

            if (j == num_of_col - 1){
                printf("%0.4f", mat_to_print[i][j]);
                if (i < num_of_row - 1)
                    printf("\n");  
            }
            else
                printf("%0.4f,", mat_to_print[i][j]);   
        }
    }
}

double **diagDegMat(){
    int i, j, k;
    double sum;
    diag_degree_mat = (double **)calloc(vector_num, vector_num * sizeof(double));
    catch_err_of_int(diag_degree_mat != NULL);
    for (i = 0; i < vector_num; i++){
        diag_degree_mat[i] = (double *)calloc(vector_num, sizeof(double));
        catch_err_of_int(diag_degree_mat[i] != NULL);
    }
    for (j = 0; j < vector_num; j++){
        sum = 0;
        for (k = 0; k < vector_num; k++)
            sum += weighted_mat[j][k];
        diag_degree_mat[j][j] = sum;
    }
    return diag_degree_mat;
}

void calcRotationMat(double **P, double c, double s, int row, int column){
    int i, j;
    for (i = 0; i < vector_num; i++){
        for (j = 0; j < vector_num; j++)
        /*id mat*/
            P[i][j] = i == j ? 1 : 0;
    }
    P[row][column] = s;
    P[column][row] = -s;
    P[row][row] = c;
    P[column][column] = c;
}

double calc_theta(double **Mat, int i, int j){
    double res, up, down;
    up = Mat[j][j] - Mat[i][i];
    down = 2 * Mat[i][j];
    res = up / down;
    return res;
}

double get_t_by_theta(double theta){
    int sign;
    double down;
    sign = theta>=0?1:-1;
    down = fabs(theta) + sqrt(theta * theta + 1);
    return sign / down;
}

double get_c(double t){
    double ret, down;
    down = sqrt(t * t + 1);
    ret = 1 / down;
    return ret;
}

double get_s(double t, double c) { return t * c; }

int isConverged(double **A, double **A_til_mat){
    int i, j, k, l;
    double ret, sumA, sum_A_tag, epsilon;
    epsilon = pow(10, -5);
    sumA = 0;
    sum_A_tag = 0;
    for (i = 0; i < vector_num; i++){
        for (j = 0; j < vector_num; j++){
            if (i != j)
                sumA += pow(A[i][j], 2);
        }
    }
    for (k = 0; k < vector_num; k++){
        for (l = 0; l < vector_num; l++){
            if (k != l)
                sum_A_tag += (A_til_mat[k][l]) * (A_til_mat[k][l]);   
        }
    }
    ret = sumA - sum_A_tag;

    return ret <= epsilon ? 1 : 0;
}

void calc_A_mat(double **A, double **A_til_mat, int i, int j, double c, double s){
    int k;
    for (k = 0; k < vector_num; k++){
        if ((k != i) && (k != j)){
            A_til_mat[k][i] = c * A[k][i] - s * A[k][j];
            A_til_mat[i][k] = A_til_mat[k][i];
            A_til_mat[k][j] = c * A[k][j] + s * A[k][i];
            A_til_mat[j][k] = A_til_mat[k][j];
        }
    }
    A_til_mat[i][i] = pow(c, 2) * A[i][i] + pow(s, 2) * A[j][j] - 2 * s * c * A[i][j];
    A_til_mat[j][j] = pow(s, 2) * A[i][i] + pow(c, 2) * A[j][j] + 2 * s * c * A[i][j];
    A_til_mat[i][j] = 0;
    A_til_mat[j][i] = A_til_mat[i][j];
}

double **calc_jacobi_patterm(double **A){
    double **A_til_mat;
    double **prime_mat;
    double **tmp;
    int i, j, l, row, col, is_con, count_var;
    double c, t, s, theta;
    row = 0, col = 1, is_con = 0, count_var = 0;
    A_til_mat = init_double_mat(vector_num);
    prime_mat = init_double_mat(vector_num);

   vectors_mat = init_double_mat(vector_num);
   
    for (l = 0; l < vector_num; l++)
    /* init to id mat*/
        vectors_mat[l][l] = 1;
    

    for (i = 0; i < vector_num; i++){ /*copy A*/
        for (j = 0; j < vector_num; j++)
            A_til_mat[i][j] = A[i][j];
        
    }
    while ((is_con == 0) && (count_var < 100)){
        for (i = 0; i < vector_num; i++){
            for (j = i + 1; j < vector_num; j++){
                if (fabs(A[i][j]) > fabs(A[row][col])){
                    row = i;
                    col = j;
                }
            }
        }

        if (A[row][col] == 0)
            break;
        
        theta = calc_theta(A, row, col);
        t = get_t_by_theta(theta);
        c = get_c(t);
        s = get_s(t, c);
        s = t*c;

        calcRotationMat(prime_mat, c, s, row, col);
        tmp = vectors_mat;
        vectors_mat = matMult(vectors_mat, prime_mat);
        freearray(tmp, vector_num);
        calc_A_mat(A, A_til_mat, row, col, c, s);
        is_con = isConverged(A, A_til_mat);

        for (i = 0; i < vector_num; i++) { /*copy A_til_mat to A*/
            for (j = 0; j < vector_num; j++)
                A[i][j] = A_til_mat[i][j];
        }
        count_var++;
    }
    eigenValues = (double *)calloc(vector_num, sizeof(double));
    catch_err_of_int(eigenValues != NULL);
    for (i = 0; i < vector_num; i++)
        eigenValues[i] = A_til_mat[i][i];
    freearray(A_til_mat, vector_num);
    freearray(prime_mat, vector_num);
    return A;
}

double **transpMat(double **Mat){
    int i, j;
    double temp;
    for (i = 1; i < vector_num; i++){
        for (j = 0; j < i; j++){
            temp = vectors_mat[i][j];
            Mat[i][j] = Mat[j][i];
            Mat[j][i] = temp;
        }
    }
    return Mat;
}

void freearray(double **array, int length){
    int i;
    for (i = 0; i < length; i++)
        free(array[i]);
    free(array);
}

void assert_double_mat(double **mat){
    if (mat == NULL){
        printf("An Error Has Occured");
        exit(0);
    }
}

void assert_double_arr(const double *arr){
    if (arr == NULL){
        printf("An Error Has Occured");
        exit(0);
    }
}

double **create_I_mat(int dim){
    /*just create the I matrix in dimentions of (dim x dim) */
    int i, j;
    double **id_mat;
    double *temp;

    temp = calloc(dim * dim, sizeof(double));
    assert_double_arr(temp);
    id_mat = calloc(dim, sizeof(double *));
    assert_double_mat(id_mat);
    for (i = 0; i < dim; i++){
        id_mat[i] = temp + i * dim;
        for (j = 0; j < dim; j++){
            if (i == j)
                id_mat[i][j] = 1.0;
            /*  
            not neccecary  
            else
                id_mat[i][j] = 0.0;
                */
        }
    }
    return id_mat;
}


double off(double **A, int dim){
    /*just calcuate A^2*/
    int i, j;
    double res;
    res = 0;
    for (i = 0; i < dim; i++){
        for (j = 0; j < dim; j++){
            if (i != j)
                res += pow(A[i][j], 2);    
        }
    }
    return res;
}

int sign(double num){return num < 0 ? -1 : 1;}

void V_kaful_P(double **mat, double s, double c, int num_of_iter, int i, int j){
   
    int r;
    double vector_cor_i, vector_cor_j;
    for (r = 0; r < num_of_iter; r++){
        vector_cor_i = mat[r][i];
        vector_cor_j = mat[r][j];
        mat[r][i] = (c * vector_cor_i) - (s * vector_cor_j);
        mat[r][j] = (s * vector_cor_i) + (c * vector_cor_j);
    }
}

double *get_ith_column(double **mat, int col_ind, int num_of_iter){
    int i;
    double *col;
    col = calloc(num_of_iter, sizeof(double));
    assert_double_arr(col);
    for (i = 0; i < num_of_iter; i++)
        col[i] = mat[i][col_ind];
    return col;
}

void assert_int_arr(const int *arr){
   
    if (arr == NULL){
        printf("An Error Has Occured");
        exit(0);
    }
}

int *max_indices_off_diag(double **A, int num_of_iter){
    
    double val;
    int i, j, max_i, max_j;
    int *arr;
    val = -1;
    max_i = 0, max_j = 0;
    arr = calloc(2, sizeof(int));
    assert_int_arr(arr);

    for (i = 0; i < num_of_iter; i++){
        for (j = 0; j < num_of_iter; j++){
            if (i != j){
                if (fabs(A[i][j]) > val){
                    val = fabs(A[i][j]);
                    max_i = i;
                    max_j = j;
                }
            }
        }
    }
    arr[0] = max_i;
    arr[1] = max_j;
    return arr;
}

#include "jacobi.c"

void print_double(double num){
    /*I don't print -0.0000 , I print 0.0000*/
    if ((num < 0) && (num > -0.00005))
        printf("0.0000");
    else
        printf("%.4f", num);
}

void print_row(double *row, int len){
    /*Prints thw array (the row is the array) according to given dimension*/
    int i;
    for (i = 0; i < len; i++){
        print_double(row[i]);
        if (i != len - 1)
            printf(",");
    }
}
double *get_diag(double **mat, int N){
    
    int i;
    double *diag;
    diag = calloc(N, sizeof(double));
    assert_double_arr(diag);

    for (i = 0; i < N; i++)
        diag[i] = mat[i][i];
    return diag;
}

void free_mat(double **mat){
    free(mat[0]);
    free(mat);
}


int main(int argc, char *argv[]){
    FILE *file;
    if (!(argc == 3)){
        printf("Invalid Input!\n");
        assert(argc == 3);
        exit(0);
    }
    goal = argv[1];
    file = fopen(argv[2], "r");
    if(file == NULL){
        printf("Invalid Input!\n");
        assert(file!=NULL);
        exit(0);
    }
    read_file(file);
    if (strcmp(goal, "wam") == 0){
        print_mat(weightedAdjMat(), vector_num, vector_num);
        freearray(weighted_mat, vector_num);
        freearray(vector_list, vector_num);
    }
    /* ddg */
    else if (strcmp(goal, "ddg") == 0){
        weightedAdjMat();
        print_mat(diagDegMat(), vector_num, vector_num);
        freearray(weighted_mat, vector_num);
        freearray(diag_degree_mat, vector_num);
        freearray(vector_list, vector_num);
    }
    /*lnorm*/
    else if (strcmp(goal, "lnorm") == 0){
        print_mat(Lnorm(), vector_num, vector_num);
        freearray(weighted_mat, vector_num);
        freearray(diag_degree_mat, vector_num);
        freearray(norm_mat, vector_num);
        freearray(vector_list, vector_num);
    }
    /*jacobi*/
    else if (strcmp(goal, "jacobi") == 0){
        int i;
        double **first_mat;
        first_mat = calc_jacobi_patterm(vector_list);
        for (i = 0; i < vector_num; i++){
            if ((first_mat[i][i] < 0) && (first_mat[i][i] > -0.00005))
                first_mat[i][i] = 0;
            if (i == vector_num - 1)
                printf("%0.4f \n", first_mat[i][i]);
            else
                printf("%0.4f,", first_mat[i][i]);
        }

        print_mat(vectors_mat, vector_num, vector_num);
        freearray(vectors_mat, vector_num);
        free(eigenValues);
        freearray(vector_list, vector_num);  
    }
    else
        catch_err_of_int(0 != 0);
    printf("\n");
    return 0;
}
