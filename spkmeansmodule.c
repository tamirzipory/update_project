#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <math.h>
#include "spkmeans.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

/*I want to say that in this section (connect python and C) 
we got help from our friends that show us exactly how to do it.
We try to do it alone during 3 weeks but we don't success.
*/

/*duplicate the struct.. it's not work us without it*/
typedef struct pair
{
    double eigenValue;
    int index;
} pair;

/* Global variables */
int clusters_num; /* make sure clusters_num == K*/
int k;
int vector_num;
int vector_len;
int max_iter;
double **vector_list;
double **weighted_mat;
double **diag_degree_mat;
double **norm_mat;
double **vectors_mat;
double **the_U_mat;
double **the_T_mat;
double *eigenValues;
int vector_num;
int vector_len;
pair *pairs;
double **centroids;
double **clusters;
int *clustersindexes;

static PyObject *fit(PyObject *self, PyObject *args){
    PyObject *centroids_list;
    PyObject *origin_vector_list;
    Py_ssize_t m, n, i, j;
    int max_iter;

    if (!PyArg_ParseTuple(args, "O!O!ii", &PyList_Type, &centroids_list, &PyList_Type, &origin_vector_list, &max_iter, &clusters_num))
        return NULL;
    n = PyList_Size(PyList_GetItem(centroids_list, 0)); /* vector_len*/
    vector_len = (int)n;

    centroids = (double **)calloc(clusters_num, n * sizeof(double));
    catch_err_of_int(centroids != NULL);
    for (i = 0; i < clusters_num; i++){
        centroids[i] = (double *)calloc(n, sizeof(double));
        catch_err_of_int(centroids[i] != NULL);
    }

    for (i = 0; i < clusters_num; i++){
        for (j = 0; j < n; j++)
            centroids[i][j] = PyFloat_AsDouble(PyList_GetItem(PyList_GetItem(centroids_list, i), j)); /*CONVERSION*/
    }

    m = PyList_Size(origin_vector_list); /* vector_num */
    vector_num = (int)m;

    vector_list = (double **)calloc(m, n * sizeof(double));
    catch_err_of_int(vector_list != NULL);
    for (i = 0; i < m; i++){
        vector_list[i] = (double *)calloc(n, sizeof(double));
        catch_err_of_int(vector_list[i] != NULL);
    }

    for (i = 0; i < m; i++){
        for (j = 0; j < n; j++)
            vector_list[i][j] = PyFloat_AsDouble(PyList_GetItem(PyList_GetItem(origin_vector_list, i), j)); 
    } 
    calccentroids(max_iter);

    PyObject *output_centroids = PyList_New(0);
    for (i = 0; i < clusters_num; i++){
        PyObject *centroid = PyList_New(0);
        for (j = 0; j < n; j++)
            PyList_Append(centroid, Py_BuildValue("d", centroids[i][j]));
        PyList_Append(output_centroids, centroid);
    }
    freearray(vector_list, vector_num);
    freearray(centroids, k);
    return output_centroids;
}

static PyObject *WeightedAdjacencyMatrix(PyObject *self, PyObject *args){
    PyObject *origin_vector_list;
    Py_ssize_t m, n, i, j;

    if (!PyArg_ParseTuple(args, "O!ii", &PyList_Type, &origin_vector_list, &vector_num, &vector_len))
        return NULL;

    n = PyList_Size(PyList_GetItem(origin_vector_list, 0));
    m = PyList_Size(origin_vector_list);
    vector_list = (double **)calloc(m, n * sizeof(double));
    catch_err_of_int(vector_list != NULL);
    for (i = 0; i < m; i++){
        vector_list[i] = (double *)calloc(n, sizeof(double));
        catch_err_of_int(vector_list[i] != NULL);
    }

    for (i = 0; i < m; i++){
        for (j = 0; j < n; j++)
            vector_list[i][j] = PyFloat_AsDouble(PyList_GetItem(PyList_GetItem(origin_vector_list, i), j)); /*CONVERSION*/
    } 
    weightedAdjMat();

    PyObject *output_matrix = PyList_New(0);
    for (i = 0; i < vector_num; i++){
        PyObject *centroid = PyList_New(0);
        for (j = 0; j < vector_num; j++)
            PyList_Append(centroid, Py_BuildValue("d", weighted_mat[i][j]));
        PyList_Append(output_matrix, centroid);
    }
    freearray(vector_list, vector_num);
    freearray(weighted_mat, vector_num);

    return output_matrix;
}

static PyObject *DiagonalDegreeMatrix(PyObject *self, PyObject *args){
    PyObject *origin_matrix;
    Py_ssize_t m, n, i, j;
    if (!PyArg_ParseTuple(args, "O!ii", &PyList_Type, &origin_matrix, &vector_num, &vector_len))
        return NULL;
    n = PyList_Size(PyList_GetItem(origin_matrix, 0));
    m = PyList_Size(origin_matrix);
    weighted_mat = (double **)calloc(m, n * sizeof(double));
    catch_err_of_int(weighted_mat != NULL);
    for (i = 0; i < m; i++){
        weighted_mat[i] = (double *)calloc(n, sizeof(double));
        assert(weightedAdjMarixt[i] != NULL); /* notice assert not catch_err_of_int because of compiler */
    }

    for (i = 0; i < m; i++){
        for (j = 0; j < n; j++)
            weighted_mat[i][j] = PyFloat_AsDouble(PyList_GetItem(PyList_GetItem(origin_matrix, i), j)); /*CONVERSION*/
    } 
    diagDegMat();

    PyObject *output_matrix = PyList_New(0);
    for (i = 0; i < vector_num; i++){
        PyObject *centroid = PyList_New(0);
        for (j = 0; j < vector_num; j++)
            PyList_Append(centroid, Py_BuildValue("d", diag_degree_mat[i][j]));
        PyList_Append(output_matrix, centroid);
    }
    freearray(weighted_mat, vector_num);
    freearray(diag_degree_mat, vector_num);

    return output_matrix;
}

static PyObject *NormalizedGraphLaplacian(PyObject *self, PyObject *args){
    PyObject *origin_vector_list;
    Py_ssize_t m, n, i, j;
    if (!PyArg_ParseTuple(args, "O!ii", &PyList_Type, &origin_vector_list, &vector_num, &vector_len))
        return NULL;
    n = PyList_Size(PyList_GetItem(origin_vector_list, 0));
    m = PyList_Size(origin_vector_list);
    vector_list = (double **)calloc(m, n * sizeof(double));
    catch_err_of_int(vector_list != NULL);
    for (i = 0; i < m; i++){
        vector_list[i] = (double *)calloc(n, sizeof(double));
        catch_err_of_int(vector_list[i] != NULL);
    }

    for (i = 0; i < m; i++){
        for (j = 0; j < n; j++)
            vector_list[i][j] = PyFloat_AsDouble(PyList_GetItem(PyList_GetItem(origin_vector_list, i), j)); 
    }
    Lnorm();

    PyObject *output_matrix = PyList_New(0);
    for (i = 0; i < vector_num; i++){
        PyObject *centroid = PyList_New(0);
        for (j = 0; j < vector_num; j++)
            PyList_Append(centroid, Py_BuildValue("d", norm_mat[i][j]));
        

        PyList_Append(output_matrix, centroid);
    }
    freearray(vector_list, vector_num);
    freearray(weighted_mat, vector_num);
    freearray(diag_degree_mat, vector_num);
    freearray(norm_mat, vector_num);
    return output_matrix;
}

static PyObject *Jacobi(PyObject *self, PyObject *args){
    PyObject *origin_vector_list;
    Py_ssize_t m, n, i, j;
    double **jacobi_mat;
    if (!PyArg_ParseTuple(args, "O!ii", &PyList_Type, &origin_vector_list, &vector_num, &vector_len))
        return NULL;

    n = PyList_Size(PyList_GetItem(origin_vector_list, 0));
    m = PyList_Size(origin_vector_list);
    jacobi_mat = (double **)calloc(m, n * sizeof(double));
    catch_err_of_int(jacobi_mat != NULL);
    for (i = 0; i < m; i++){
        jacobi_mat[i] = (double *)calloc(n, sizeof(double));
        catch_err_of_int(jacobi_mat[i] != NULL);
    }

    for (i = 0; i < m; i++){
        for (j = 0; j < n; j++)
            jacobi_mat[i][j] = PyFloat_AsDouble(PyList_GetItem(PyList_GetItem(origin_vector_list, i), j)); /*CONVERSION*/
        
    }
    calc_jacobi_patterm(jacobi_mat);

    PyObject *output_values = PyList_New(0);

    for (j = 0; j < vector_num; j++)
        PyList_Append(output_values, Py_BuildValue("d", eigenValues[j]));
    

    transpMat(vectors_mat);
    PyObject *output_vectors = PyList_New(0);
    for (i = 0; i < vector_num; i++){
        PyObject *centroid = PyList_New(0);
        for (j = 0; j < vector_num; j++)
            PyList_Append(centroid, Py_BuildValue("d", vectors_mat[i][j]));
        PyList_Append(output_vectors, centroid);
    }

    freearray(jacobi_mat, (int)m);
    freearray(vectors_mat, vector_num);
    free(eigenValues);

    return Py_BuildValue("OO", output_values, output_vectors);
}

static PyObject *heuristic(PyObject *self, PyObject *args){
    int newK;
    Py_ssize_t m, n, i, j;
    PyObject *origin_vector_list;
    if (!PyArg_ParseTuple(args, "O!ii", &PyList_Type, &origin_vector_list, &vector_num, &vector_len))
        return NULL;

    n = PyList_Size(PyList_GetItem(origin_vector_list, 0));
    m = PyList_Size(origin_vector_list);
    vector_list = (double **)calloc(m, n * sizeof(double));
    catch_err_of_int(vector_list != NULL);
    for (i = 0; i < m; i++){
        vector_list[i] = (double *)calloc(n, sizeof(double));
        catch_err_of_int(vector_list[i] != NULL);
    }

    for (i = 0; i < m; i++){
        for (j = 0; j < n; j++)
            vector_list[i][j] = PyFloat_AsDouble(PyList_GetItem(PyList_GetItem(origin_vector_list, i), j)); /*CONVERSION*/
    }

    newK = eigengapHeuristic();
    freearray(vector_list, (int)m);
    freearray(weighted_mat, vector_num);
    freearray(diag_degree_mat, vector_num);
    freearray(norm_mat, vector_num);
    free(pairs);
    return Py_BuildValue("i", newK); 
}

static PyObject *fullSpectralPy(PyObject *self, PyObject *args){
    PyObject *origin_vector_list;
    Py_ssize_t m, n, i, j;
    if (!PyArg_ParseTuple(args, "iiO!", &k, &vector_num, &PyList_Type, &origin_vector_list))
        return NULL;
    clusters_num = k;
    n = PyList_Size(PyList_GetItem(origin_vector_list, 0));
    m = PyList_Size(origin_vector_list);
    vector_len = (int)n;
    vector_list = (double **)calloc(m, n * sizeof(double));
    catch_err_of_int(vector_list != NULL);
    for (i = 0; i < m; i++){
        vector_list[i] = (double *)calloc(n, sizeof(double));
        catch_err_of_int(vector_list[i] != NULL);
    }

    for (i = 0; i < m; i++){
        for (j = 0; j < n; j++)
            vector_list[i][j] = PyFloat_AsDouble(PyList_GetItem(PyList_GetItem(origin_vector_list, i), j)); /*CONVERSION*/
    }

    eigengapHeuristic();
    create_T_mat();
    freearray(vector_list, vector_num);
    vector_list = the_U_mat;

    PyObject *output_matrix = PyList_New(0);
    for (i = 0; i < vector_num; i++){
        PyObject *centroid = PyList_New(0);
        for (j = 0; j < k; j++)
            PyList_Append(centroid, Py_BuildValue("d", vector_list[i][j]));
        PyList_Append(output_matrix, centroid);
    }
    freearray(vector_list, vector_num);
    freearray(weighted_mat, vector_num);
    freearray(diag_degree_mat, vector_num);
    freearray(norm_mat, vector_num);
    freearray(vectors_mat, vector_num);
    free(eigenValues);
    free(pairs);

    return output_matrix;
}

static PyObject *kmeans(PyObject *self, PyObject *args){
    PyObject *origin_vector_list, *origin_final_centroids;
    Py_ssize_t m, n, i, j;
    int count_var, isequal;
    if (!PyArg_ParseTuple(args, "iiiO!O!", &clusters_num, &vector_num, &vector_len, &PyList_Type, &origin_vector_list, &PyList_Type, &origin_final_centroids))    
        return NULL;
    k = clusters_num;
    n = PyList_Size(PyList_GetItem(origin_vector_list, 0));
    m = PyList_Size(origin_vector_list);
    vector_len = (int)n;
    vector_list = (double **)calloc(m, n * sizeof(double));
    catch_err_of_int(vector_list != NULL);
    for (i = 0; i < m; i++){
        vector_list[i] = (double *)calloc(n, sizeof(double));
        catch_err_of_int(vector_list[i] != NULL);
    }
    for (i = 0; i < m; i++){
        for (j = 0; j < n; j++){
            vector_list[i][j] = PyFloat_AsDouble(PyList_GetItem(PyList_GetItem(origin_vector_list, i), j)); /*CONVERSION*/
        }
    }
    n = PyList_Size(PyList_GetItem(origin_final_centroids, 0));
    m = PyList_Size(origin_final_centroids);
    centroids = (double **)calloc(k, vector_len * sizeof(double));
    catch_err_of_int(centroids != NULL);
    for (i = 0; i < k; i++){
        centroids[i] = (double *)calloc(vector_len, sizeof(double));
        catch_err_of_int(centroids[i] != NULL);
    }
    for (i = 0; i < k; i++){
        for (j = 0; j < vector_len; j++)
            centroids[i][j] = PyFloat_AsDouble(PyList_GetItem(PyList_GetItem(origin_final_centroids, i), j));
    }
    clusters = (double **)calloc(k, sizeof(double *));
    /*kmeans */
    max_iter = 300, count_var = 0, isequal = 1;
    while (count_var < max_iter && isequal == 1){
        vector_to_cluster(k);
        isequal = update_centroids();
        count_var++;
    }
    PyObject *output_matrix = PyList_New(0);
    for (i = 0; i < k; i++){
        PyObject *centroid = PyList_New(0);
        for (j = 0; j < vector_len; j++)
            PyList_Append(centroid, Py_BuildValue("d", centroids[i][j]));
        PyList_Append(output_matrix, centroid);
    }
    freearray(vector_list, (int)m);
    freearray(centroids, k);
    freearray(clusters, k);
    free(clustersindexes);
    return output_matrix;
}


/*Sorry about my english... this is not my strong side :)*/
static PyMethodDef the_module_methods[] = {
    {"fit",
     (PyCFunction)fit,
     METH_VARARGS,
     PyDoc_STR("the final centroids produced by the Kmeans algorithm")},
    {"WeightedAdjacencyMatrix",
     (PyCFunction)WeightedAdjacencyMatrix,
     METH_VARARGS,
     PyDoc_STR("get weighted adjacency of the mat")},
    {"DiagonalDegreeMatrix",
     (PyCFunction)DiagonalDegreeMatrix,
     METH_VARARGS,
     PyDoc_STR("get diag deg of the matrix")},
    {"NormalizedGraphLaplacian",
     (PyCFunction)NormalizedGraphLaplacian,
     METH_VARARGS,
     PyDoc_STR("just normalized graph matrix")},
    {"Jacobi",
     (PyCFunction)Jacobi,
     METH_VARARGS,
     PyDoc_STR("calculate of the eigenvalues and pairs with the Jacobi algorithm")},
    {"heuristic",
     (PyCFunction)heuristic,
     METH_VARARGS,
     PyDoc_STR("calculate of the k value with the eigengap heuristic method")},
    {"fullSpectralPy",
     (PyCFunction)fullSpectralPy,
     METH_VARARGS,
     PyDoc_STR("Calculate of the first part of the algorithm")},
    {"kmeans",
     (PyCFunction)kmeans,
     METH_VARARGS,
     PyDoc_STR("calculate of the second part of the algorithm")},
    {NULL, NULL, 0, NULL}};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "myspkmeans",
    NULL,
    -1,
    the_module_methods};

PyMODINIT_FUNC PyInit_myspkmeans(void){
    PyObject *m;
    m = PyModule_Create(&moduledef);
    if (!m)
        return NULL;
    else
        return m;
}
