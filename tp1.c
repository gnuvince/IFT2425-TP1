/* IFT2425 - TP1 */
/* Vincent Foley-Bourgon (FOLV08078309) */
/* Eric Thivierge (THIE09016601) */

#include <stdlib.h>
#include <stdio.h>


#define SIZE 20


/*----------------------------------------------------------*/
/*  Alloue de la memoire pour une matrice 2d de float       */
/*----------------------------------------------------------*/
float** fmatrix_allocate_2d(int vsize,int hsize) {
    int i;
    float** matrix;
    float *imptr;

    matrix=(float**)calloc(vsize, sizeof(float*));
    if (matrix==NULL) printf("probleme d'allocation memoire");

    imptr=(float*)calloc(hsize*vsize, sizeof(float));
    if (imptr==NULL) printf("probleme d'allocation memoire");

    for(i=0; i<vsize; i++,imptr+=hsize) matrix[i]=imptr;
    return matrix;
}

/*----------------------------------------------------------*/
/* Libere la memoire de la matrice 2d de float              */
/*----------------------------------------------------------*/
void free_fmatrix_2d(float** pmat) {
    free(pmat[0]);
    free(pmat);
}


void MakeTridiagonalMatrix(int n, float** matrix) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i == j) matrix[i][j] = 2.0;
            else if (i == j+1 || j == i+1) matrix[i][j] = -1.0;
        }
    }
}


void MakeBVector(float (*f)(float), int n, float** vector) {
    float h = 1.0/n;

    for (int i = 0; i < n; ++i)
        vector[i][0] = (*f)(h * (i+1));
}


void PrintMatrix(int n, int m, float **matrix) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            printf("%+.2f ", matrix[i][j]);
        }
        putchar('\n');
    }
}



float force(float x) {
    return x * (x - 1);
}


int FindMaxCoefficient(float** A, int col) {
    double max = abs(A[col][col]);
    for (int row = col+1; row < SIZE; ++row) {
        double x = abs(A[row][col]);
        if (x > max) max = x;
    }
    return max;
}


void Swap(float* vect, int i, int j) {
    float tmp = vect[i];
    vect[i] = vect[j];
    vect[j] = tmp;
}


void SwapRows(float** A, int i, int j) {
    float* tmp = A[i];
    A[i] = A[j];
    A[j] = tmp;
}

/*
  A: la matrice à factoriser
  L: paramètre sortant contenant la matrice L
  U: paramètre sortant contenant la matrice U
  pvect: paramètre sortant contenant le vecteur des permutations
 */
void PLUFactorize(float** A, float** L, float** U, float* pvect) {
    /* Initialiser pvect */
    for (int i = 0; i < SIZE; ++i)
        pvect[i] = (float)i;

    for (int i = 0; i < SIZE; ++i) {
        /* Échanger la ligne courante avec celle possédant le plus
         * grand pivot (en valeur absolue). */
        int pivot_index = FindMaxCoefficient(A, i);
        Swap(pvect, i, pivot_index);
        SwapRows(A, i, pivot_index);

        /* Transcrire dans L le contenu de la colonne courante et
         * faire la division par le pivot. */
        float pivot = A[i][i];
        for (int j = i; j < SIZE; j++) {
            L[j][i] = A[j][i] / pivot;
        }


        /* Appliquer Gauss aux autres lignes. */
    }
}


int main(void) {
    float** matrix = fmatrix_allocate_2d(SIZE, SIZE);
    float** vector = fmatrix_allocate_2d(SIZE, 1);
    float (*f)(float) = &force;


    MakeTridiagonalMatrix(SIZE, matrix);
    MakeBVector(f, SIZE, vector);

    PrintMatrix(SIZE, SIZE, matrix);
    PrintMatrix(SIZE, 1, vector);

    free_fmatrix_2d(matrix);
    free_fmatrix_2d(vector);

    return 0;
}
