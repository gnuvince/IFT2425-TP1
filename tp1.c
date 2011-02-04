/* IFT2425 - TP1 */
/* Vincent Foley-Bourgon (FOLV08078309) */
/* Eric Thivierge (THIE09016601) */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define N 20
#define SIZE (N-1)


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


void MakeTridiagonalMatrix(float** matrix) {
    for (int i = 0; i < SIZE; ++i) {
        for (int j = 0; j < SIZE; ++j) {
            if (i == j) matrix[i][j] = 2.0;
            else if (i == j+1 || j == i+1) matrix[i][j] = -1.0;
        }
    }
}


void MakeBVector(float (*f)(float), float** vector) {
    float h = 1.0/N;

    for (int i = 0; i < SIZE; ++i)
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


int FindMaxCoefficient(float** A, int dim, int col) {
    int max_row = col;
    float max = fabsf(A[col][col]);
    for (int row = col+1; row < dim; ++row) {
        float x = fabsf(A[row][col]);
        if (x > max) {
            max = x;
            max_row = row;
        }
    }
    return max_row;
}


void Swap(int* vect, int i, int j) {
    int tmp = vect[i];
    vect[i] = vect[j];
    vect[j] = tmp;
}


void SwapRows(float** A, int i, int j) {
    float* tmp = A[i];
    A[i] = A[j];
    A[j] = tmp;
}

void ReplaceLine(float** A, int dim, int pivot_line, int replaced_line, float pivot, float k) {
    for (int col = pivot_line; col < dim; ++col) {
        A[replaced_line][col] = pivot*A[replaced_line][col] - k*A[pivot_line][col];
    }
}

/*
  A: la matrice à factoriser; A va devenir la matrice triangulaire supérieure U.
  L: paramètre sortant contenant la matrice triangulaire inférieure L.
  pvect: paramètre sortant contenant le vecteur des permutations
 */
void PLUFactorize(float** A, int dim, float** L, int* pvect) {
    /* Initialiser pvect */
    for (int i = 0; i < dim; ++i)
        pvect[i] = i;

    for (int i = 0; i < dim; ++i) {
        /* Échanger la ligne courante avec celle possédant le plus
         * grand pivot (en valeur absolue). */
        int pivot_index = FindMaxCoefficient(A, dim, i);
        Swap(pvect, i, pivot_index);
        SwapRows(A, i, pivot_index);

        /* Transcrire dans L le contenu de la colonne courante et
         * faire la division par le pivot. */
        float pivot = A[i][i];
        for (int j = i; j < dim; j++) {
            L[j][i] = A[j][i] / pivot;
        }

        /* Appliquer Gauss aux autres lignes. */
        for (int row = i+1; row < dim; ++row) {
            ReplaceLine(A, dim, i, row, pivot, A[row][i]);
        }
    }
}

/* P must be zeroed out. */
void MakePermutationMatrix(int* pvect, int dim, float** P) {
    for (int i = 0; i < dim; ++i) {
        int j = pvect[i];
        P[i][j] = 1;
    }
}


int main(void) {
    /*
    float** matrix = fmatrix_allocate_2d(SIZE, SIZE);
    float** vector = fmatrix_allocate_2d(SIZE, 1);
    float (*f)(float) = &force;


    MakeTridiagonalMatrix(matrix);
    MakeBVector(f, vector);

    PrintMatrix(SIZE, SIZE, matrix);
    PrintMatrix(SIZE, 1, vector);

    free_fmatrix_2d(matrix);
    free_fmatrix_2d(vector);
    */


    float** A = fmatrix_allocate_2d(3, 3);
    float** L = fmatrix_allocate_2d(3, 3);
    float** P = fmatrix_allocate_2d(3, 3);
    int pvect[3];

    A[0][0] = 1;
    A[0][1] = 3;
    A[0][2] = 6;
    A[1][0] = 2;
    A[1][1] = 4;
    A[1][2] = 4;
    A[2][0] = 3;
    A[2][1] = 3;
    A[2][2] = 3;

    PrintMatrix(3, 3, A);
    PrintMatrix(3, 3, L);

    putchar('\n');

    PLUFactorize(A, 3, L, pvect);
    MakePermutationMatrix(pvect, 3, P);
    PrintMatrix(3, 3, A);
    PrintMatrix(3, 3, L);
    PrintMatrix(3, 3, P);
    for (int i = 0; i < 3; ++i)
        printf("%d ", pvect[i]);
    putchar('\n');

    return 0;
}
