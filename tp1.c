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

    matrix=(float**)malloc(sizeof(float*)*vsize);
    if (matrix==NULL) printf("probleme d'allocation memoire");

    imptr=(float*)malloc(sizeof(float)*hsize*vsize);
    if (imptr==NULL) printf("probleme d'allocation memoire");

    for(i=0;i<vsize;i++,imptr+=hsize) matrix[i]=imptr;
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
            else matrix[i][j] = 0.0;
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
