/* IFT2425 - TP1 */
/* Vincent Foley-Bourgon (FOLV08078309) */
/* Eric Thivierge (THIE09016601) */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define N 20
#define SIZE (N-1)


typedef struct {
    float** elems;
    float*  start;
    int     rows, cols;
} matrix_t;



matrix_t* NewMatrix(int rows, int cols) {
    matrix_t* m = malloc(sizeof(matrix_t));
    if (m == NULL) {
        fprintf(stderr, "not enough memory\n");
        abort();
    }

    float** elems = malloc(sizeof(float*) * rows);
    if (elems == NULL) {
        fprintf(stderr, "not enough memory\n");
        abort();
    }

    float* data = calloc(rows*cols, sizeof(float));
    if (data == NULL) {
        fprintf(stderr, "not enough memory\n");
        abort();
    }

    for (int i = 0; i < rows; ++i) {
        elems[i] = data + (i*cols);
    }

    m->elems = elems;
    m->start = data;
    m->rows = rows;
    m->cols = cols;

    return m;
}


void FreeMatrix(matrix_t* m) {
    free(m->start);
    free(m->elems);
    free(m);
}



void MakeTridiagonalMatrix(matrix_t* matrix) {
    // Do nothing if the matrix isn't square.
    if (matrix->rows != matrix->cols)
        return;

    for (int i = 0; i < matrix->rows; ++i) {
        for (int j = 0; j < matrix->cols; ++j) {
            if (i == j) matrix->elems[i][j] = 2.0;
            else if (i == j+1 || j == i+1) matrix->elems[i][j] = -1.0;
        }
    }
}


void MakeBVector(float (*f)(float), matrix_t* vector) {
    float h = 1.0/N;

    for (int i = 0; i < vector->rows; ++i)
        vector->elems[i][0] = (*f)(h * (i+1));
}


void PrintMatrix(matrix_t* matrix) {
    for (int i = 0; i < matrix->rows; ++i) {
        for (int j = 0; j < matrix->cols; ++j) {
            printf("%+.2f ", matrix->elems[i][j]);
        }
        putchar('\n');
    }
}



float force(float x) {
    return x * (x - 1);
}


int FindMaxCoefficient(matrix_t* A, int col) {
    int max_row = col;
    float max = fabsf(A->elems[col][col]);
    for (int row = col+1; row < A->rows; ++row) {
        float x = fabsf(A->elems[row][col]);
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


void SwapRows(matrix_t* A, int i, int j) {
    float* tmp = A->elems[i];
    A->elems[i] = A->elems[j];
    A->elems[j] = tmp;
}

void ReplaceLine(matrix_t* A, int pivot_line, int replaced_line, float pivot, float k) {
    for (int col = pivot_line; col < A->cols; ++col) {
        A->elems[replaced_line][col] = A->elems[replaced_line][col] - (k/pivot)*A->elems[pivot_line][col];
    }
}

/*
  A: la matrice à factoriser; A va devenir la matrice triangulaire supérieure U.
  L: paramètre sortant contenant la matrice triangulaire inférieure L.
  pvect: paramètre sortant contenant le vecteur des permutations
 */
void PLUFactorize(matrix_t* A, matrix_t* L, int* pvect) {
    /* Initialiser pvect */
    for (int i = 0; i < A->rows; ++i)
        pvect[i] = i;

    for (int i = 0; i < A->rows; ++i) {
        /* Échanger la ligne courante avec celle possédant le plus
         * grand pivot (en valeur absolue). */
        int pivot_index = FindMaxCoefficient(A, i);
        Swap(pvect, i, pivot_index);
        SwapRows(A, i, pivot_index);

        /* Transcrire dans L le contenu de la colonne courante et
         * faire la division par le pivot. */
        float pivot = A->elems[i][i];
        for (int j = i; j < A->rows; j++) {
            L->elems[j][i] = A->elems[j][i] / pivot;
        }

        /* Appliquer Gauss aux autres lignes. */
        for (int row = i+1; row < A->rows; ++row) {
            ReplaceLine(A, i, row, pivot, A->elems[row][i]);
        }
    }
}

/* P must be zeroed out. */
void MakePermutationMatrix(int* pvect, matrix_t* P) {
    for (int i = 0; i < P->rows; ++i) {
        int j = pvect[i];
        P->elems[i][j] = 1;
    }
}


void SolveForward(matrix_t* A, matrix_t* x, matrix_t* b) {
    for (int i = 0; i < A->rows; ++i) {
        float sum = 0.0;
        for (int j = 0; j < i; ++j) {
            sum += A->elems[i][j] * x->elems[j][0];
        }
        x->elems[i][0] = (b->elems[i][0] - sum) / A->elems[i][i];
    }
}


void SolveBackward(matrix_t* A, matrix_t* x, matrix_t* b) {
    for (int i = A->rows - 1; i >= 0; --i) {
        float sum = 0.0;
        for (int j = i+1; j < A->cols; j++) {
            sum += A->elems[i][j] * x->elems[j][0];
        }
        x->elems[i][0] = (b->elems[i][0] - sum) / A->elems[i][i];
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


    matrix_t* A = NewMatrix(3, 3);
    matrix_t* L = NewMatrix(3, 3);
    matrix_t* P = NewMatrix(3, 3);
    int pvect[3];


    A->elems[0][0] = 1;
    A->elems[0][1] = 3;
    A->elems[0][2] = 6;
    A->elems[1][0] = 2;
    A->elems[1][1] = 4;
    A->elems[1][2] = 4;
    A->elems[2][0] = 3;
    A->elems[2][1] = 3;
    A->elems[2][2] = 3;

    PrintMatrix(A);
    PrintMatrix(L);

    putchar('\n');

    PLUFactorize(A, L, pvect);
    MakePermutationMatrix(pvect, P);
    PrintMatrix(A);
    PrintMatrix(L);
    PrintMatrix(P);
    for (int i = 0; i < 3; ++i)
        printf("%d ", pvect[i]);
    putchar('\n');

    printf("----\n");

    matrix_t* b = NewMatrix(3, 1);
    matrix_t* x = NewMatrix(3, 1);
    matrix_t* y = NewMatrix(3, 1);

    b->elems[0][0] = 18;
    b->elems[1][0] = 18;
    b->elems[2][0] = 6;

    SolveForward(L, y, b);

    PrintMatrix(y);

    SolveBackward(A, x, y);

    PrintMatrix(x);


    FreeMatrix(A);
    FreeMatrix(L);
    FreeMatrix(P);
    FreeMatrix(b);
    FreeMatrix(x);
    FreeMatrix(y);

    return 0;
}
