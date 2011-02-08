/*  IFT2425 - TP1 */
/* Vincent Foley-Bourgon (FOLV08078309) */
/* Eric Thivierge (THIE09016601) */

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define N 20
#define SIZE (N-1)
#define PI 3.14159265358979323846

/*
  Structure pour représenter une matrice:
  - elems: pointeurs vers les rangées
  - start: pointeur vers le début des données
           (pour libérer la mémoire et copier les données)
  - rows, cols: dimensions de la matrice.
 */
typedef struct {
    float** elems;
    float*  start;
    int     rows, cols;
} matrix_t;


/* Allocation dynamique d'une nouvelle matrice.  Tous les éléments
 * sont initialisés à 0. */
matrix_t* NewMatrix(int rows, int cols) {
    matrix_t* m = malloc(sizeof(matrix_t));
    if (m == NULL) {
        fprintf(stderr, "memoire insuffisante\n");
        abort();
    }

    float** elems = malloc(sizeof(float*) * rows);
    if (elems == NULL) {
        fprintf(stderr, "memoire insuffisante\n");
        abort();
    }

    float* data = calloc(rows*cols, sizeof(float));
    if (data == NULL) {
        fprintf(stderr, "memoire insuffisante\n");
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

/* Libération de la mémoire utilisée par une matrice. */
void FreeMatrix(matrix_t* m) {
    free(m->start);
    free(m->elems);
    free(m);
}


/* Créer une matrice ayant le format tri-diagonal:
  2  -1  0 ... 0
  -1  2 -1 ... 0
  0  -1  2 ... 0
  .
  .
  0 ...   0 -1 2
 */
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


/* Populer un vecteur en appliquant le pointeur de fonction f aux
 * éléments de l'intervalle {1/N, 2/N, ..., (N-1)/N} */
void MakeBVector(float (*f)(float), matrix_t* vector) {
    float h = 1.0/N;

    for (int i = 0; i < vector->rows; ++i)
        vector->elems[i][0] = (*f)(h * (i+1));
}


/* Afficher une matrice. */
void PrintMatrix(matrix_t* matrix) {
    for (int i = 0; i < matrix->rows; ++i) {
        for (int j = 0; j < matrix->cols; ++j) {
            printf("%+.2f ", matrix->elems[i][j]);
        }
        putchar('\n');
    }
    putchar('\n');
}



/* Trouver le plus grand coefficient d'une colonne en valeur absolue.
 * Le coefficient est recherché à partir de la ligne correspondant à
 * la colonne (ex. ligne 3, colonne 3). */
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


/* Échanger deux entiers dans un int[] */
void Swap(int* vect, int i, int j) {
    int tmp = vect[i];
    vect[i] = vect[j];
    vect[j] = tmp;
}


/* Échanger deux rangées dans une matrices en échangeant leurs
 * pointeurs. */
void SwapRows(matrix_t* A, int i, int j) {
    float* tmp = A->elems[i];
    A->elems[i] = A->elems[j];
    A->elems[j] = tmp;
}


/* Faire un remplacement de ligne sans faire de cadrage (nécessaire
 * pour la factorisation LU). */
void ReplaceLine(matrix_t* M, int pivot_line, int replaced_line, float pivot, float k) {
    for (int col = pivot_line; col < M->cols; ++col) {
        M->elems[replaced_line][col] = M->elems[replaced_line][col] - (k/pivot)*M->elems[pivot_line][col];
    }
}


/* Copier les données d'une matrice dans une autre. */
void CopyMatrix(matrix_t* dst, matrix_t* src) {
    if (src->rows != dst->rows || src->cols != dst->cols)
        return;

    memcpy(dst->start, src->start, sizeof(float) * src->rows * src->cols);
}

/* Faire une factorisation PLU d'une matrice A.
  A: la matrice à factoriser
  U: paramètre sortant contenant la matrice triangulaire supérieure U.
  L: paramètre sortant contenant la matrice triangulaire inférieure L.
  pvect: paramètre sortant contenant le vecteur des permutations
 */
void PLUFactorize(matrix_t* A, matrix_t* L, matrix_t* U, int* pvect) {
    /* Initialiser pvect */
    for (int i = 0; i < A->rows; ++i)
        pvect[i] = i;

    /* Copier les éléments de A dans U. */
    CopyMatrix(U, A);

    for (int i = 0; i < U->rows; ++i) {
        /* Échanger la ligne courante avec celle possédant le plus
         * grand pivot (en valeur absolue). */
        int pivot_index = FindMaxCoefficient(U, i);
        Swap(pvect, i, pivot_index);
        SwapRows(U, i, pivot_index);

        /* Transcrire dans L le contenu de la colonne courante et
         * faire la division par le pivot. */
        float pivot = U->elems[i][i];
        for (int j = i; j < U->rows; j++) {
            L->elems[j][i] = U->elems[j][i] / pivot;
        }

        /* Appliquer Gauss aux autres lignes. */
        for (int row = i+1; row < U->rows; ++row) {
            ReplaceLine(U, i, row, pivot, U->elems[row][i]);
        }
    }
}


/* À l'aide d'un vecteur des permutations, créer la matrice des
 * permutations.
 * N.B.: P doit contenur que des 0 initialement. */
void MakePermutationMatrix(int* pvect, matrix_t* P) {
    for (int i = 0; i < P->rows; ++i) {
        int j = pvect[i];
        P->elems[i][j] = 1;
    }
}


/* Résolution par le bas d'un système triangulaire inférieur. */
void SolveForward(matrix_t* A, matrix_t* x, matrix_t* b) {
    for (int i = 0; i < A->rows; ++i) {
        float sum = 0.0;
        for (int j = 0; j < i; ++j) {
            sum += A->elems[i][j] * x->elems[j][0];
        }
        x->elems[i][0] = (b->elems[i][0] - sum) / A->elems[i][i];
    }
}

/* Résolution par le haut d'un système triangulaire supérieur. */
void SolveBackward(matrix_t* A, matrix_t* x, matrix_t* b) {
    for (int i = A->rows - 1; i >= 0; --i) {
        float sum = 0.0;
        for (int j = i+1; j < A->cols; j++) {
            sum += A->elems[i][j] * x->elems[j][0];
        }
        x->elems[i][0] = (b->elems[i][0] - sum) / A->elems[i][i];
    }
}


/* Résoudre Ax = b à l'aide de la factorisation PLU. x, L, U et pvect seront modifiés. */
void SolvePLU(matrix_t* A, matrix_t* x, matrix_t* b, matrix_t* L, matrix_t* U, int* pvect) {
    matrix_t* y = NewMatrix(x->rows, 1);

    PLUFactorize(A, L, U, pvect);
    SolveForward(L, y, b);
    SolveBackward(U, x, y);

    FreeMatrix(y);
}



/* Fonction de la question 1. */
float force(float x) {
    return x * (x - 1);
}


/* Fonction de la question 2. */
float force2(float x) {
    float y = 2 * PI * x;
    return x * sin(y * y);
}

void MatrixMult(matrix_t* A, matrix_t* B, matrix_t* C) {
    if ((A->cols != B->rows) || (C->rows != A->rows) || (C->cols != B->cols))
        return;

    for (int i = 0; i < A->rows; i++)
        for (int j = 0; j < B->cols; j++) {
            float temp = 0;
            for (int k = 0; k < B->rows; k++)
                temp += A->elems[i][k] * B->elems[k][j];
            C->elems[i][j] = temp;
        }
}

matrix_t* MakeI(int n) {
  matrix_t* I = NewMatrix(n, n);

  for (int i = 0; i < n; i++)
    I->elems[i][i] = 1;

  return I;
}

int MatrixEq(matrix_t* A, matrix_t* B) {
  if ((A->rows != B-> rows) || (A->cols != B->cols))
    return 0;

  for (int i = 0; i < A->rows; i++)
    for (int j = 0; j < A->cols; j++)
      if (A->elems[i][j] != B->elems[i][j])
        return 0;

  return 1;
}


int main(void) {
    matrix_t* A = NewMatrix(SIZE, SIZE);
    matrix_t* LU = NewMatrix(SIZE, SIZE);
    matrix_t* b = NewMatrix(SIZE, 1);
    float (*f)(float) = &force;

    MakeTridiagonalMatrix(A);
    MakeBVector(f, b);


    PrintMatrix(A);
/*     PrintMatrix(b); */

    matrix_t* L = NewMatrix(SIZE, SIZE);
    matrix_t* U = NewMatrix(SIZE, SIZE);
    matrix_t* P = NewMatrix(SIZE, SIZE);
    matrix_t* x = NewMatrix(SIZE, 1);
    int pvect[SIZE];

    SolvePLU(A, x, b, L, U, pvect);
    PrintMatrix(L);
    PrintMatrix(U);


    MakePermutationMatrix(pvect, P);

    PrintMatrix(L);
    PrintMatrix(U);
    MatrixMult(L, U, LU);
    PrintMatrix(LU);
    printf("%d\n", MatrixEq(A, LU));
    exit(0);

    PrintMatrix(P);
    PrintMatrix(x);

    FreeMatrix(A);
    FreeMatrix(b);
    FreeMatrix(L);
    FreeMatrix(U);
    FreeMatrix(LU);
    FreeMatrix(P);
    FreeMatrix(x);

    return 0;
}
