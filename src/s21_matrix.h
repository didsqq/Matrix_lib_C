#ifndef S21_MATRIX_H
#define S21_MATRIX_H

#define OK 0
#define INCORRECT_ERROR 1
#define CALC_ERROR 2

#define SUCCESS 1
#define FAILURE 0

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

typedef struct matrix_struct {
  double **matrix;
  int rows;
  int columns;
} matrix_t;

int s21_create_matrix(int rows, int columns, matrix_t *result);
int s21_eq_matrix(matrix_t *A, matrix_t *B);
int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result);
int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result);
int s21_mult_number(matrix_t *A, double number, matrix_t *result);
int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result);
int s21_transpose(matrix_t *A, matrix_t *result);
int s21_calc_complements(matrix_t *A, matrix_t *result);
int s21_determinant(matrix_t *A, double *result);
int s21_inverse_matrix(matrix_t *A, matrix_t *result);
void s21_remove_matrix(matrix_t *A);

int s21_memory_allocation(int rows, int columns, matrix_t *result);
int s21_check_matrix(matrix_t *A);
double s21_gauss(matrix_t *A);
void s21_swap_columns(matrix_t *A, int row, int col_index);
double s21_minor(matrix_t *A, int row_delete, int column_delete);

#endif