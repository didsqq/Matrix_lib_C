#include "s21_matrix.h"

int s21_create_matrix(int rows, int columns, matrix_t *result) {
  int error = 0;
  if (rows < 1 || columns < 1 || result == NULL) {
    error = INCORRECT_ERROR;
  } else if (s21_memory_allocation(rows, columns, result) != 0) {
    error = CALC_ERROR;
  } else {
    result->rows = rows;
    result->columns = columns;
    for (int i = 0; i < rows; i++) {
      for (int j = 0; j < columns; j++) {
        result->matrix[i][j] = 0;
      }
    }
  }
  return error;
}

void s21_remove_matrix(matrix_t *A) {
  if (A != NULL) {
    if (A->matrix != NULL) {
      for (int j = 0; j < A->rows; j++) {
        free(A->matrix[j]);
      }
      free(A->matrix);
    }
    A->columns = 0;
    A->rows = 0;
    A->matrix = NULL;
  }
}

int s21_eq_matrix(matrix_t *A, matrix_t *B) {
  int result = 1;
  if (s21_check_matrix(A) != OK || s21_check_matrix(B) != OK) {
    result = 0;
  } else if (A->rows != B->rows || A->columns != B->columns) {
    result = 0;
  } else {
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < A->columns; j++) {
        if (fabs(A->matrix[i][j] - B->matrix[i][j]) > 1e-7) {
          result = 0;
        }
      }
    }
  }
  return result;
}

int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int error = OK;
  if (s21_check_matrix(A) != OK || s21_check_matrix(B) != OK) {
    error = INCORRECT_ERROR;
  } else if (A->rows != B->rows || A->columns != B->columns) {
    error = CALC_ERROR;
  } else {
    error = s21_create_matrix(A->rows, A->columns, result);
    if (error == OK) {
      for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->columns; j++) {
          result->matrix[i][j] = A->matrix[i][j] + B->matrix[i][j];
        }
      }
    }
  }
  return error;
}

int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int error = OK;
  if (s21_check_matrix(A) != OK || s21_check_matrix(B) != OK) {
    error = INCORRECT_ERROR;
  } else if (A->rows != B->rows || A->columns != B->columns) {
    error = CALC_ERROR;
  } else {
    error = s21_create_matrix(A->rows, A->columns, result);
    if (error == OK) {
      for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->columns; j++) {
          result->matrix[i][j] = A->matrix[i][j] - B->matrix[i][j];
        }
      }
    }
  }
  return error;
}

int s21_mult_number(matrix_t *A, double number, matrix_t *result) {
  int error = OK;
  if (s21_check_matrix(A) != OK) {
    error = INCORRECT_ERROR;
  } else {
    error = s21_create_matrix(A->rows, A->columns, result);
    if (error == OK) {
      for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->columns; j++) {
          result->matrix[i][j] = A->matrix[i][j] * number;
        }
      }
    }
  }
  return error;
}

int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int error = OK;
  if (s21_check_matrix(A) != OK || s21_check_matrix(B) != OK)
    return INCORRECT_ERROR;
  if (A->columns != B->rows) return CALC_ERROR;

  error = s21_create_matrix(A->rows, B->columns, result);
  if (error == OK) {
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < B->columns; j++) {
        for (int k = 0; k < A->columns; k++) {
          result->matrix[i][j] += A->matrix[i][k] * B->matrix[k][j];
        }
      }
    }
  }

  return error;
}

int s21_transpose(matrix_t *A, matrix_t *result) {
  int error = OK;

  if (s21_check_matrix(A) != OK) {
    error = INCORRECT_ERROR;
  } else {
    error = s21_create_matrix(A->columns, A->rows, result);
    if (error == OK) {
      for (int i = 0; i < A->columns; i++) {
        for (int j = 0; j < A->rows; j++) {
          result->matrix[i][j] = A->matrix[j][i];
        }
      }
    }
  }
  return error;
}

int s21_determinant(matrix_t *A, double *result) {
  int error = OK;
  *result = 0.0;
  if (s21_check_matrix(A) != OK) {
    error = INCORRECT_ERROR;
  } else if (A->rows != A->columns) {
    error = CALC_ERROR;
  } else {
    matrix_t buf;
    error = s21_create_matrix(A->rows, A->columns, &buf);
    if (error == OK) {
      for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->columns; j++) {
          buf.matrix[i][j] = A->matrix[i][j];
        }
      }
      *result = s21_gauss(&buf);
    }
    s21_remove_matrix(&buf);
  }
  return error;
}

int s21_calc_complements(matrix_t *A, matrix_t *result) {
  int error = OK;
  if (s21_check_matrix(A) != OK) {
    error = INCORRECT_ERROR;
  } else if (A->rows != A->columns) {
    error = CALC_ERROR;
  } else if (A->rows == 1) {
    error = s21_create_matrix(A->rows, A->columns, result);
    if (error == OK) {
      result->matrix[0][0] = A->matrix[0][0];
    }
  } else {
    error = s21_create_matrix(A->rows, A->columns, result);
    if (error == OK) {
      for (int i = 0; i < A->columns; i++) {
        for (int j = 0; j < A->rows; j++) {
          result->matrix[i][j] = s21_minor(A, i, j);
        }
      }
    }
  }
  return error;
}

int s21_inverse_matrix(matrix_t *A, matrix_t *result) {
  int error = OK;
  if (s21_check_matrix(A) != OK) return INCORRECT_ERROR;

  if (A->rows != A->columns) return CALC_ERROR;

  double det = 0.0;
  error = s21_determinant(A, &det);

  if (A->rows == 1 && A->columns == 1) {
    if (fabs(det) < DBL_EPSILON)
      error = CALC_ERROR;
    else if (error == OK) {
      error = s21_create_matrix(A->rows, A->columns, result);
      if (error == OK) result->matrix[0][0] = 1 / A->matrix[0][0];
    }
  } else if (error == OK) {
    matrix_t minor, transpose_minor;
    if (fabs(det) < DBL_EPSILON) {
      error = CALC_ERROR;
    } else if (error == OK) {
      if (error == OK) {
        error = s21_calc_complements(A, &minor);
        if (error == OK) error = s21_transpose(&minor, &transpose_minor);
        if (error == OK) {
          s21_mult_number(&transpose_minor, 1 / det, result);
        }
        s21_remove_matrix(&minor);
        s21_remove_matrix(&transpose_minor);
      }
    }
  }
  return error;
}

double s21_minor(matrix_t *A, int row_delete, int column_delete) {
  double result = 0.0;
  matrix_t buf;
  s21_create_matrix(A->rows - 1, A->columns - 1, &buf);

  for (int i = 0, minor_i = 0; i < A->rows; i++) {
    if (i != row_delete) {
      for (int j = 0, minor_j = 0; j < A->columns; j++) {
        if (j != column_delete) {
          buf.matrix[minor_i][minor_j] = A->matrix[i][j];
          minor_j++;
        }
      }
      minor_i++;
    }
  }

  int error = s21_determinant(&buf, &result);
  if ((row_delete + column_delete) % 2 != 0 && error == OK) result *= -1;
  s21_remove_matrix(&buf);
  return result;
}

double s21_gauss(matrix_t *A) {
  int sign_determinant = 1;
  for (int i = 0; i < A->rows - 1; i++) {
    for (int k = i + 1; k < A->rows; k++) {
      if (fabs(A->matrix[i][i]) < DBL_EPSILON) {
        int swap_col_value = A->matrix[i][i];
        int swap_col_index = i;
        while (swap_col_index < A->columns - 1 && swap_col_value == 0) {
          swap_col_index++;
          swap_col_value = A->matrix[i][swap_col_index];
        }
        s21_swap_columns(A, i, swap_col_index);
        sign_determinant *= -1;
      }
      if (A->matrix[i][i] != 0.0) {
        double multiplier = A->matrix[k][i] / A->matrix[i][i];
        for (int j = i; j < A->columns; j++) {
          A->matrix[k][j] -= A->matrix[i][j] * multiplier;
        }
      }
    }
  }
  double result = A->matrix[0][0];
  for (int i = 1; i < A->rows; i++) {
    result *= A->matrix[i][i];
  }
  result *= sign_determinant;
  return result;
}

void s21_swap_columns(matrix_t *A, int row_index, int col_index) {
  matrix_t buf;
  s21_create_matrix(A->rows, 1, &buf);

  for (int i = 0; i < A->rows; i++) {
    buf.matrix[i][0] = A->matrix[i][row_index];
  }
  for (int i = 0; i < A->rows; i++) {
    A->matrix[i][row_index] = A->matrix[i][col_index];
  }
  for (int i = 0; i < A->rows; i++) {
    A->matrix[i][col_index] = buf.matrix[i][0];
  }
  s21_remove_matrix(&buf);
}

int s21_check_matrix(matrix_t *A) {
  int error = OK;
  if (A == NULL) {
    error = INCORRECT_ERROR;
  } else if (A->rows < 1 || A->columns < 1 || A->matrix == NULL) {
    error = INCORRECT_ERROR;
  }
  return error;
}

int s21_memory_allocation(int rows, int columns, matrix_t *result) {
  int error = OK;

  result->matrix = calloc(rows, sizeof(double *));
  if (result->matrix == NULL) error = CALC_ERROR;
  if (error == OK) {
    for (int i = 0; i < rows && error == OK; i++) {
      result->matrix[i] = calloc(columns, sizeof(double));
      if (result->matrix[i] == NULL) {
        error = CALC_ERROR;
        s21_remove_matrix(result);
      }
    }
  }

  return error;
}