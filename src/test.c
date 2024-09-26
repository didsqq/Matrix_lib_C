#include <check.h>
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <unistd.h>

#include "s21_matrix.h"

#define OK 0
#define ERR_MAT 1
#define ERR_CAL 2

#define MAX_DOUBLE 1.79769e+308
#define MIN_DOUBLE 2.22507e-308

void testcases(Suite *testcase);
double get_rand(double min, double max);
matrix_t test;
matrix_t test2;
matrix_t result_sum;
matrix_t result_sub;
matrix_t result_mult_num;
matrix_t result_mult;
matrix_t result_trans;
matrix_t result_calc;
matrix_t result_inv;
double deter_result = 0;
int res = 0;

double calc_matrix[3][3] = {{0, 10, -20}, {4, -14, 8}, {-8, -2, 4}};

double inv_matrix[2][2] = {{-2.25, 1.25}, {1.75, -0.75}};

double trans_matrix[3][3] = {{1, 3, 5}, {2, 4, 6}};

double mult_num_matrix[3][2] = {{2, 4}, {6, 8}, {10, 12}};

void s21_fill_matrix(matrix_t *A, double val) {
  for (int i = 0; i < A->rows; i++) {
    for (int j = 0; j < A->columns; j++) {
      A->matrix[i][j] = val;
      val++;
    }
  }
}

START_TEST(create_1) {
  int rows = rand() % 100 + 1;
  int cols = rand() % 100 + 1;
  matrix_t m = {0};
  s21_create_matrix(rows, cols, &m);
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      m.matrix[i][j] = 0;
      ck_assert_ldouble_eq_tol(0, m.matrix[i][j], 1e-07);
    }
  }
  ck_assert_int_eq(m.rows, rows);
  ck_assert_int_eq(m.columns, cols);
  s21_remove_matrix(&m);
}
END_TEST

START_TEST(create_2) {
  int rows = 0;
  int cols = 10;
  matrix_t m = {0};
  ck_assert_int_eq(s21_create_matrix(rows, cols, &m), ERR_MAT);
  s21_remove_matrix(&m);
}
END_TEST

START_TEST(create_3) {
  int rows = 10;
  int cols = 0;

  matrix_t m = {0};
  ck_assert_int_eq(s21_create_matrix(rows, cols, &m), ERR_MAT);
  s21_remove_matrix(&m);
}
END_TEST

START_TEST(create1) {
  res = s21_create_matrix(0, 5, &test);
  ck_assert_int_eq(1, res);
}
END_TEST

START_TEST(create2) {
  res = s21_create_matrix(-2, 5, &test);
  ck_assert_int_eq(1, res);
}
END_TEST

START_TEST(create3) {
  res = s21_create_matrix(3, 3, NULL);
  ck_assert_int_eq(1, res);
}
END_TEST

START_TEST(create4) {
  res = s21_create_matrix(3, 3, &test);
  ck_assert_int_eq(0, res);
  s21_remove_matrix(&test);
}
END_TEST

START_TEST(create5) {
  s21_create_matrix(3, 3, &test);
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      test.matrix[i][j] = 0;
    }
  }
  s21_remove_matrix(&test);
}
END_TEST

START_TEST(transpose_matrix) {
  int rows = rand() % 100 + 1;
  int cols = rand() % 100 + 1;
  matrix_t m = {0};
  s21_create_matrix(rows, cols, &m);

  matrix_t check = {0};
  s21_create_matrix(cols, rows, &check);

  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      double rand_val = get_rand(-10e10, 10e10);
      m.matrix[i][j] = rand_val;
      check.matrix[j][i] = rand_val;
    }
  }

  matrix_t res = {0};
  ck_assert_int_eq(s21_transpose(&m, &res), OK);
  ck_assert_int_eq(s21_eq_matrix(&check, &res), SUCCESS);

  s21_remove_matrix(&m);
  s21_remove_matrix(&res);
  s21_remove_matrix(&check);
}
END_TEST

START_TEST(transpose_matrix2) {
  int rows = rand() % 100 + 1;
  int cols = rand() % 100 + 1;
  rows = -rows;
  cols = -cols;
  matrix_t m = {0};
  s21_create_matrix(rows, cols, &m);
  matrix_t check = {0};
  s21_create_matrix(cols, rows, &check);

  matrix_t res = {0};
  ck_assert_int_eq(s21_transpose(&m, &res), ERR_MAT);

  s21_remove_matrix(&m);
  s21_remove_matrix(&res);
  s21_remove_matrix(&check);
}
END_TEST

START_TEST(trans1) { ck_assert_int_eq(s21_transpose(NULL, &result_trans), 1); }
END_TEST

START_TEST(trans2) {
  s21_create_matrix(3, 2, &test);
  s21_fill_matrix(&test, 1);
  ck_assert_int_eq(s21_transpose(&test, &result_trans), 0);
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 3; j++) {
      ck_assert_double_eq_tol(trans_matrix[i][j], result_trans.matrix[i][j],
                              0.0000001);
    }
  }
  s21_remove_matrix(&test);
  s21_remove_matrix(&result_trans);
}
END_TEST

START_TEST(eq_matrix) {
  int rows = rand() % 100 + 1;
  int cols = rand() % 100 + 1;
  matrix_t m = {0};
  s21_create_matrix(rows, cols, &m);
  matrix_t mtx = {0};
  s21_create_matrix(rows, cols, &mtx);
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      double rand_val = get_rand(DBL_MIN, DBL_MAX);
      m.matrix[i][j] = rand_val;
      mtx.matrix[i][j] = rand_val;
    }
  }
  ck_assert_int_eq(s21_eq_matrix(&m, &mtx), SUCCESS);
  s21_remove_matrix(&m);
  s21_remove_matrix(&mtx);
}
END_TEST

START_TEST(not_eq) {
  matrix_t m = {0};
  int rows = rand() % 100 + 1;
  int cols = rand() % 100 + 1;
  s21_create_matrix(rows, cols, &m);
  matrix_t mtx = {0};
  int rows1 = rand() % 100 + 1;
  int cols1 = rand() % 100 + 1;
  s21_create_matrix(rows1, cols1, &mtx);

  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      m.matrix[i][j] = get_rand(DBL_MIN, DBL_MAX) + 1;
    }
  }
  for (int i = 0; i < rows1; i++) {
    for (int j = 0; j < cols1; j++) {
      mtx.matrix[i][j] = get_rand(DBL_MIN, DBL_MAX);
    }
  }
  ck_assert_int_eq(s21_eq_matrix(&m, &mtx), 0);
  s21_remove_matrix(&m);
  s21_remove_matrix(&mtx);
}
END_TEST

START_TEST(not_eq1) {
  matrix_t m = {0};
  int rows = rand() % 100 + 1;
  int cols = rand() % 100 + 1;
  rows = -rows;
  cols = -cols;
  s21_create_matrix(rows, cols, &m);
  matrix_t mtx = {0};
  int rows1 = rand() % 100 + 1;
  int cols1 = rand() % 100 + 1;
  s21_create_matrix(rows1, cols1, &mtx);
  ck_assert_int_eq(s21_eq_matrix(&m, &mtx), FAILURE);
  s21_remove_matrix(&m);
  s21_remove_matrix(&mtx);
}
END_TEST

START_TEST(zero_matrix) {
  matrix_t A = {0};
  matrix_t B = {0};
  int result = s21_eq_matrix(&A, &B);
  ck_assert_int_eq(0, result);
}

START_TEST(zero_matrix_1) {
  matrix_t A = {0};
  matrix_t B = {0};
  s21_create_matrix(0, 0, &A);
  s21_create_matrix(0, 0, &B);
  int result = s21_eq_matrix(&A, &B);
  ck_assert_int_eq(0, result);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}

START_TEST(casual_matrix_1) {
  matrix_t A = {0};
  matrix_t B = {0};
  s21_create_matrix(1, 1, &A);
  s21_create_matrix(1, 1, &B);
  A.matrix[0][0] = 1;
  B.matrix[0][0] = 1;
  int result = s21_eq_matrix(&A, &B);
  ck_assert_int_eq(1, result);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}

START_TEST(casual_matrix_2) {
  matrix_t A = {0};
  matrix_t B = {0};
  s21_create_matrix(1, 1, &A);
  s21_create_matrix(1, 1, &B);
  A.matrix[0][0] = 1;
  B.matrix[0][0] = 2;
  int result = s21_eq_matrix(&A, &B);
  ck_assert_int_eq(0, result);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}

START_TEST(casual_matrix_3) {
  matrix_t A = {0};
  matrix_t B = {0};
  s21_create_matrix(2, 2, &A);
  s21_create_matrix(2, 2, &B);
  A.matrix[0][0] = 1;
  A.matrix[0][1] = 2;
  A.matrix[1][0] = 3;
  A.matrix[1][1] = 4;
  B.matrix[0][0] = 1;
  B.matrix[0][1] = 2;
  B.matrix[1][0] = 3;
  B.matrix[1][1] = 4;
  int result = s21_eq_matrix(&A, &B);
  ck_assert_int_eq(1, result);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}

START_TEST(casual_matrix_4) {
  matrix_t A = {0};
  matrix_t B = {0};
  s21_create_matrix(2, 2, &A);
  s21_create_matrix(2, 2, &B);
  A.matrix[0][0] = 1;
  A.matrix[0][1] = 2;
  A.matrix[1][0] = 3;
  A.matrix[1][1] = 4;
  B.matrix[0][0] = 1;
  B.matrix[0][1] = 2;
  B.matrix[1][0] = 3;
  B.matrix[1][1] = 4;
  int result = s21_eq_matrix(&A, &B);
  ck_assert_int_eq(1, result);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}

START_TEST(casual_matrix_5) {
  matrix_t A = {0};
  matrix_t B = {0};
  s21_create_matrix(2, 2, &A);
  s21_create_matrix(2, 2, &B);
  A.matrix[0][0] = 1.01;
  A.matrix[0][1] = 2;
  A.matrix[1][0] = 3.05;
  A.matrix[1][1] = 4;
  B.matrix[0][0] = 1.01;
  B.matrix[0][1] = 2;
  B.matrix[1][0] = 3.05;
  B.matrix[1][1] = 4;
  int result = s21_eq_matrix(&A, &B);
  ck_assert_int_eq(1, result);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}

START_TEST(casual_matrix_6) {
  matrix_t A = {0};
  matrix_t B = {0};
  s21_create_matrix(2, 2, &A);
  s21_create_matrix(2, 2, &B);
  A.matrix[0][0] = 1.01;
  A.matrix[0][1] = -2;
  A.matrix[1][0] = 3.05;
  A.matrix[1][1] = -4;
  B.matrix[0][0] = 1.01;
  B.matrix[0][1] = -2;
  B.matrix[1][0] = 3.05;
  B.matrix[1][1] = -4;
  int result = s21_eq_matrix(&A, &B);
  ck_assert_int_eq(1, result);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}

START_TEST(casual_matrix_7) {
  matrix_t A = {0};
  matrix_t B = {0};
  s21_create_matrix(2, 2, &A);
  s21_create_matrix(2, 2, &B);
  A.matrix[0][0] = 1.00000000234;
  A.matrix[0][1] = 2;
  A.matrix[1][0] = 3.05;
  A.matrix[1][1] = 4;
  B.matrix[0][0] = 1;
  B.matrix[0][1] = 2;
  B.matrix[1][0] = 3.05;
  B.matrix[1][1] = 4;
  int result = s21_eq_matrix(&A, &B);
  ck_assert_int_eq(1, result);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}

START_TEST(casual_matrix_8) {
  matrix_t A = {0};
  matrix_t B = {0};
  s21_create_matrix(2, 2, &A);
  s21_create_matrix(2, 2, &B);
  A.matrix[0][0] = 1.0001;
  A.matrix[0][1] = 2;
  A.matrix[1][0] = 3.05;
  A.matrix[1][1] = 4;
  B.matrix[0][0] = 1;
  B.matrix[0][1] = 2;
  B.matrix[1][0] = 3.05;
  B.matrix[1][1] = 4;
  int result = s21_eq_matrix(&A, &B);
  ck_assert_int_eq(0, result);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}

START_TEST(sum_matrix) {
  int rows = rand() % 100 + 1;
  int cols = rand() % 100 + 1;
  matrix_t m = {0};
  s21_create_matrix(rows, cols, &m);
  matrix_t mtx = {0};
  s21_create_matrix(rows, cols, &mtx);
  matrix_t check = {0};
  s21_create_matrix(rows, cols, &check);

  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      m.matrix[i][j] = get_rand(DBL_MIN, DBL_MAX);
      mtx.matrix[i][j] = get_rand(DBL_MIN, DBL_MAX);
      check.matrix[i][j] = m.matrix[i][j] + mtx.matrix[i][j];
    }
  }
  matrix_t res = {0};

  ck_assert_int_eq(s21_sum_matrix(&m, &mtx, &res), OK);
  ck_assert_int_eq(s21_eq_matrix(&check, &res), SUCCESS);
  s21_remove_matrix(&m);
  s21_remove_matrix(&mtx);
  s21_remove_matrix(&res);
  s21_remove_matrix(&check);
}
END_TEST

START_TEST(sum_matrix1) {
  matrix_t m = {0};
  int rows = rand() % 100 + 1;
  int cols = rand() % 100 + 1;
  s21_create_matrix(rows, cols, &m);
  matrix_t mtx = {0};
  int rows1 = rand() % 100 + 1;
  int cols1 = rand() % 100 + 1;
  s21_create_matrix(rows1, cols1, &mtx);

  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      m.matrix[i][j] = get_rand(DBL_MIN, DBL_MAX) + 1;
    }
  }
  for (int i = 0; i < rows1; i++) {
    for (int j = 0; j < cols1; j++) {
      mtx.matrix[i][j] = get_rand(DBL_MIN, DBL_MAX);
    }
  }

  matrix_t res = {0};
  ck_assert_int_eq(s21_sum_matrix(&m, &mtx, &res), ERR_CAL);
  s21_remove_matrix(&m);
  s21_remove_matrix(&mtx);
  s21_remove_matrix(&res);
}
END_TEST

START_TEST(sum_matrix2) {
  int rows = rand() % 100 + 1;
  int cols = rand() % 100 + 1;
  rows = -rows;
  cols = -cols;
  matrix_t m = {0};
  s21_create_matrix(rows, cols, &m);
  matrix_t mtx = {0};
  s21_create_matrix(rows, cols, &mtx);

  matrix_t res = {0};

  ck_assert_int_eq(s21_sum_matrix(&m, &mtx, &res), ERR_MAT);
  s21_remove_matrix(&m);
  s21_remove_matrix(&mtx);
  s21_remove_matrix(&res);
}
END_TEST

START_TEST(sub_matrix) {
  int rows = rand() % 100 + 1;
  int cols = rand() % 100 + 1;
  matrix_t m = {0};
  s21_create_matrix(rows, cols, &m);
  matrix_t mtx = {0};
  s21_create_matrix(rows, cols, &mtx);
  matrix_t check = {0};
  s21_create_matrix(rows, cols, &check);

  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      m.matrix[i][j] = get_rand(DBL_MIN, DBL_MAX);
      mtx.matrix[i][j] = get_rand(DBL_MIN, DBL_MAX);
      check.matrix[i][j] = m.matrix[i][j] - mtx.matrix[i][j];
    }
  }
  matrix_t res = {0};
  ck_assert_int_eq(s21_sub_matrix(&m, &mtx, &res), OK);
  ck_assert_int_eq(s21_eq_matrix(&check, &res), SUCCESS);

  s21_remove_matrix(&m);
  s21_remove_matrix(&mtx);
  s21_remove_matrix(&res);
  s21_remove_matrix(&check);
}
END_TEST

START_TEST(sub_matrix1) {
  matrix_t m = {0};
  int rows = rand() % 100 + 1;
  int cols = rand() % 100 + 1;
  s21_create_matrix(rows, cols, &m);
  matrix_t mtx = {0};
  int rows1 = rand() % 100 + 1;
  int cols1 = rand() % 100 + 1;
  s21_create_matrix(rows1, cols1, &mtx);

  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      m.matrix[i][j] = get_rand(DBL_MIN, DBL_MAX) + 1;
    }
  }
  for (int i = 0; i < rows1; i++) {
    for (int j = 0; j < cols1; j++) {
      mtx.matrix[i][j] = get_rand(DBL_MIN, DBL_MAX);
    }
  }

  matrix_t res = {0};
  ck_assert_int_eq(s21_sub_matrix(&m, &mtx, &res), ERR_CAL);
  s21_remove_matrix(&m);
  s21_remove_matrix(&mtx);
  s21_remove_matrix(&res);
}
END_TEST

START_TEST(sub_matrix2) {
  int rows = rand() % 100 + 1;
  int cols = rand() % 100 + 1;
  rows = -rows;
  cols = -cols;
  matrix_t m = {0};
  s21_create_matrix(rows, cols, &m);
  matrix_t mtx = {0};
  s21_create_matrix(rows, cols, &mtx);

  matrix_t res = {0};
  ck_assert_int_eq(s21_sub_matrix(&m, &mtx, &res), ERR_MAT);

  s21_remove_matrix(&m);
  s21_remove_matrix(&mtx);
  s21_remove_matrix(&res);
}
END_TEST

START_TEST(simple_mult) {
  int rows = 2;
  int cols = 3;
  matrix_t m = {0};
  s21_create_matrix(rows, cols, &m);
  matrix_t mtx = {0};
  s21_create_matrix(cols, rows, &mtx);

  for (int i = 0, c = 1; i < rows; i++)
    for (int j = 0; j < cols; j++) m.matrix[i][j] = c++;

  for (int i = 0, c = 7; i < cols; i++)
    for (int j = 0; j < rows; j++) mtx.matrix[i][j] = c++;

  matrix_t check = {0};
  s21_create_matrix(m.rows, mtx.columns, &check);
  check.matrix[0][0] = 58;
  check.matrix[0][1] = 64;
  check.matrix[1][0] = 139;
  check.matrix[1][1] = 154;

  matrix_t res = {0};
  ck_assert_int_eq(s21_mult_matrix(&m, &mtx, &res), OK);
  int eq = s21_eq_matrix(&check, &res);

  ck_assert_int_eq(eq, SUCCESS);

  s21_remove_matrix(&m);
  s21_remove_matrix(&mtx);
  s21_remove_matrix(&res);
  s21_remove_matrix(&check);
}
END_TEST

START_TEST(mult_matrix2) {
  int rows = rand() % 100 + 1;
  int cols = rand() % 100 + 1;
  rows = -rows;
  cols = -cols;
  matrix_t m = {0};
  s21_create_matrix(rows, cols, &m);
  matrix_t mtx = {0};
  s21_create_matrix(cols, rows, &mtx);

  matrix_t check = {0};
  s21_create_matrix(m.rows, mtx.columns, &check);

  matrix_t res = {0};
  ck_assert_int_eq(s21_mult_matrix(&m, &mtx, &res), ERR_MAT);

  s21_remove_matrix(&m);
  s21_remove_matrix(&mtx);
  s21_remove_matrix(&res);
  s21_remove_matrix(&check);
}
END_TEST

START_TEST(mult_matrix3) {
  matrix_t m = {0};
  int rows = 2;
  int cols = 3;
  s21_create_matrix(rows, cols, &m);
  matrix_t mtx = {0};
  int rows1 = 4;
  int cols1 = 5;
  s21_create_matrix(rows1, cols1, &mtx);

  matrix_t check = {0};
  s21_create_matrix(m.rows, mtx.columns, &check);

  matrix_t res = {0};
  ck_assert_int_eq(s21_mult_matrix(&m, &mtx, &res), ERR_CAL);

  s21_remove_matrix(&m);
  s21_remove_matrix(&mtx);
  s21_remove_matrix(&res);
  s21_remove_matrix(&check);
}
END_TEST

START_TEST(mult_number_matrix) {
  int rows = rand() % 100 + 1;
  int cols = rand() % 100 + 1;
  matrix_t m = {0};
  s21_create_matrix(rows, cols, &m);
  matrix_t check = {0};
  s21_create_matrix(rows, cols, &check);
  double mult_number = get_rand(-10e5, 10e5);
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      m.matrix[i][j] = get_rand(DBL_MIN, DBL_MAX);
      check.matrix[i][j] = m.matrix[i][j] * mult_number;
    }
  }
  matrix_t res = {0};
  ck_assert_int_eq(s21_mult_number(&m, mult_number, &res), OK);
  ck_assert_int_eq(s21_eq_matrix(&check, &res), SUCCESS);
  s21_remove_matrix(&m);
  s21_remove_matrix(&res);
  s21_remove_matrix(&check);
}
END_TEST

START_TEST(mult_number_matrix2) {
  int rows = rand() % 100 + 1;
  int cols = rand() % 100 + 1;
  rows = -rows;
  cols = -cols;
  matrix_t m = {0};
  s21_create_matrix(rows, cols, &m);
  matrix_t check = {0};
  s21_create_matrix(rows, cols, &check);
  double mult_number = get_rand(-10e5, 10e5);

  matrix_t res = {0};
  ck_assert_int_eq(s21_mult_number(&m, mult_number, &res), ERR_MAT);

  s21_remove_matrix(&m);
  s21_remove_matrix(&res);
  s21_remove_matrix(&check);
}
END_TEST

START_TEST(mult_num1) {
  ck_assert_int_eq(s21_mult_number(NULL, 0.25, &result_mult_num), 1);
}
END_TEST

START_TEST(mult_num2) {
  s21_create_matrix(3, 2, &test);
  s21_fill_matrix(&test, 1);
  ck_assert_int_eq(s21_mult_number(&test, 2, &result_mult_num), 0);
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 2; j++) {
      ck_assert_double_eq_tol(mult_num_matrix[i][j],
                              result_mult_num.matrix[i][j], 0.0000001);
    }
  }
  s21_remove_matrix(&test);
  s21_remove_matrix(&result_mult_num);
}
END_TEST

START_TEST(test_incorrect) {
  matrix_t m = {0};
  matrix_t result = {0};
  int code = s21_calc_complements(&m, &result);
  ck_assert_int_eq(code, ERR_MAT);
}
END_TEST

START_TEST(test_not_sqare) {
  matrix_t m = {0};
  matrix_t result = {0};
  int codec = s21_create_matrix(3, 4, &m);
  if (codec == OK) {
    int code = s21_calc_complements(&m, &result);
    ck_assert_int_eq(code, ERR_CAL);
    s21_remove_matrix(&m);
  }
}
END_TEST

START_TEST(calc1) {
  ck_assert_int_eq(s21_calc_complements(NULL, &result_calc), 1);
}
END_TEST

START_TEST(calc3) {
  s21_create_matrix(3, 3, &test);
  test.matrix[0][0] = 1;
  test.matrix[0][1] = 2;
  test.matrix[0][2] = 3;
  test.matrix[1][0] = 0;
  test.matrix[1][1] = 4;
  test.matrix[1][2] = 2;
  test.matrix[2][0] = 5;
  test.matrix[2][1] = 2;
  test.matrix[2][2] = 1;
  ck_assert_int_eq(s21_calc_complements(&test, &result_calc), 0);
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      ck_assert_double_eq_tol(calc_matrix[i][j], result_calc.matrix[i][j],
                              0.0000001);
    }
  }
  s21_remove_matrix(&test);
  s21_remove_matrix(&result_calc);
}
END_TEST

START_TEST(test_s21_calc_complements_1) {
  matrix_t A;
  matrix_t B;
  s21_create_matrix(3, 3, &A);

  A.matrix[0][0] = 1;
  A.matrix[0][1] = 2;
  A.matrix[0][2] = 3;
  A.matrix[1][0] = 0;
  A.matrix[1][1] = 4;
  A.matrix[1][2] = 2;
  A.matrix[2][0] = 5;
  A.matrix[2][1] = 2;
  A.matrix[2][2] = 1;

  int res = s21_calc_complements(&A, &B);
  ck_assert_int_eq(B.matrix[0][0], 0);
  ck_assert_int_eq(B.matrix[0][1], 10);
  ck_assert_int_eq(B.matrix[0][2], -20);
  ck_assert_int_eq(B.matrix[1][0], 4);
  ck_assert_int_eq(B.matrix[1][1], -14);
  ck_assert_int_eq(B.matrix[1][2], 8);
  ck_assert_int_eq(B.matrix[2][0], -8);
  ck_assert_int_eq(B.matrix[2][1], -2);
  ck_assert_int_eq(B.matrix[2][2], 4);
  ck_assert_int_eq(res, 0);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}
END_TEST

START_TEST(test_s21_calc_complements_2) {
  matrix_t A;
  matrix_t B;
  s21_create_matrix(3, 3, &A);

  A.matrix[0][0] = 3;
  A.matrix[0][1] = 2;
  A.matrix[0][2] = 2;
  A.matrix[1][0] = 2;
  A.matrix[1][1] = 2;
  A.matrix[1][2] = 8;
  A.matrix[2][0] = 3;
  A.matrix[2][1] = 2;
  A.matrix[2][2] = 2;

  int res = s21_calc_complements(&A, &B);
  ck_assert_int_eq(B.matrix[0][0], -12);
  ck_assert_int_eq(B.matrix[0][1], 20);
  ck_assert_int_eq(B.matrix[0][2], -2);
  ck_assert_int_eq(B.matrix[1][0], 0);
  ck_assert_int_eq(B.matrix[1][1], 0);
  ck_assert_int_eq(B.matrix[1][2], 0);
  ck_assert_int_eq(B.matrix[2][0], 12);
  ck_assert_int_eq(B.matrix[2][1], -20);
  ck_assert_int_eq(B.matrix[2][2], 2);
  ck_assert_int_eq(res, 0);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}
END_TEST

START_TEST(test_s21_calc_complements_3) {
  matrix_t A;
  matrix_t B;
  s21_create_matrix(3, 1, &A);

  int res = s21_calc_complements(&A, &B);
  ck_assert_int_eq(res, 2);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(determinant1) {
  int size = 5;
  matrix_t m = {0};
  s21_create_matrix(size, size, &m);

  for (int i = 0; i < size; i++)
    for (int j = 0; j < size; j++) m.matrix[i][j] = j;
  double res = 0;
  int code = s21_determinant(&m, &res);
  ck_assert_double_eq_tol(res, 0, 1e-6);
  ck_assert_int_eq(code, OK);

  s21_remove_matrix(&m);
}
END_TEST

START_TEST(determinant2) {
  int size = 4;
  matrix_t m = {0};
  s21_create_matrix(size, size, &m);

  for (int i = 0; i < size; i++)
    for (int j = 0; j < size; j++) m.matrix[i][j] = j + i;

  double res = 0;
  int code = s21_determinant(&m, &res);
  ck_assert_double_eq_tol(res, 0, 1e-6);
  ck_assert_int_eq(code, OK);

  s21_remove_matrix(&m);
}
END_TEST

START_TEST(determinant3) {
  int size = 5;
  matrix_t m = {0};
  s21_create_matrix(size, size, &m);
  m.matrix[0][1] = 6;
  m.matrix[0][2] = -2;
  m.matrix[0][3] = -1;
  m.matrix[0][4] = 5;
  m.matrix[1][3] = -9;
  m.matrix[1][4] = -7;
  m.matrix[2][1] = 15;
  m.matrix[2][2] = 35;
  m.matrix[3][1] = -1;
  m.matrix[3][2] = -11;
  m.matrix[3][3] = -2;
  m.matrix[3][4] = 1;
  m.matrix[4][0] = -2;
  m.matrix[4][1] = -2;
  m.matrix[4][2] = 3;
  m.matrix[4][4] = -2;

  double res = 0;
  int code = s21_determinant(&m, &res);
  ck_assert_double_eq_tol(res, 2480, 1e-6);
  ck_assert_int_eq(code, OK);

  s21_remove_matrix(&m);
}
END_TEST

START_TEST(determinant4) {
  int size = 3;
  matrix_t m = {0};
  s21_create_matrix(size, size, &m);
  m.matrix[0][0] = 2;
  m.matrix[0][1] = 3;
  m.matrix[0][2] = 1;
  m.matrix[1][0] = 7;
  m.matrix[1][1] = 4;
  m.matrix[1][2] = 1;
  m.matrix[2][0] = 9;
  m.matrix[2][1] = -2;
  m.matrix[2][2] = 1;

  double res = 0;
  int code = s21_determinant(&m, &res);
  ck_assert_double_eq_tol(res, -32, 1e-6);
  ck_assert_int_eq(code, OK);

  s21_remove_matrix(&m);
}
END_TEST

START_TEST(determinant5) {
  int size = 2;
  matrix_t m = {0};
  s21_create_matrix(size, size, &m);
  m.matrix[0][0] = -5;
  m.matrix[0][1] = -4;
  m.matrix[1][0] = -2;
  m.matrix[1][1] = -3;

  double res = 0;
  int code = s21_determinant(&m, &res);
  ck_assert_double_eq_tol(res, 7, 1e-6);
  ck_assert_int_eq(code, OK);

  s21_remove_matrix(&m);
}
END_TEST

START_TEST(determinant6) {
  int size = 1;
  matrix_t m = {0};
  s21_create_matrix(size, size, &m);
  m.matrix[0][0] = -5;

  double res = 0;
  int code = s21_determinant(&m, &res);
  ck_assert_double_eq_tol(res, -5, 1e-6);
  ck_assert_int_eq(code, OK);

  s21_remove_matrix(&m);
}
END_TEST

START_TEST(determinant7) {
  matrix_t m = {0};
  int rows = rand() % 100 + 1;
  rows = -rows;
  s21_create_matrix(rows, rows, &m);
  double res = 0;
  int code = s21_determinant(&m, &res);
  ck_assert_int_eq(code, ERR_MAT);

  s21_remove_matrix(&m);
}
END_TEST

START_TEST(determinant8) {
  matrix_t m = {0};
  int rows = 4;
  int cols = 5;
  s21_create_matrix(rows, cols, &m);

  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      m.matrix[i][j] = get_rand(DBL_MIN, DBL_MAX) + 1;
    }
  }
  double res = 0;
  int code = s21_determinant(&m, &res);
  ck_assert_int_eq(code, ERR_CAL);

  s21_remove_matrix(&m);
}
END_TEST

START_TEST(deter1) {
  ck_assert_int_eq(s21_determinant(NULL, &deter_result), 1);
}
END_TEST

START_TEST(deter2) {
  s21_create_matrix(3, 2, &test);
  ck_assert_int_eq(s21_determinant(&test, &deter_result), 2);
  s21_remove_matrix(&test);
}
END_TEST

START_TEST(deter3) {
  s21_create_matrix(1, 1, &test);
  s21_fill_matrix(&test, 1);
  s21_determinant(&test, &deter_result);
  ck_assert_double_eq_tol(1, deter_result, 0.0000001);
  s21_remove_matrix(&test);
}
END_TEST

START_TEST(deter4) {
  s21_create_matrix(2, 2, &test);
  s21_fill_matrix(&test, 1);
  ck_assert_int_eq(s21_determinant(&test, &deter_result), 0);
  ck_assert_double_eq_tol(-2, deter_result, 0.0000001);
  s21_remove_matrix(&test);
}
END_TEST

START_TEST(deter5) {
  s21_create_matrix(3, 3, &test);
  s21_fill_matrix(&test, 3);
  s21_determinant(&test, &deter_result);
  ck_assert_double_eq_tol(0, deter_result, 0.0000001);
  s21_remove_matrix(&test);
}
END_TEST

START_TEST(s21_inverse_1) {
  matrix_t A = {0};
  matrix_t C = {0};
  s21_create_matrix(3, 3, &A);
  s21_create_matrix(3, 3, &C);
  C.matrix[0][0] = 1.0;
  C.matrix[0][1] = -1.0;
  C.matrix[0][2] = 1.0;
  C.matrix[1][0] = -38.0;
  C.matrix[1][1] = 41.0;
  C.matrix[1][2] = -34.0;
  C.matrix[2][0] = 27.0;
  C.matrix[2][1] = -29.0;
  C.matrix[2][2] = 24.0;
  A.matrix[0][0] = 2.0;
  A.matrix[0][1] = 5.0;
  A.matrix[0][2] = 7.0;
  A.matrix[1][0] = 6.0;
  A.matrix[1][1] = 3.0;
  A.matrix[1][2] = 4.0;
  A.matrix[2][0] = 5.0;
  A.matrix[2][1] = -2.0;
  A.matrix[2][2] = -3.0;
  matrix_t B = {0};
  s21_inverse_matrix(&A, &B);
  int res = s21_eq_matrix(&B, &C);
  ck_assert_int_eq(res, 1);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&C);
}
END_TEST

START_TEST(determinant) {
  int size = 2;
  matrix_t m = {0};
  s21_create_matrix(size, size, &m);
  m.matrix[0][0] = 1;
  m.matrix[0][1] = 1;
  m.matrix[1][0] = 1;
  m.matrix[1][1] = 1;

  matrix_t result = {0};
  int code = s21_inverse_matrix(&m, &result);
  ck_assert_int_eq(code, ERR_CAL);

  s21_remove_matrix(&m);
}
END_TEST

START_TEST(inverse_one) {
  int size = 1;
  matrix_t m = {0};
  s21_create_matrix(size, size, &m);

  m.matrix[0][0] = 2;

  matrix_t res = {0};
  s21_inverse_matrix(&m, &res);

  matrix_t expected = {0};
  s21_create_matrix(size, size, &expected);
  expected.matrix[0][0] = 0.5;
  ck_assert_int_eq(s21_eq_matrix(&expected, &res), SUCCESS);
  s21_remove_matrix(&expected);
  s21_remove_matrix(&res);
  s21_remove_matrix(&m);
}

START_TEST(inverse) {
  int size = 3;
  matrix_t m = {0};
  s21_create_matrix(size, size, &m);

  m.matrix[0][0] = 2;
  m.matrix[0][1] = 5;
  m.matrix[0][2] = 7;
  m.matrix[1][0] = 6;
  m.matrix[1][1] = 3;
  m.matrix[1][2] = 4;
  m.matrix[2][0] = 5;
  m.matrix[2][1] = -2;
  m.matrix[2][2] = -3;

  matrix_t res = {0};
  s21_inverse_matrix(&m, &res);

  matrix_t expected = {0};
  s21_create_matrix(size, size, &expected);
  expected.matrix[0][0] = 1;
  expected.matrix[0][1] = -1;
  expected.matrix[0][2] = 1;
  expected.matrix[1][0] = -38;
  expected.matrix[1][1] = 41;
  expected.matrix[1][2] = -34;
  expected.matrix[2][0] = 27;
  expected.matrix[2][1] = -29;
  expected.matrix[2][2] = 24;
  ck_assert_int_eq(s21_eq_matrix(&expected, &res), SUCCESS);

  s21_remove_matrix(&expected);
  s21_remove_matrix(&res);
  s21_remove_matrix(&m);
}
END_TEST

START_TEST(inv1) { ck_assert_int_eq(s21_inverse_matrix(NULL, &result_inv), 1); }
END_TEST

START_TEST(inv2) {
  s21_create_matrix(1, 2, &test);
  ck_assert_int_eq(s21_inverse_matrix(&test, &result_inv), 2);
  s21_remove_matrix(&test);
}
END_TEST

START_TEST(inv3) {
  s21_create_matrix(2, 2, &test);
  s21_fill_matrix(&test, 1.5);
  ck_assert_int_eq(s21_inverse_matrix(&test, &result_inv), 0);
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      ck_assert_double_eq_tol(inv_matrix[i][j], result_inv.matrix[i][j],
                              0.0000001);
    }
  }
  s21_remove_matrix(&test);
  s21_remove_matrix(&result_inv);
}
END_TEST

START_TEST(inv4) {
  s21_create_matrix(3, 3, &test);
  s21_fill_matrix(&test, 1);
  ck_assert_int_eq(s21_inverse_matrix(&test, &result_inv), 2);
  s21_remove_matrix(&test);
}
END_TEST

START_TEST(inv5) {
  matrix_t m = {0};
  matrix_t result = {0};
  int codec = s21_create_matrix(1, 4, &m);
  if (!codec) {
    int code = s21_inverse_matrix(&m, &result);
    ck_assert_int_eq(code, 2);
    s21_remove_matrix(&m);
  }
}
END_TEST

START_TEST(inv6) {
  matrix_t m = {0};
  matrix_t result = {0};
  int code = s21_inverse_matrix(&m, &result);
  ck_assert_int_eq(code, 1);
}
END_TEST

START_TEST(inv7) {
  matrix_t m = {0};
  matrix_t result = {0};
  int codec = s21_create_matrix(1, 1, &m);
  if (!codec) {
    int code = s21_inverse_matrix(&m, &result);
    ck_assert_int_eq(code, 2);
    s21_remove_matrix(&m);
  }
}
END_TEST

Suite *suite_inverse_matrix(void) {
  Suite *s = suite_create("suite_inverse_matrix");
  TCase *tc = tcase_create("case_inverse_matrix");

  tcase_add_test(tc, inverse_one);
  tcase_add_test(tc, inverse);
  tcase_add_test(tc, s21_inverse_1);
  tcase_add_test(tc, determinant);
  tcase_add_test(tc, test_incorrect);
  tcase_add_test(tc, inv1);
  tcase_add_test(tc, inv2);
  tcase_add_test(tc, inv3);
  tcase_add_test(tc, inv4);
  tcase_add_test(tc, inv5);
  tcase_add_test(tc, inv6);
  tcase_add_test(tc, inv7);

  suite_add_tcase(s, tc);
  return s;
}

Suite *suite_determinant(void) {
  Suite *s = suite_create("suite_determinant");
  TCase *tc = tcase_create("case_determinant");

  tcase_add_test(tc, determinant1);
  tcase_add_test(tc, determinant2);
  tcase_add_test(tc, determinant3);
  tcase_add_test(tc, determinant4);
  tcase_add_test(tc, determinant5);
  tcase_add_test(tc, determinant6);
  tcase_add_loop_test(tc, determinant7, 0, 100);
  tcase_add_test(tc, determinant8);
  tcase_add_test(tc, deter1);
  tcase_add_test(tc, deter2);
  tcase_add_test(tc, deter3);
  tcase_add_test(tc, deter4);
  tcase_add_test(tc, deter5);

  suite_add_tcase(s, tc);
  return s;
}

Suite *suite_calc_complements(void) {
  Suite *s = suite_create("suite_calc_complements");
  TCase *tc = tcase_create("case_calc_complements");

  tcase_add_test(tc, test_not_sqare);
  tcase_add_test(tc, test_s21_calc_complements_1);
  tcase_add_test(tc, test_s21_calc_complements_2);
  tcase_add_test(tc, test_s21_calc_complements_3);
  tcase_add_test(tc, calc1);
  tcase_add_test(tc, calc3);
  tcase_add_test(tc, test_incorrect);

  suite_add_tcase(s, tc);
  return s;
}

Suite *suite_mult_number_matrix(void) {
  Suite *s = suite_create("suite_mult_number_matrix");
  TCase *tc = tcase_create("case_mult_number_matrix");

  tcase_add_loop_test(tc, mult_number_matrix, 0, 100);
  tcase_add_loop_test(tc, mult_number_matrix2, 0, 100);
  tcase_add_test(tc, mult_num1);
  tcase_add_test(tc, mult_num2);

  suite_add_tcase(s, tc);
  return s;
}

Suite *suite_mult_matrix(void) {
  Suite *s = suite_create("suite_mult_matrix");
  TCase *tc = tcase_create("case_mult_matrix");

  tcase_add_loop_test(tc, mult_matrix2, 0, 100);
  tcase_add_test(tc, mult_matrix3);
  tcase_add_test(tc, simple_mult);

  suite_add_tcase(s, tc);
  return s;
}

Suite *suite_sub_matrix(void) {
  Suite *s = suite_create("suite_sub_matrix");
  TCase *tc = tcase_create("case_sub_matrix");

  tcase_add_loop_test(tc, sub_matrix, 0, 100);
  tcase_add_loop_test(tc, sub_matrix1, 0, 100);
  tcase_add_loop_test(tc, sub_matrix2, 0, 100);

  suite_add_tcase(s, tc);
  return s;
}

Suite *suite_sum_matrix(void) {
  Suite *s = suite_create("suite_sum_matrix");
  TCase *tc = tcase_create("case_sum_matrix");

  tcase_add_loop_test(tc, sum_matrix, 0, 100);
  tcase_add_loop_test(tc, sum_matrix1, 0, 100);
  tcase_add_loop_test(tc, sum_matrix2, 0, 100);

  suite_add_tcase(s, tc);
  return s;
}

Suite *suite_eq_matrix(void) {
  Suite *s = suite_create("suite_eq_matrix");
  TCase *tc = tcase_create("case_eq_matrix");

  tcase_add_test(tc, not_eq);
  tcase_add_loop_test(tc, eq_matrix, 0, 100);
  tcase_add_test(tc, not_eq1);

  tcase_add_test(tc, zero_matrix);
  tcase_add_test(tc, zero_matrix_1);
  tcase_add_test(tc, casual_matrix_1);
  tcase_add_test(tc, casual_matrix_2);
  tcase_add_test(tc, casual_matrix_3);
  tcase_add_test(tc, casual_matrix_4);
  tcase_add_test(tc, casual_matrix_5);
  tcase_add_test(tc, casual_matrix_6);
  tcase_add_test(tc, casual_matrix_7);
  tcase_add_test(tc, casual_matrix_8);

  suite_add_tcase(s, tc);
  return s;
}

Suite *suite_transpose_matrix(void) {
  Suite *s = suite_create("suite_transpose_matrix");
  TCase *tc = tcase_create("case_transpose_matrix");

  tcase_add_loop_test(tc, transpose_matrix, 0, 100);
  tcase_add_loop_test(tc, transpose_matrix2, 0, 100);
  tcase_add_test(tc, trans1);
  tcase_add_test(tc, trans2);

  suite_add_tcase(s, tc);
  return s;
}

Suite *suite_create_matrix(void) {
  Suite *s = suite_create("suite_create_matrix");
  TCase *tc = tcase_create("case_create_matrix");

  tcase_add_test(tc, create_1);
  tcase_add_test(tc, create_2);
  tcase_add_test(tc, create_3);
  tcase_add_test(tc, create1);
  tcase_add_test(tc, create2);
  tcase_add_test(tc, create3);
  tcase_add_test(tc, create4);
  tcase_add_test(tc, create5);

  suite_add_tcase(s, tc);
  return s;
}

void testcases(Suite *testcase) {
  static int counter_testcase = 1;

  if (counter_testcase > 0) putchar('\n');
  counter_testcase++;
  SRunner *sr = srunner_create(testcase);

  srunner_set_fork_status(sr, CK_NOFORK);
  srunner_run_all(sr, CK_NORMAL);

  srunner_free(sr);
}

double get_rand(double min, double max) {
  double val = (double)rand() / RAND_MAX;
  return min + val * (max - min);
}

int main(void) {
  Suite *list_cases[] = {suite_create_matrix(),
                         suite_eq_matrix(),
                         suite_sub_matrix(),
                         suite_sum_matrix(),
                         suite_mult_matrix(),
                         suite_determinant(),
                         suite_mult_number_matrix(),
                         suite_transpose_matrix(),
                         suite_calc_complements(),
                         suite_inverse_matrix(),
                         NULL};
  for (Suite **current_testcase = list_cases; *current_testcase != NULL;
       current_testcase++) {
    testcases(*current_testcase);
  }
  return 0;
}