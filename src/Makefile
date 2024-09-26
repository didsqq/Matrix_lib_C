CC=gcc
CFLAGS=-c -Wall -Werror -Wextra -std=c11
SOURCE = s21_*.c

all: s21_matrix.a

s21_matrix.a:
	gcc $(CFLAGS) $(SOURCE)
	ar rc s21_matrix.a s21_*.o
	ranlib s21_matrix.a

test: test.o s21_matrix.a
	$(CC) test.o s21_matrix.a -o test -lcheck
	./test

gcov_report:
	$(CC) --coverage s21_matrix.c test.c -o s21_test -lcheck
	./s21_test
	lcov -t "s21_test" -o s21_test.info -c -d .
	genhtml -o report s21_test.info
	open ./report/index.html

s21_matrix.o: s21_matrix.c
	$(CC) $(CFLAGS) s21_matrix.c -o s21_matrix.o

test.o: test.c
	$(CC) $(CFLAGS) test.c -o test.o

clean: 
	-rm -rf *.o *.out test s21_matrix.a *.gcda *.gcno s21_test s21_test.info report

.PHONY: test
.PHONY: s21_matrix.a