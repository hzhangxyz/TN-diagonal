gcc main.c -std=c99 -lm -O3
./a.out 2 2 | tee 2.2
./a.out 2 3 | tee 2.3
./a.out 2 4 | tee 2.4
./a.out 2 5 | tee 2.5
./a.out 2 6 | tee 2.6
./a.out 3 3 | tee 3.3
./a.out 3 4 | tee 3.4

