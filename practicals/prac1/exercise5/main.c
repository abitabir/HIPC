#include <stdio.h>

int sum(int);
int sum_squares(int);
int sum_primes(int);

int main() {
    int n;

    printf("Enter a positive integer:");
    scanf("%d", &n);

    printf("Sum of first %d positive integers is %d\n", n, sum(n));
    printf("Sum of squares of first %d positive integers is %d\n", n, sum_squares(n));
    printf("Sum of all prime numbers from the first %d positive integers is %d\n", n, sum_primes(n));

    return 0;
}