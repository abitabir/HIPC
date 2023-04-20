#include <stdio.h>

int is_prime(int n) {
    int i;
    for (i = 2; i < n; i++) {
        if (n % i == 0)
            return 0;
    }
    return 1;
}

int prime(int n) {
    int sum_of_primes = 0;

    int i;
    for (i = 2; i <= n; i++) {
        if (is_prime(i))
            sum_of_primes += i;
    }
    return sum_of_primes;
}
