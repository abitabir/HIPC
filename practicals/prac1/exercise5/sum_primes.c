#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>

bool is_prime(int o) {  // checks if input is prime
    int p = 2;
    bool isprime = true;  // true until contradiction
    while (isprime == true && p <= o/2) {
        if (o % p == 0) isprime = false;
        p++;
    }
    return isprime;
}

int sum_primes(int m) {  // sums the primes up to input positive integer m
    int i, sum = 0;

    for (i = 2; i <= m; i++) {
        if (is_prime(i)) {
            sum += i;
        }
    }
    return sum;
}
