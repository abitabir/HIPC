#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>

bool isPrime(int o) {  // checks if input is prime
    int p = 2;
    bool isprime = true;  // true until contradiction
    while (isprime == true && p <= o/2) {
        if (o % p == 0) isprime = false;
        p++;
    }
    return isprime;
}

int primeFind(int m, int *primes_addr) {  // computes the primes up to given positive integer m
    int i = 2;
    int curr_count;

    for (i = 2; i <= m; i++) {
        if (isPrime(i)) {
            primes_addr[0]++;  // first index is of prime count i.e. array length
            curr_count = primes_addr[0];
            int *new_addr = (int *) realloc(primes_addr, sizeof(int) * primes_addr[0]);

            if (new_addr != NULL) primes_addr = new_addr;
            else perror("Failed to reallocate memory.");

            primes_addr[primes_addr[0]] = i;
        }
    }
    return 0;
}

int main() {
    int n;
    int i;
    printf("Enter a positive integer under 100:");
    scanf("%d", &n); // scanf() scans the integer entered by the user

    int *primes = (int *) malloc(sizeof(int) * 1);
    primes[0] = 0;
    primeFind(n, primes);  // arrays are already pointers so do not need to reference the address with &

    printf("The primes up to %d are:\n", n);
    for(i=1;i<=primes[0];i++) printf(" %d ", primes[i]);
    
    free(primes);  // wrapup
    primes = NULL;

    return 0;
}
