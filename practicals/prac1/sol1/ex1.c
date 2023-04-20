#include <stdio.h>

int is_prime(int n) {
    int i;
    for (i = 2; i < n; i++) {
        if (n % i == 0)
            return 0;
    }
    return 1;
}

int main(int argc, char *argv[]) {
    int n;
    printf("Enter n: ");
    scanf("%d", &n);

    int i;
    for (i = 2; i <= n; i++) {
        if (is_prime(i))
            printf("%d, ", i);
    }
    printf("\n");
    return 0;
}
