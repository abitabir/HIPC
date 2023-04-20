#include <stdio.h>
#include <stdlib.h>

// write the code for the sort() function that sorts an integer array in ascending order

void swap(int *a, int *b) {
    int tmp = *a;
    *a = *b;
    *b = tmp;
}

void sort(int *arr, int len) {
    int swapped = 1;
    int n = len;
    while (swapped == 1) {
        swapped = 0;
        int i;
        for (i = 1; i < len; i++) {
            if (arr[i-1] > arr[i]) {
                swap(&arr[i-1], &arr[i]);
                swapped = 1;
            }
        }
        n--;
    }
}

int main() {
    int i;

    int count = 0;
    printf("How many elements would you like? ");
    scanf("%d", &count);

    int *x = (int *) malloc(sizeof(int) * count);

    printf("Enter the values you'd like to store: ");
    for (i = 0; i < count; i++) {
        scanf("%d", &x[i]);
    }

    sort(x, count); // sort() function sorts the array x in ascending order

    printf("The sorted array is as follows:\n");

    for (i=0; i < count; i++){
        printf("%d ", x[i]);
    }

    printf("\n");
    return 0;
}
