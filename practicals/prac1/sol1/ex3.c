#include <stdio.h>

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
    int x[]= { 4, 1, 4, 3, 10, 5 };
    int i;

    sort(x, 6); // sort() function sorts the array x in ascending order

    printf("The sorted array is as follows:\n");

    for (i=0; i < 6; i++){
        printf("%d ", x[i]);
    }

    printf("\n");
    return 0;
}
