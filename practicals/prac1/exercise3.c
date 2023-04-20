#include <stdio.h>
#include <stdbool.h>

// write the code for the sort() function that sorts an integer array in ascending order

void sort(int n, int arr[]) {
   int i;
   int tmp;
   bool swappings = true;
   while (swappings) {
      swappings = false;
      for (i = 0; i < n - 1; i++) {
         if (arr[i] > arr[i + 1]) {
            swappings = true;  // at least one swappings
            tmp = arr[i];
            arr[i] = arr[i + 1];
            arr[i + 1] = tmp;
         }
      }
   }
}

int main() {
   int x[]= { 4, 1, 4, 3, 10, 5 };
   int i;

   sort(6, x); // sort() function sorts the array x in ascending order

   printf("The sorted array is as follows:\n");

   for (i=0; i < 6; i++){
      printf("%d ", x[i]);
   }

   printf("\n");
   return 0;
}