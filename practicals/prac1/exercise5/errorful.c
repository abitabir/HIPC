#include <stdlib.h>
#include <stdio.h>

int main(int argc, char *argv[]) {
   int *my_array = (int *) malloc(5 * sizeof(int)); 

   int i;

   for (i = 0; i < 1000000; i++) { 
      my_array[i] = i;
   }

   free(my_array);

   printf("done");
}