#include <stdio.h>

// write the code for the function swap() that would swap the values of two integers

void swap(int *a, int *b) {
    int tmp = *a;
    *a = *b;
    *b = tmp;
}

int main(int argc, char *argv[]) {
   int val_a = 50;
   int val_b = 20;

   swap(&val_a, &val_b);
   // call the swap() function correctly, swap() should swap the values stored in val_a and val_b so
   // that the swapped values of val_a and val_b are printed in the next two lines

   printf("val_a is %d (should be 20)\n", val_a);
   printf("val_b is %d (should be 50)\n", val_b);

  return 0;
}
