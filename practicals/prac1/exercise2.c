#include <stdio.h>

// write code that swaps the values of two integers

int main(int argc, char *argv[]) {
   int val_a = 50;
   int val_b = 20;

   printf("val_a was %d (should be 50)\n", val_a);
   printf("val_b was %d (should be 20)\n", val_b);

   int* addr_a = &val_a;
   int* addr_b = &val_b;

   int temp_a = * addr_a;
   * addr_a = * addr_b;
   // * addr_a is (the value stored in (the memory address that addr_a is pointing to))
   // * addr_b is (the value stored in (the memory address that addr_b is pointing to))
   * addr_b = temp_a;
   
   printf("val_a is %d (should be 20)\n", val_a);
   printf("val_b is %d (should be 50)\n", val_b);

  return 0;
}