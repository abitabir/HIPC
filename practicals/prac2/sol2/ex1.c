#include <stdio.h>
#include <string.h>

int count_words(char *str) {
    // this function should return the number of words in str
    int words = 0;
    int i;
    char *my_str_ptr;
    // this loop starts from the beginning of the string and uses strstr to jump to
    // the next space character each time.
    // Words are only counted if the character straight after the space is not a space
    // (otherwise the loop continues for another iteration to compress the spaces)
    for (my_str_ptr = str; my_str_ptr != NULL; my_str_ptr = strstr(my_str_ptr, " ")) {
        my_str_ptr++;
        if (my_str_ptr[0] != ' ') words++;
    }
    return words;
}

int main() {
    char str[100];

    printf("Enter a string:");
    fgets(str, 100, stdin);    // this function reads a line or at most 99 bytes from stdin file stream that represents the keyboard

    printf("Number of words in the entered string is %d\n", count_words(str));

    return 0;
}
