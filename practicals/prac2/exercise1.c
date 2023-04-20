#include <stdio.h>
#include <string.h>
#include <stdbool.h>

int count_words(char *str) {  // this function returns the number of words in str
    // assumes all strings end with words
    int count_space = 0;
    bool has_space = true;
    char * space = " ";

    while (has_space) {
        count_space++;
        str = strstr(str, space);
        if (str==NULL) {
            has_space = false;
        } else {
            str = &(str[1]);
        }
    }

    return count_space;
}

int main() {
    char str[100];

    printf("Enter a string:");
    fgets(str, 100, stdin);    // this function reads a line or at most 99 bytes from stdin file stream that represents the keyboard

    printf("Number of words in the entered string is %d\n", count_words(str));

    return 0;

}