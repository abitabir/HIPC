#include <stdio.h>
#include <stdlib.h>

struct element {
  struct element * next;
  int data;
};

struct linked_list {
  struct element * head;
};

void append_int(struct linked_list * list, int val) {
  struct element * elem = malloc(sizeof(struct element));

  elem->data = val;
  elem->next = NULL; // Really important to explicitly set this to null. Malloc does not zero memory

  if (list->head == NULL) {
    // Empty list, we need to append to head
    list->head = elem;
  } else {
    // List has some elements, find the end and append to that
    struct element * tail = list->head;
    while (tail->next != NULL) {
      tail = tail->next;
    }
    tail->next = elem;
  }
}

void prepend_int(struct linked_list * list, int val) {
    struct element * elem = malloc(sizeof(struct element));

    elem->data = val;
    elem->next = list->head;
    list->head = elem;
}

int remove_from_head(struct linked_list * list) {
    struct element *elem = list->head;
    list->head = elem->next;
    return elem->data;
}

int main() {
    struct linked_list my_list = { NULL };

    append_int(&my_list, 1);
    append_int(&my_list, 2);
    append_int(&my_list, 3);
    append_int(&my_list, 4);
    append_int(&my_list, 5);

    int ret = remove_from_head(&my_list); // will remove first element from list (1)
    printf("Removed %d\n", ret);

    prepend_int(&my_list, 6); // will add 6 to start of list

    struct element *my_element;
    for (my_element = my_list.head; my_element != NULL; my_element = my_element->next) {
        printf("%d, ", my_element->data);
    }
    printf("\n");

}
