#include <stdio.h>
#include <stdlib.h>

struct element {
  struct element * next;
  int data;
};

struct queue {
  struct element * head;
};

void enqueue(struct queue * q, int val) {
  struct element * elem = malloc(sizeof(struct element));

  elem->data = val;
  elem->next = NULL; // Really important to explicitly set this to null. Malloc does not zero memory

  if (q->head == NULL) {
    // Empty list, we need to append to head
    q->head = elem;
  } else {
    // List has some elements, find the end and append to that
    struct element * tail = q->head;
    while (tail->next != NULL) {
      tail = tail->next;
    }
    tail->next = elem;
  }
}

int dequeue(struct queue * q) {
    struct element *elem = q->head;
    q->head = elem->next;
    return elem->data;
}

int isempty(struct queue * q) {
    if (q->head == NULL) return 1;
    return 0;
}

int main() {
    struct queue my_q = { NULL };

    enqueue(&my_q, 1);
    enqueue(&my_q, 2);
    enqueue(&my_q, 3);
    enqueue(&my_q, 4);
    enqueue(&my_q, 5);

    while (isempty(&my_q) == 0) {
        int num = dequeue(&my_q);
        printf("%d, ", num);
    }
    printf("\n");

}
