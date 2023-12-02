
#include "./line.h"

// serial implementation of scan
unsigned int scan_serial(const unsigned int* restrict vec, unsigned int size, unsigned int* restrict result);

struct MyListNode {
    Line** data;
    unsigned int data_size;
    // unsigned int layer;
    struct MyListNode* next;
    struct MyListNode* prev;
};
typedef struct MyListNode MyListNode;
struct MyList {
    MyListNode* head;
    MyListNode* tail;
    unsigned int size;
};
typedef struct MyList MyList;

MyList* MyList_create();

// insert data to the end of the list
void MyList_push_back(MyList* list, Line** data, const unsigned int data_size);
// remove the last element of the list
void MyList_pop_back(MyList* list);
void MyList_destroy(MyList* list);