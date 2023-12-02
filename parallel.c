#include "./parallel.h"
// serial implementation of scan
unsigned int scan_serial(const unsigned int* restrict vec, unsigned int size, unsigned int* restrict result){
    result[0] = 0;
    for(int i = 1; i < size; i++){
        result[i] = result[i - 1] + vec[i - 1];
    }
    unsigned int sum = result[size - 1] + vec[size - 1];
    return sum;
}

MyList* MyList_create(){
    MyList* list = malloc(sizeof(MyList));
    list->head = malloc(sizeof(MyListNode));
    list->tail = malloc(sizeof(MyListNode));
    list->size = 0;

    list->head->next = list->tail;
    list->head->prev = NULL;
    list->tail->next = NULL;
    list->tail->prev = list->head;
    return list;
}
void MyList_push_back(MyList* list, Line** data, const unsigned int data_size){
    MyListNode* node = malloc(sizeof(MyListNode));
    node->data = data;
    node->data_size = data_size;
    node->next = list->tail;
    node->prev = list->tail->prev;
    list->tail->prev->next = node;
    list->tail->prev = node;
    list->size++;
}
void MyList_pop_back(MyList* list){
    MyListNode* node = list->tail->prev;
    node->prev->next = list->tail;
    list->tail->prev = node->prev;
    free(node);
    list->size--;
}
void MyList_destroy(MyList* list){
    MyListNode* node = list->head;
    while(node != NULL){
        MyListNode* next = node->next;
        free(node);
        node = next;
    }
    free(list);
}
