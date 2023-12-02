
#include "./quad_tree.h"
#include "./line.h"
#include <assert.h>
#include "./vec.h"
#include "./intersection_detection.h"
#include "./parallel.h"
#include <stdlib.h>
#include <string.h>

static void Quadtree_build_dfs(Line** lines, unsigned int numLines, double timestep, QuadtreeNode* node){
    assert(numLines > 0);
    if(numLines < MAX_LINES_PER_NODE){
        node->lines = lines;
        node->linesSize = numLines;
        node->children[0] = NULL;
        node->children[1] = NULL;
        node->children[2] = NULL;
        node->children[3] = NULL;
        return;
    }
    unsigned int* belong_sub[4];
    unsigned int* visited = malloc(sizeof(unsigned int) * numLines);
    unsigned int* belong_this = malloc(sizeof(unsigned int) * numLines);
    memset(visited, 0, sizeof(unsigned int) * numLines);
    memset(belong_this, 0, sizeof(unsigned int) * numLines);
    for(int i = 0; i < 4; i++) {
        belong_sub[i] = malloc(sizeof(unsigned int) * numLines);
        memset(belong_sub[i], 0, sizeof(unsigned int) * numLines);
    }

    Vec** children_min_max = Node_split(node);
    for(int i = 0; i < 4; i++){
        Vec min = children_min_max[i][0];
        Vec max = children_min_max[i][1];
        for(int index = 0; index < numLines; index++){
            if(visited[index] == 1) continue;
            Line* line = lines[index];

            // if(i == 0){
            //     assert(line->p1.x >= node->min.x && line->p1.x <= node->max.x);
            //     assert(line->p1.y >= node->min.y && line->p1.y <= node->max.y);
            //     assert(line->p2.x >= node->min.x && line->p2.x <= node->max.x);
            //     assert(line->p2.y >= node->min.y && line->p2.y <= node->max.y);
            //     assert(line->p1.x + line->velocity.x * timestep >= node->min.x && line->p1.x + line->velocity.x * timestep <= node->max.x);
            //     assert(line->p1.y + line->velocity.y * timestep >= node->min.y && line->p1.y + line->velocity.y * timestep <= node->max.y);
            //     assert(line->p2.x + line->velocity.x * timestep >= node->min.x && line->p2.x + line->velocity.x * timestep <= node->max.x);
            //     assert(line->p2.y + line->velocity.y * timestep >= node->min.y && line->p2.y + line->velocity.y * timestep <= node->max.y);
            // }

            ParallelogramIntersectionType type = intersectParallelograms(line->p1, line->p2, 
                                    Vec_add(line->p1, Vec_multiply(line->velocity, timestep)),
                                    Vec_add(line->p2, Vec_multiply(line->velocity, timestep)), 
                                    min, Vec_make(max.x, min.y), 
                                    Vec_make(min.x, max.y), max);
            if(type == NO_INTERSECTION_PARALLELOGRAM){
                ;;
            }
            else if(type == P1_IN_P2){
                belong_sub[i][index] = 1;
                visited[index] = 1;
            }
            else if(type == INTERSECTION){
                belong_this[index] = 1;
                visited[index] = 1;
            }       
        }
    }
    for(int i = 0; i < numLines; ++i){
        if(visited[i] == 0){
            // printf("error: line %d not visited\n", i);
            // printf("line %d: (%lf, %lf) -> (%lf, %lf)\n", i, lines[i]->p1.x, lines[i]->p1.y, lines[i]->p2.x, lines[i]->p2.y);
            // printf("line %d: (%lf, %lf) -> (%lf, %lf)\n", i, lines[i]->p1.x + lines[i]->velocity.x * timestep, lines[i]->p1.y + lines[i]->velocity.y * timestep, lines[i]->p2.x + lines[i]->velocity.x * timestep, lines[i]->p2.y + lines[i]->velocity.y * timestep);
            belong_this[i] = 1;
        }
    }
    unsigned int* belong_sub_dst_pos[4];
    unsigned int* belong_this_dst_pos = malloc(sizeof(unsigned int) * numLines);

    for(int i = 0; i < 4; i++) {
        belong_sub_dst_pos[i] = malloc(sizeof(unsigned int) * numLines);
    }
    unsigned int lines_size[4];
    for(int i = 0; i < 4; i++){
        lines_size[i] = scan_serial(belong_sub[i], numLines, belong_sub_dst_pos[i]);
    }
    // 初始化各个Line*数组
    Line** lines_sub[4];
    for(int i = 0; i < 4; i++){
        lines_sub[i] = malloc(sizeof(Line*) * lines_size[i]);
    }
    for(int i = 0; i < 4; i++){
        for(int index = 0; index < numLines; index++){
            if(belong_sub[i][index]){
                lines_sub[i][belong_sub_dst_pos[i][index]] = lines[index];
            }
        }
    }
    unsigned int lines_size_this = scan_serial(belong_this, numLines, belong_this_dst_pos);
    Line** line_this = malloc(sizeof(Line*) * lines_size_this);
    for(int index = 0; index < numLines; index++){
        if(belong_this[index]){
            line_this[belong_this_dst_pos[index]] = lines[index];
        }
    }
    assert(lines_size_this + lines_size[0] + lines_size[1] + lines_size[2] + lines_size[3] == numLines);
    
    node->lines = line_this;
    node->linesSize = lines_size_this;

    for (int i = 0; i < 4; i++) {
        if(lines_size[i] == 0){
            node->children[i] = NULL;
            continue;
        }
        node->children[i] = malloc(sizeof(QuadtreeNode));
        node->children[i]->min = children_min_max[i][0];
        node->children[i]->max = children_min_max[i][1];
        Quadtree_build_dfs(lines_sub[i], lines_size[i], timestep, node->children[i]);
    }
}
static int compare_with_self(Line** group, const unsigned int size,
                            IntersectionEventList* intersectionEventList,
                            const double timeStep){
    int numLineLineCollisions = 0;
    for(int i = 0; i < size; i++){
        Line* l1 = group[i];
        for(int j = i + 1; j < size; j++){
            Line* l2 = group[j];
            // intersect expects compareLines(l1, l2) < 0 to be true.
            // Swap l1 and l2, if necessary.
            if (compareLines(l1, l2) >= 0) {
                Line *temp = l1;
                l1 = l2;
                l2 = temp;
            }

            IntersectionType intersectionType =
                intersect(l1, l2, timeStep);
            if (intersectionType != NO_INTERSECTION) {
                IntersectionEventList_appendNode(intersectionEventList, l1, l2, intersectionType);
                numLineLineCollisions++;
            }
        }
    }
    return numLineLineCollisions;
}

static int compare_with_other(Line** restrict group1, Line** restrict group2, 
                            const unsigned int size1, const unsigned int size2, 
                            IntersectionEventList* intersectionEventList,
                            const double timeStep){
    int numLineLineCollisions = 0;
    for(int i = 0; i < size1; i++){
        Line* l1 = group1[i];
        for(int j = 0; j < size2; j++){
            Line* l2 = group2[j];
            // intersect expects compareLines(l1, l2) < 0 to be true.
            // Swap l1 and l2, if necessary.
            if (compareLines(l1, l2) >= 0) {
                Line *temp = l1;
                l1 = l2;
                l2 = temp;
            }
            
            IntersectionType intersectionType =
                intersect(l1, l2, timeStep);
            if (intersectionType != NO_INTERSECTION) {
                IntersectionEventList_appendNode(intersectionEventList, l1, l2,
                                                intersectionType);
                numLineLineCollisions++;
            }
        }
    }
    return numLineLineCollisions;
}
static int num = 0;

static int QuadtreeNode_traverse_dfs(QuadtreeNode* node, IntersectionEventList* intersectionEventList, 
                                MyList* path, const double timestep, int depth){
    if(node == NULL){
        return 0;
    }
    int count = 0;
    Line** group = node->lines;
    unsigned int size = node->linesSize;
    MyList_push_back(path, group, size);

    // printf("depth: %d, box: (%lf, %lf) (%lf, %lf)\n", depth, node->min.x, node->min.y, node->max.x, node->max.y);
    // for(int i = 0; i < 4; ++i){
    //     assert(node->children[i] == NULL || (node->children[i]->min.x >= node->min.x && node->children[i]->min.y >= node->min.y));
    //     assert(node->children[i] == NULL || (node->children[i]->max.x <= node->max.x && node->children[i]->max.y <= node->max.y));
    //     if(node->children[i] != NULL)
    //         printf("child box: (%lf, %lf) (%lf, %lf)\n", node->children[i]->min.x, node->children[i]->min.y, node->children[i]->max.x, node->children[i]->max.y);
    // }
    

    for(int i = 0; i < 4; i++){
        count += QuadtreeNode_traverse_dfs(node->children[i], intersectionEventList, path, timestep, depth + 1);
    }
    MyList_pop_back(path);
    // num += size;

    // printf("size_all: %d\n", num);
    //     printf("size: %d\n : box: (%lf, %lf) (%lf, %lf)\n", size, node->min.x, node->min.y, node->max.x, node->max.y);
        for(int i = 0; i < size; i++){
            // printf(" (%lf, %lf) (%lf, %lf) - ", group[i]->p1.x, group[i]->p1.y, group[i]->p2.x, group[i]->p2.y);
            // printf("(%lf, %lf) (%lf, %lf) - ", group[i]->p1.x + group[i]->velocity.x * timestep, group[i]->p1.y + group[i]->velocity.y * timestep, 
            //                                     group[i]->p2.x + group[i]->velocity.x * timestep, group[i]->p2.y + group[i]->velocity.y * timestep);
            // if(i % 5 == 4) printf("\n");     

#define min(a, b) ((a) < (b) ? (a) : (b))
#define max(a, b) ((a) < (b) ? (b) : (a))
            assert(min(group[i]->p1.x, BOX_XMAX) <= node->max.x && max(group[i]->p1.x, BOX_XMIN) >= node->min.x);
            assert(min(group[i]->p1.y, BOX_YMAX) <= node->max.y && max(group[i]->p1.y, BOX_YMIN) >= node->min.y);
            assert(min(group[i]->p2.x, BOX_XMAX) <= node->max.x && max(group[i]->p2.x, BOX_XMIN) >= node->min.x);
            assert(min(group[i]->p2.y, BOX_YMAX) <= node->max.y && max(group[i]->p2.y, BOX_YMIN) >= node->min.y);

            assert(min(group[i]->p1.x + group[i]->velocity.x * timestep, BOX_XMAX) <= node->max.x && 
                    max(group[i]->p1.x + group[i]->velocity.x * timestep, BOX_XMIN) >= node->min.x);
            assert(min(group[i]->p1.y + group[i]->velocity.y * timestep, BOX_YMAX) <= node->max.y &&
                    max(group[i]->p1.y + group[i]->velocity.y * timestep, BOX_YMIN) >= node->min.y);
            assert(min(group[i]->p2.x + group[i]->velocity.x * timestep, BOX_XMAX) <= node->max.x &&
                    max(group[i]->p2.x + group[i]->velocity.x * timestep, BOX_XMIN) >= node->min.x);
            assert(min(group[i]->p2.y + group[i]->velocity.y * timestep, BOX_YMAX) <= node->max.y &&
                    max(group[i]->p2.y + group[i]->velocity.y * timestep, BOX_YMIN) >= node->min.y);

        }
    count += compare_with_self(group, size, intersectionEventList, timestep);
    assert(depth == path->size);
    MyListNode* p = path->head->next;

    while(p != path->tail){
        count += compare_with_other(group, p->data,
                                    size, p->data_size,
                                    intersectionEventList, timestep);
        if(count != 0)
            printf("count: %d\n", count);
        // printf("size1: %d, size2: %d\n", size, p->data_size);
        // printf("box1: (%lf, %lf) (%lf, %lf)\n", node->min.x, node->min.y, node->max.x, node->max.y);
        p = p->next;
    }
    return count;
}

//     printf("\n");

Vec** Node_split(QuadtreeNode* node){
    Vec parallelogram_min = node->min;
    Vec parallelogram_max = node->max;
    Vec** children = malloc(sizeof(Vec*) * 4);
    for(int i = 0; i < 4; i++){
        children[i] = malloc(sizeof(Vec) * 2);
    }
    children[0][0] = parallelogram_min;
    children[0][1] = Vec_make((parallelogram_min.x + parallelogram_max.x) / 2, (parallelogram_min.y + parallelogram_max.y) / 2);
    children[1][0] = Vec_make((parallelogram_min.x + parallelogram_max.x) / 2, parallelogram_min.y);
    children[1][1] = Vec_make(parallelogram_max.x, (parallelogram_min.y + parallelogram_max.y) / 2);
    children[2][0] = Vec_make(parallelogram_min.x, (parallelogram_min.y + parallelogram_max.y) / 2);
    children[2][1] = Vec_make((parallelogram_min.x + parallelogram_max.x) / 2, parallelogram_max.y);
    children[3][0] = Vec_make((parallelogram_min.x + parallelogram_max.x) / 2, (parallelogram_min.y + parallelogram_max.y) / 2);
    children[3][1] = parallelogram_max;
    return children;
}

void QuadtreeNode_destroy(QuadtreeNode* node){
    // if(node->linesSize > 0){
    //     free(node->lines);
    // }
    for(int i = 0; i < 4; i++){
        if(node->children[i] != NULL){
            QuadtreeNode_destroy(node->children[i]);
        }
    }
    free(node);
}



Quadtree* Quadtree_build(Line** lines, unsigned int numLines, double timestep){
   
    Quadtree* qt = malloc(sizeof(Quadtree));
    QuadtreeNode* node = malloc(sizeof(QuadtreeNode));
    qt->root = node;
    node->min = Vec_make(BOX_XMIN, BOX_YMIN);
    node->max = Vec_make(BOX_XMAX, BOX_YMAX);
    Quadtree_build_dfs(lines, numLines, timestep, node);
    return qt;
}

void Quadtree_destroy(Quadtree* qt){
    QuadtreeNode* node = qt->root;
    QuadtreeNode_destroy(node);
    free(qt);
}


unsigned int Quadtree_traverse(Quadtree* tree, IntersectionEventList* intersectionEventList, double timestep){
    QuadtreeNode* node = tree->root;

    MyList* list = MyList_create();
    unsigned int res = QuadtreeNode_traverse_dfs(node, intersectionEventList, list, timestep, 0);
    MyList_destroy(list);
    
    return res;
}



