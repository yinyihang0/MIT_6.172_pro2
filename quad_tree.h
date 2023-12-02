#include "./line.h"
#include "./intersection_event_list.h"

#define MAX_LINES_PER_NODE 10

struct QuadtreeNode {
    Line** lines;
    unsigned int linesSize;
    struct QuadtreeNode* children[4];
    Vec min;
    Vec max;
};
typedef struct QuadtreeNode QuadtreeNode;

void QuadtreeNode_destroy(QuadtreeNode* node);

// split the area into 4 sub areas, setting children's min and max
Vec** Node_split(QuadtreeNode* node);

struct Quadtree {
    QuadtreeNode* root;
};
typedef struct Quadtree Quadtree;
Quadtree* Quadtree_build(Line** lines, unsigned int numLines, double timestep);

void Quadtree_destroy(Quadtree* qt);

unsigned int Quadtree_traverse(Quadtree* tree, IntersectionEventList* intersectionEventList, double timestep);

