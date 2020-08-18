#pragma once

// List
#define LIST_STRUCT \
	struct node* prev, * next

typedef struct node {
	LIST_STRUCT;
} node_t, list;

/* Initialize a list */
void list_init(list* list);
/* Add item to the back of the list */
void list_append(list* list, node_t* item);
/* Remove item from the list */
void list_erase(list* list, node_t* item);
/* Empty a list */
void list_empty(list* list);
/* Destroy a list */
void list_destroy(list* list);

// Graph (directed graph)
typedef struct Point Point;
typedef struct Edge Edge;
#define MAXEDGES  40

#define POINT_STRUCT \
	LIST_STRUCT; \
	unsigned char type; \
	Edge *edges[MAXEDGES]; \
	int	edges_num

#define EDGE_STRUCT \
	LIST_STRUCT; \
	Point *endpoints[2]

#define GRAPH_STRUCT \
	Point* points; \
	Edge* edges

struct Point {
	POINT_STRUCT;
};

struct Edge {
	EDGE_STRUCT;
};

typedef struct Graph {
	GRAPH_STRUCT;

}Graph;


void points_init(Point* points);
void edges_init(Edge* edges);
void graph_init(Graph* graph);
void points_destroy(Point* points);
void edges_destroy(Edge* edges);
void graph_destroy(Graph* graph);
void graph_refresh(Graph* graph);
