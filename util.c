#include<stdio.h>
#include<stdlib.h>
#include "util.h"



// LIST
void list_init(list* list) {
	list->prev = list->next = list;
}

void list_append(list* list, node_t* item) {
	item->prev = list->prev;
	item->next = list;
	item->prev->next = item;
	item->next->prev = item;
}

void list_erase(list* list, node_t* item) {

	if (list != item) {
		item->next->prev = item->prev;
		item->prev->next = item->next;
	}
}

void list_empty(list* list) {
	node_t *li, *l;
	li = list->next;
	for (; li != list;) {
		l = li;
		li = li->next;
		free(l);
	}
	list_init(list);
}

void list_destroy(list* list) {
	node_t *li;
	list->prev->next = NULL;
	for (; list;) {
		li = list;
		list = li->next;
		free(li);
	}
}


// GRAPH
void graph_init(Graph* graph) {
	//graph->points = (Point*)malloc(sizeof(Point));
	points_init(graph->points);
	//graph->edges = (Edge*)malloc(sizeof(Edge));
	edges_init(graph->edges);
}

void graph_destroy(Graph* graph) {
	points_destroy(graph->points);
	edges_destroy(graph->edges);
}

void graph_refresh(Graph* graph) {
	Point* p;
	for (p = graph->points->next; p != graph->points; p = p->next) {
		p->edges_num = 0;
		list_empty(graph->edges);
	}
}

void points_init(Point* points) {
	list_init((list*)points);
}

void edges_init(Edge* edges) {
	list_init((list*)edges);
}

void points_destroy(Point* points) {
	list_destroy((list*)points);
}

void edges_destroy(Edge* edges) {
	list_destroy((list*)edges);
}

//int point_not_in(Point* point, int id) {
//	Point* p;
//	for (p = point; p->next; p = p->next) {
//		if (p->id == id)
//			return p->index;
//	}
//	return -1;
//}