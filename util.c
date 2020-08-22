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

void list_delete(list* list, node_t* item) {

	if (list != item) {
		item->next->prev = item->prev;
		item->prev->next = item->next;
	}
}

void list_erase(list* list, node_t* item) {

	if (list != item) {
		item->next->prev = item->prev;
		item->prev->next = item->next;
	}
	free(item);
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
	graph->points_num = 0;
	//graph->points = (Point*)malloc(sizeof(Point));
	points_init((Point*)graph->points);
	//graph->edges = (Edge*)malloc(sizeof(Edge));
	edges_init((Edge*)graph->edges);
}

void graph_destroy(Graph* graph) {
	points_destroy((Point*)graph->points);
	edges_destroy((Edge*)graph->edges);
}

void graph_refresh(Graph* graph) {
	Point* p;
	for (p = (Point*)graph->points->next; p != graph->points; p = (Point*)p->next) {
		p->in_num = 0;
		p->out_num = 0;
		list_init((list*)&p->in_edges);
		list_init((list*)&p->out_edges);
	}
	list_empty((list*)graph->edges);
}

void points_init(Point* points) {
	list_init((list*)points);
}

void point_append(Graph* graph, Point* point) {
	list_append((list*)graph->points, (node_t*)point);
	list_init((list*)&point->in_edges);
	list_init((list*)&point->out_edges);
	point->index = graph->points_num;
	graph->points_num++;
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

void print_graph(Graph* graph){
	Point* p;
	EdgeList* el;
	Edge* e;
	for (p = (Point*)graph->points->next; p != (Point*)graph->points; p = (Point*)p->next){
		printf("Point %d has %d ins and %d outs. they are:\n", p->index, p->in_num, p->out_num);
		printf("IN\n");
		for (el = (EdgeList*)p->in_edges.next; el != &p->in_edges; el = (EdgeList*)el->next){
			e = el->edge;
			printf("endpoints %d %d\n",e->endpoints[0]->index,e->endpoints[1]->index);
		}
		printf("OUT\n");
		for (el = (EdgeList*)p->out_edges.next; el != &p->out_edges; el = (EdgeList*)el->next){
			e = el->edge;
			printf("endpoints %d %d\n",e->endpoints[0]->index,e->endpoints[1]->index);
		}
	}
}