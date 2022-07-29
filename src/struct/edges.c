//
// Created by lun on 3/9/21.
//


#include <stdio.h>
#include "edges.h"

void show_mutation(struct Mutation *pMut, int numOfIndel);
void show_date(struct Date * date);
//#define DEBUG
void edges_show(struct Edge * edges, int numOfEdges){
    if(NULL == edges){
        return;
    }
    printf("edges:\n");
    for (int i = 0; i < numOfEdges; ++i) {
        printf("\n");
        printf("Hap1: ");
        show_mutation(edges[i].pHap1->mutations, edges[i].pHap1->numOfMut );
        printf("Hap2: ");
        show_mutation(edges[i].pHap2->mutations, edges[i].pHap2->numOfMut );
        printf("source: %d\ttarget: %d\tdist: %d\tsubset %d\n", edges[i].source, edges[i].target, edges[i].distance, edges[i].subset);

#ifdef DEBUG
        if(SOURCEISSUBSET == edges[i].subset && edges[i].pHap1->numOfMut >10){
            getchar();
        }
#endif
    }
}

void edges_show_new(struct Edge * edges, int numOfEdges){
    if(NULL == edges){
        return;
    }
    printf("edges:\n");
    for (int i = 0; i < numOfEdges; ++i) {
        printf("\n");
        printf("Hap1: ");
        show_mutation(edges[i].pHap1->mutations, edges[i].pHap1->numOfMut );
        printf("Hap2: ");
        show_mutation(edges[i].pHap2->mutations, edges[i].pHap2->numOfMut );
        printf("mintime = %d\n", edges[i].mintime);
        printf("source: %d\ttarget: %d\tdist: %d\tsubset %d\n", edges[i].source, edges[i].target, edges[i].distance, edges[i].subset);
        printf("source: %d\tsize:%d\n", edges[i].source, edges[i].pHap1->numOfVirus);
        show_date((edges[i].pHap1)->date);
        printf("target: %d\tsize:%d\n", edges[i].target, edges[i].pHap2->numOfVirus);
        show_date((edges[i].pHap2)->date);
    }
}