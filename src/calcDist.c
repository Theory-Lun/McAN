//
// Created by lun on 3/9/21.
//

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#include "edges.h"
#include "haplotype.h"
#include "virus.h"

bool isSameMut(struct Mutation * mut1, struct Mutation * mut2);
int cmp_time_virus(const void* _a, const void* _b);
int cmp_time_date(const void* _a, const void* _b);
void edges_show_new(struct Edge * edges, int numOfEdges);
void reverse_direct(struct Edge * edge);

int findset(int r, int *pr);

#ifdef MUTILTHREADING
void calcDist_thread(struct Hap* haps, int numOfHaps, struct Edge* edges, int nTask, int nThread);
#endif

//#define DEBUG
#define ONLY_SUBSET
void calcMintime(struct Hap * pHap1, struct Hap * pHap2, struct Edge * edge){
#ifdef DEBUG
    printf("\n");
    printf("%d.%d.%d\n", pHap1->date->year, pHap1->date->month, pHap1->date->day);
    printf("%d.%d.%d\n", pHap2->date->year, pHap2->date->month, pHap2->date->day);
#endif

    int diff = cmp_time_date(pHap1->date, pHap2->date);
    if(diff < 0){
        edge->mintime = SOURCE_TIME;
    }
    else if( diff > 0 ){
        edge->mintime = TARGET_TIME;
    }
    else{
        edge->mintime = EQ_TIME;
    }
#ifdef DEBUG
    //printf("mintime = %d\n", edge->mintime);
#endif

}

void calcSubset_linear_forsorted(struct Hap * pHap_source, struct Hap * pHap_target, struct Edge * edge){
    int l1 = pHap_source->numOfMut;
    int l2 = pHap_target->numOfMut;

    if(l1 == l2 || l1 > l2){
        edge->subset = NOSUBSET;
        return;
    }
    if(0 == l1){
        edge->subset = SOURCEISSUBSET;
        return;
    }
    if(0 == l2){
        edge->subset = TARGETISSUBSET;
        return;
    }

    int j_old=0;
    for (int i = 0; i < l1; ++i) {  //i, index of mut in source
        bool contain = false;
        for (int j = j_old; j < l2; ++j) { //j, index of mut in target
            if( true == isSameMut( &((pHap_source->mutations)[i]), &((pHap_target->mutations)[j]) ) ){
                contain = true;
                j_old = j+1;
                break;
            }
        }
        if( false == contain){
            edge->subset = NOSUBSET;
            return;
        }
    }
    edge->subset = SOURCEISSUBSET;
}

int hap_forAnc_cmp(const void* _a, const void* _b){
    if(NULL == _a || NULL == _b){
        fprintf(stderr, "no time infor\n");
        exit(1);
    }
    struct Hap * a = (struct Hap*) _a;
    struct Hap * b = (struct Hap*) _b;
    int diff = b->numOfMut - a->numOfMut; //length of mut, from large to small
    if(0 != diff){
        return diff;
    }
    else{
        int diff_size = b->numOfVirus - a->numOfVirus;//large to small
        if(0 != diff_size){
            return diff_size;
        }
        else{
#ifdef TEST_REVERSE_CHRONOLOGICAL_ORDER
            int diff_time = cmp_time_date(b->date, a->date);//late to early
#else
            int diff_time = cmp_time_date(a->date, b->date);//early to late
#endif
            return diff_time;
        }
    }
}

#ifdef MUTILTHREADING
void calcDistDirected_sorthap(struct Hap * haps, int numOfHaps, struct Edge * edges, int * numOfEdges, int nThread) {
#else
void calcDistDirected_sorthap(struct Hap * haps, int numOfHaps, struct Edge * edges, int * numOfEdges ){
#endif
    if(NULL == edges){
        edges = (struct Edge * )malloc(sizeof(struct Edge) * *numOfEdges);
    }
    if(NULL == haps){
        printf("error in calcDistAll!\n");
        exit(1);
    }

    //init, all haplotype point to itself
    for (int i = 0; i < *numOfEdges; ++i) {
        edges[i].source = i;
        edges[i].target = i;
        edges[i].pHap1 = &(haps[i]);
        edges[i].pHap2 = &(haps[i]);
        edges[i].distance = 0;
    }

    qsort(haps, numOfHaps, sizeof(struct Hap), hap_forAnc_cmp);//good!!


#ifdef MUTILTHREADING
    int nTask = numOfHaps;
    if (nThread > nTask)
        nThread = nTask;

    if (nThread > 1) {
        calcDist_thread(haps, numOfHaps, edges, nTask, nThread);
    }
    else {
#endif
        for (int i = 0; i < numOfHaps; ++i) {//target, for each target
            struct Edge *tmpE = NULL;
            tmpE = &(edges[i]);
            tmpE->target = i;
            tmpE->pHap2 = &(haps[tmpE->target]);
            for (int j = i + 1; j < numOfHaps; ++j) {//source, determine if j is i's ancestor
                tmpE->source = j;
                tmpE->pHap1 = &(haps[tmpE->source]);
                calcSubset_linear_forsorted(tmpE->pHap1, tmpE->pHap2, tmpE);//2021.3.28
                if( NOSUBSET == tmpE->subset ){//only keep subset
                    continue;
                }
                else{
                    tmpE->distance = tmpE->pHap2->numOfMut - tmpE->pHap1->numOfMut;
                    calcMintime(tmpE->pHap1, tmpE->pHap2, tmpE );
                    break;//find ancestor for next hap
                }
            }
            printf("Percent completed:%5.1f%%\r", (float )i / (float)numOfHaps * 100.0);
            fflush(stdout);
        }
#ifdef MUTILTHREADING
    }
#endif
#ifdef DEBUG
    //edges_show_new(edges, *numOfEdges);
#endif
}
