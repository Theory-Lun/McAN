//
// Created by lun on 3/19/21.
//

#include "edges.h"
#include "haplotype.h"

void treeHeight(struct Hap * haps, int numOfHaps, int * height){
    *height = 1;
    int jump = 0;
    for (int i = 0; i < numOfHaps; ++i) {
        jump = haps[i].jump;
        if( *height < jump){
            *height = jump;
        }
    }
}

void calcDisAll(int *network, int numOfEdgesInNet, struct Edge *edges, int *disAll) {
    int sum = 0;
    for (int i = 0; i < numOfEdgesInNet; ++i) {
        int edgeIndex = network[i];
        sum += edges[edgeIndex].distance;
    }
    *disAll = sum;
}

void calcSrcSizeAve(int *network, int numOfEdgesInNet, struct Edge *edges, double *aveOfSourceSize) {
    double sum = 0;
    for (int i = 0; i < numOfEdgesInNet; ++i) {
        int edgeIndex = network[i];
        sum += ( (double)edges[edgeIndex].pHap1->numOfVirus / (double)numOfEdgesInNet );
    }
    *aveOfSourceSize = sum;
}

void calcNumOfEdgesTimeGood(int *network, int numOfEdgesInNet, struct Edge *edges, int *numOfEdgesTimeGood){
    int sum = 0;
    for (int i = 0; i < numOfEdgesInNet; ++i) {
        int edgeIndex = network[i];
        if( SOURCEISSUBSET == edges[edgeIndex].mintime){
            sum ++;
        }
    }
    *numOfEdgesTimeGood = sum;
}