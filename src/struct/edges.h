//
// Created by lun on 3/9/21.
//

#include "haplotype.h"

#ifndef HAPNET_V0_2_3_EDGES_H
#define HAPNET_V0_2_3_EDGES_H

#define NOSUBSET 0
#define SOURCEISSUBSET 1
#define TARGETISSUBSET 2

#define EQ_TIME 0
#define SOURCE_TIME 1
#define TARGET_TIME 2

struct Edge{
    int source;
    int target;
    int distance;
    struct Hap * pHap1;
    struct Hap * pHap2;
    short subset;//0: no subset relation, 1: source is subset of target, 2: target is subset of source
    short mintime;
};


#endif //HAPNET_V0_2_3_EDGES_H