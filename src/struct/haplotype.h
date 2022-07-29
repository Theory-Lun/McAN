//
// Created by lun on 2021/3/8.
//
#include <stdbool.h>
#include <hashmap.h>
#include "mutation.h"
#include "indel.h"

#ifndef HAPNET_V0_2_2_HAPLOTYPE_H
#define HAPNET_V0_2_2_HAPLOTYPE_H

struct Hap{
    int index;
    int numOfVirus;
    int * allViruses;
    int numOfMut;
    struct Mutation * mutations;
    HASHMAP(char, int) country;
    int jump;
    struct Date * date;//earliest appearance time
};

#endif //HAPNET_V0_2_2_HAPLOTYPE_H