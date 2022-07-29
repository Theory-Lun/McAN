//
// Created by lun on 5/18/21.
//
#include <stdbool.h>
#include <hashmap.h>
#include "virus.h"

#ifndef HAPNET_V0_2_5_HASHMAP_ACC2LINE_H
#define HAPNET_V0_2_5_HASHMAP_ACC2LINE_H
struct Hashmap_acc2line{
    HASHMAP(char, char) acc2line;
};

struct Hashmap_acc2Date{
    HASHMAP(char, struct Date) acc2date;
};

#endif //HAPNET_V0_2_5_HASHMAP_ACC2LINE_H
