//
// Created by lun on 3/19/21.
//

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>


#include "virus.h"
#include "haplotype.h"
#include "edges.h"

int cmp_time_virus(const void* _a, const void* _b){
    if(NULL == _a || NULL == _b){
        fprintf(stderr, "no time infor\n");
        exit(1);
    }
    struct Virus* a = (struct Virus*) _a;
    struct Virus* b = (struct Virus*) _b;
    int diff;


    diff = a->date->year - b->date->year;
    if(0 != diff){
        return diff;
    }

    if(-1 == a->date->month ){
        return 1;
    }
    if(-1 == b->date->month ){
        return -1;
    }
    diff = a->date->month - b->date->month;
    if(0 != diff){
        return diff;
    }


    if(-1 == a->date->day ){
        return 1;
    }
    if(-1 == b->date->day ){
        return -1;
    }
    diff = a->date->day - b->date->day;
    return diff;
}

int cmp_time_date(const void* _a, const void* _b){
    if(NULL == _a || NULL == _b){
        fprintf(stderr, "no time infor\n");
        exit(1);
    }
    struct Date* a = (struct Date*) _a;
    struct Date* b = (struct Date*) _b;
    int diff;


    diff = a->year - b->year;
    if(0 != diff){
        return diff;
    }


    if(-1 == a->month && -1 != b->month){
        return 1;
    }
    if(-1 == b->month && -1 != a->month){
        return -1;
    }

    diff = a->month - b->month;
    if(0 != diff){
        return diff;
    }

    if(a->day == b->day){
        return 0;
    }
    if(-1 == a->day ){
        return 1;
    }
    if(-1 == b->day ){
        return -1;
    }
    diff = a->day - b->day;
    return diff;
}