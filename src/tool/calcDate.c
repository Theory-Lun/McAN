//
// Created by lun on 2021/5/26.
//
#include "virus.h"

void NextMonth(short *year, short *month){
    (*month) ++;
    if(*month > 12){
        *month = 1;
        (*year) ++;
    }
}

void NextMonth_struct(struct Date *structDate){
    NextMonth(&(structDate->year), &(structDate->month));
}