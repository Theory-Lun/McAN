//
// Created by lun on 2021/3/3.
//
#include "mutation.h"
#include "indel.h"
#ifndef HAPNET_V0_2_VIRUS_H
#define HAPNET_V0_2_VIRUS_H

struct Date{
    short int year; //2019, 2020, 2021, ...
    short int month; //1, 2, ..., 11, 12
    short int day; //1, 2, 3,..., 30, 31
};

struct Location{
    char * country;
    char * state;
    char * city;
};

struct Virus{
    //the char * need end with '\0'

    char * name;//virus strain name
    char * acc;//accession ID
    struct Location * location;//location
    struct Date * date; // sample collection date
    struct Mutation * mutations;
    int numOfMut;
    int indexOfHap;
};

#endif //HAPNET_V0_2_VIRUS_H
