//
// Created by lun on 3/12/21.
//
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <hashmap.h>
#include "virus.h"
#include "haplotype.h"
#include "edges.h"
//#define DEBUG

void countryDistribution(struct Hap * haps, int numOfHaps, struct Virus * viruses, int numOfViruses){
    for (int i = 0; i < numOfHaps; ++i) {
        int r;
        hashmap_init(&haps[i].country, hashmap_hash_string, strcmp);
        for (int j = 0; j < haps[i].numOfVirus; ++j) {
            int virusIndex = haps[i].allViruses[j];
            if(NULL == viruses[virusIndex].location){
                fprintf(stderr,"error: %s no location", viruses[virusIndex].acc);
                exit(1);
            }
            if(NULL != viruses[virusIndex].location->country){
                int * numOfCountry = (int *)malloc(sizeof(int));
                *numOfCountry = 1;
                r = hashmap_put(&haps[i].country, viruses[virusIndex].location->country, numOfCountry);
                if(r < 0){
                    free(numOfCountry);
                    numOfCountry = NULL;
                    numOfCountry = hashmap_get(&haps[i].country, viruses[virusIndex].location->country);
                    (*numOfCountry) ++;
                }
            }
        }
#ifdef DEBUG
        const char * tmpcountry;
        int *value;
        printf("\n%d\n", i);
        hashmap_foreach(tmpcountry, value, &haps[i].country) {
            printf("acc=%s, count = %d\n", tmpcountry, *value);
        }
#endif
    }
}