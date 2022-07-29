//
// Created by lun on 2021/3/8.
//

#include <stdio.h>
#include <stdlib.h>

#include "haplotype.h"

void showMut(struct Mutation * mut, int numOfMut);

void hap_init(struct Hap * pHaps, int numOfHaps){
    if(NULL == pHaps){
        return;
    }
    else{
        for (int i = 0; i < numOfHaps; ++i) {
            pHaps[i].allViruses = NULL;
            pHaps[i].index = -1;
            pHaps[i].mutations = NULL;
            pHaps[i].numOfMut = -1;
            pHaps[i].numOfVirus = -1;
            pHaps[i].jump = -1;
            pHaps[i].date = NULL;
        }
    }
}

void hap_show(struct Hap * pHaps, int numOfHaps){
    if(NULL == pHaps){
        return;
    }
    else{
        printf("Haplotype:\n");
        for (int i = 0; i < numOfHaps; ++i) {
            printf("\n%d\n", i);
            printf("index = %d\n", pHaps[i].index);

            printf("numOfMut = %d\n", pHaps[i].numOfMut);
            showMut(pHaps[i].mutations, pHaps[i].numOfMut);

            printf("numOfVirus = %d\n", pHaps[i].numOfVirus);
            if(NULL != pHaps[i].allViruses){
                printf("viruses: \n");
                for (int j = 0; j < pHaps[i].numOfVirus; ++j) {
                    printf("%d\t", pHaps[i].allViruses[j]);
                }
            }
            printf("\n");
        }
    }
}

void free_haps2(struct Hap ** pHaps, int numOfHaps){
    if(NULL == *pHaps){
        return;
    }

    //free allViruses date mutations
    for (int i = 0; i < numOfHaps; ++i) {
        struct Hap *phaptmp = &((*pHaps)[i]);

        if(NULL != phaptmp->allViruses){
            free(phaptmp->allViruses);
            phaptmp->allViruses = NULL;
        }

        if(NULL != phaptmp->mutations){
            free(phaptmp->mutations);
            phaptmp->mutations = NULL;
        }

        if(NULL != phaptmp->date){
            free(phaptmp->date);
            phaptmp->date = NULL;
        }

    }

    //free country

    int *value;
    hashmap_foreach_data(value, &((*pHaps)->country)) {
        free(value);
    }
    hashmap_cleanup(&((*pHaps)->country));

    free(*pHaps);
    *pHaps = NULL;
}