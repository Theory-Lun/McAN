//
// Created by lun on 3/5/21.
//
#include <stdio.h>
#include "string.h"

#include "virus.h"
void show_location(struct Location * pLocation);
void show_date(struct Date * date);
void show_mutation(struct Mutation *pMut, int numOfIndel);

void showMut(struct Mutation * mut, int numOfMut){
    if(NULL == mut || 0 == numOfMut){
        return;
    }
    char str[20];

    for (int i = 0; i < numOfMut; ++i) {
        if(MUTATION_INDEL == mut[i].mutationType){
            strcpy(str, "Indel");
        }
        else if(MUTATION_DELETION == mut[i].mutationType){
            strcpy(str, "Deletion");
        }
        else if(MUTATION_INSERTION == mut[i].mutationType){
            strcpy(str, "Insertion");
        }
        else if(MUTATION_SNP == mut[i].mutationType){
            strcpy(str, "SNP");
        }
        printf("%d(%s:%s->%s);", mut[i].pos, str, mut[i].ref, mut[i].alt);
    }
    printf("\n");
}

void showVirus(struct Virus * viruses, int numOfViruses){
    for (int i = 0; i < numOfViruses; ++i) {
        printf("\n");
        printf("%d\n", i);
        if(NULL != viruses[i].name){
            printf("name:   %s\n", viruses[i].name);
        }
        if(NULL != viruses[i].acc){
            printf("acc:    %s\n", viruses[i].acc);
        }
        if(viruses[i].numOfMut >= 0){
            printf("mutation number:   %d\n", viruses[i].numOfMut);
        }
        if(NULL != viruses[i].mutations){
            show_mutation(viruses[i].mutations, viruses[i].numOfMut);
        }

        if(NULL != viruses[i].location){
            show_location(viruses[i].location);
        }
        if(NULL != viruses[i].date){
            show_date(viruses[i].date);
        }
        printf("numOfIndel = %d\n",viruses[i].numOfMut);
    }
}