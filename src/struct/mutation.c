//
// Created by lun on 3/5/21.
//
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include "mutation.h"

void show_mutation(struct Mutation *pMut, int numOfIndel);

//#define DEBUG
bool isSameMut(struct Mutation * mut1, struct Mutation * mut2){
    if(NULL == mut1 || NULL == mut2){
        fprintf(stderr, "mut = NULL\n");
        exit(1);
    }
    if(NULL == mut1->alt || NULL == mut1->ref || NULL == mut2->alt || NULL == mut2->ref){
        fprintf(stderr, "ref or alt error\n");
        exit(1);
    }
#ifdef DEBUG
    show_mutation(mut1, 1);
    show_mutation(mut2, 1);
#endif
    if( mut1->pos == mut2->pos && 0 == strcmp(mut1->alt, mut2->alt) && 0 == strcmp(mut1->ref, mut2->ref)){
        return true;
    }
    else{
        return false;
    }
}

void show_mutation(struct Mutation *pMut, int numOfIndel){
    char tmp[100];
    printf("Mutation(s):\n");
    for (int i = 0; i < numOfIndel; ++i) {
        if(MUTATION_INDEL == pMut[i].mutationType ){
            strcpy(tmp, "Indel");
        }
        else if(MUTATION_DELETION == pMut[i].mutationType ){
            strcpy(tmp, "Deletion");
        }
        else if(MUTATION_INSERTION == pMut[i].mutationType ){
            strcpy(tmp, "Insertion");
        }
        else if(MUTATION_SNP == pMut[i].mutationType ){
            strcpy(tmp, "SNP");
        }
        printf("%d(%s:%s->%s);", pMut[i].pos, tmp, pMut[i].ref, pMut[i].alt);
    }
    printf("\n");
}