//
// Created by lun on 2021/3/8.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include "haplotype.h"


#include "virus.h"
#include "haplotype.h"

#define DEBUG

void hap_init(struct Hap * pHaps, int numOfHaps);

int intLength(int n){
    int i = 0;
    while(n!=0)    {
        n = n / 10;//每次除以10
        i++;//统计循环次数
    }
    return i;
}

void muts2str(struct Mutation *muts, int numOfMuts, char **str){
    if(numOfMuts == 0){
        *str = (char *)malloc(sizeof(char) * (3+1));
        strcpy(*str, "N/A");
        return;
    }
    //calc length
    int l = 0;
    for (int i = 0; i < numOfMuts; ++i) {
        l += intLength(muts[i].pos)+1;
        l += (int )strlen(muts[i].ref)+1;
        l += (int )strlen(muts[i].alt)+1;
    }
    (*str) = (char *)malloc(sizeof(char) * (l+1));
    //char *str_tmp = (char *)malloc(sizeof(char) * (l+1));
    for (int i = 0; i < numOfMuts; ++i) {
        sprintf(*str, "%d,%s,%s,", muts[0].pos, muts[0].ref, muts[0].alt);
    }
    for (int i = 1; i < numOfMuts; ++i) {
        int l_curent = 0;
        sprintf( *str , "%s%d,%s,%s,", *str, muts[i].pos, muts[i].ref, muts[i].alt);
    }
}

void findHap_hashmap(struct Virus * viruses, int numOfViruses, struct Hap ** haps, int * numOfHaps){
    int * virus2hap = (int *)malloc(sizeof(int) * numOfViruses);//hap ID of sequence ID is seq2hap(seqID)
    int * hap2numOfVirus = (int *)malloc(sizeof(int) * numOfViruses);
    int * hap2virus = (int *)malloc(sizeof(int) * numOfViruses);
    //init
    for (int i = 0; i < numOfViruses; ++i) {
        virus2hap[i] = -1;
    }
    for (int i = 0; i < numOfViruses; ++i) {
        hap2numOfVirus[i] = 0;
    }
    for (int i = 0; i < numOfViruses; ++i) {
        hap2virus[i] = 0;
    }
    *numOfHaps = 0;

    //find hap

    HASHMAP(char, struct Virus) MutVirus;
    hashmap_init(&MutVirus, hashmap_hash_string, strcmp);
    for (int i = 0; i < numOfViruses; ++i) { //ith virus
        char *str = NULL;
        muts2str(viruses[i].mutations, viruses[i].numOfMut, &str);
        int r = hashmap_put(&MutVirus, str, &viruses[i]);
        if(r < 0){
            struct Virus *pVirus = hashmap_get(&MutVirus, str);
            free(str);
            str = NULL;
            virus2hap[i] = virus2hap[pVirus - viruses];
        }
        else{
            virus2hap[i] = *numOfHaps;
            (*numOfHaps) ++;
            hap2virus[virus2hap[i]] = i;
        }
        hap2numOfVirus[virus2hap[i]] ++;
    }

    *haps = (struct Hap * )malloc(sizeof(struct Hap) * (*numOfHaps));
    hap_init(*haps, *numOfHaps);

    for (int i = 0; i < *numOfHaps; ++i) {
        (*haps)[i].numOfVirus = hap2numOfVirus[i];
        (*haps)[i].allViruses = (int *)malloc(sizeof(int) * ((*haps)[i].numOfVirus));
        (*haps)[i].numOfMut = viruses[hap2virus[i]].numOfMut;
        (*haps)[i].mutations = (struct Mutation*)malloc(sizeof(struct Mutation)*(*haps)[i].numOfMut);
        for (int j = 0; j < (*haps)[i].numOfMut; ++j) {
            (*haps)[i].mutations[j] = viruses[hap2virus[i]].mutations[j];
        }
        (*haps)[i].index = i;
        (*haps)[i].date = (struct Date *)malloc(sizeof(struct Date));
        *((*haps)[i].date) = *(viruses[hap2virus[i]].date);
    }

    int *hap2currentNumOfVirus = (int * )malloc(sizeof(int) * (*numOfHaps));
    for (int i = 0; i < *numOfHaps; ++i) {
        hap2currentNumOfVirus[i] = 0;
    }
    for (int i = 0; i < numOfViruses; ++i) {
        int hapIndex = virus2hap[i];
        if(hap2currentNumOfVirus[hapIndex] < (*haps)[hapIndex].numOfVirus){
            ((*haps)[hapIndex].allViruses)[hap2currentNumOfVirus[hapIndex]] = i;
            hap2currentNumOfVirus[hapIndex]++;
        }
    }

    //free
    const char * key;
    struct Virus * value;
    hashmap_foreach_key(key,&MutVirus){
        free((char*)key);
        key = NULL;
    }
    hashmap_cleanup(&MutVirus);

    free(hap2currentNumOfVirus);
    hap2currentNumOfVirus = NULL;
    free(virus2hap);
    virus2hap = NULL;
    free(hap2numOfVirus);
    hap2numOfVirus = NULL;
    free(hap2virus);
    hap2virus = NULL;
}