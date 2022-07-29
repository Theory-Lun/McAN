//
// Created by lun on 3/3/21.
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
//#include <unistd.h>

#include "virus.h"
#include "haplotype.h"
#include "edges.h"

//#define DEBUG
#define DEBUG_TIME

void help();
void readNumOfVirusesFromMut(const char * fin, int * numOfViruses);
void readMutaions(const char * fin, bool * mutKeep, const int maxLengthOfSeq, int numOfViruses, struct Virus * viruses);
void showVirus(struct Virus * viruses, int numOfViruses);
void virus_init(struct Virus * viruses, int numOfviruses);
void hap_show(struct Hap * pHaps, int numOfHaps);
void edges_show(struct Edge * edges, int numOfEdges);
char * creat_network(struct Virus * viruses, int numOfViruses, struct Hap * haps, int numOfHaps, struct Edge * edges, int numOfEdges, int * network, int numOfEdgesInNet, int heightOfTree);
void countryDistribution(struct Hap * haps, int numOfHaps, struct Virus * viruses, int numOfViruses);
void output(const char * str_folder, char * json, long memoryUsed, double time);
void free_viruses(struct Virus ** viruses, int numOfViruses);
void readMeta_hash(const char * fin, int numOfViruses, struct Virus * viruses);
int cmp_time_virus(const void* _a, const void* _b);
void treeHeight(struct Hap * haps, int numOfHaps, int * height);
void show_mutation(struct Mutation *pMut, int numOfIndel);

inline bool isFileExist(const char* fileName) {
    FILE*  fp = fopen(fileName, "r");
    bool ret = (fp != NULL);
    if (ret) fclose(fp);
    return ret;
}

bool checkUsageAll(const char * str_vcf, const char * str_mut, const char * str_sta_anno, const char * str_meta, const char * str_minFreq, const char * str_maxvirus, const char * str_folder,  double * minFreq, int * maxvirus ){
    //some arguments do not exists
    if( NULL == str_vcf && NULL == str_mut || NULL == str_meta || NULL == str_folder ){
        printf("error: no %s, %s or %s\n", "vcf|mutation", "meta", "folder");
        return false;
    }
    //check if format of folder is good
    if('/' == str_folder[strlen(str_folder)-1]){
        printf("error: folder end with '/'\n");
        return false;
    }
    //files can not read
    if(str_vcf && !(str_vcf[0] == '-' && str_vcf[1] == 0) && !isFileExist(str_vcf)){
        printf("error: can not read vcf file: \n%s\n", str_vcf);
        return false;
    }
    if(str_mut && !isFileExist(str_mut)){
        printf("error: can not read mutation file: \n%s\n", str_mut);
        return false;
    }
    if(!isFileExist(str_meta)){
        printf("error: can not read meta file: \n%s\n", str_meta);
        return false;
    }
    //output folder do not exists
    //if(-1 == access(str_folder, F_OK)){
    //    printf("error: output folder do not exists: \n%s\n", str_folder);
    //    return false;
    //}
    //check if numbers are good
    if(NULL != str_minFreq){
        char *endptr = NULL; // for strtol
        *minFreq = strtod(str_minFreq, &endptr);
        if('\0' != endptr[0] || '\0' == str_minFreq[0]){
            printf("error: minfreq is not a real number:\n%s\n", str_minFreq);
            return false;
        }
        endptr = NULL;

        if(*minFreq < 0 || *minFreq >= 1.0){
            printf("error: minfreq = %f\n", *minFreq);
            printf("minfreq must belong to [0,1)\n");
            return false;
        }
    }

    if(NULL != str_maxvirus){
        char *endptr = NULL; // for strtol
        *maxvirus = (int)strtol(str_maxvirus, &endptr, 10);
        if(endptr[0] != '\0' || str_maxvirus[0] == '\0'){
            printf("error: maxnsample is not a integer:\n%s\n", str_maxvirus);
            return false;
        }
        endptr = NULL;
        if(*maxvirus <= 0){
            printf("error: maxnsample = %d\n", *maxvirus);
            printf("maxnsample must > 0\n");
            return false;
        }
    }
#ifdef DEBUG
    printf("\n###CHECK AUG START###\n");
    printf("minFreq = %f\n", *minFreq);
    printf("maxvirus = %d\n", *maxvirus);
    printf("###CHECK AUG END###\n");
#endif
    return true;
}
