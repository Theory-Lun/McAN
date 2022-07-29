//
// Created by lun on 3/12/21.
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//#define DEBUG

void output_json(const char * str_folder, char * filename, char * json){
    char * str_json = NULL;
    size_t l = strlen(str_folder) + strlen(filename) + 2;//2 for '\0' and '/'
    str_json = (char *)malloc(sizeof(char)*l);
    if(NULL == str_json){
        exit(1);
    }
    strcpy(str_json, str_folder);
    strcat(str_json, "/");
    strcat(str_json, filename);
#ifdef DEBUG
    printf("%s\n", str_json);
#endif

    FILE * fpout = fopen(str_json, "w");
    if(NULL == fpout){
        printf("open file error!\n");
        exit(EXIT_FAILURE);
    }

    fprintf(fpout, "%s", json);

    fclose(fpout);
    fpout = NULL;
    free(str_json);
    str_json = NULL;
}

void output_time(const char * str_folder, char * filename, double time){
    char * str_time = NULL;
    size_t l = strlen(str_folder) + strlen(filename) + 2;//2 for '\0' and '/'
    str_time = (char *)malloc(sizeof(char)*l);
    if(NULL == str_time){
        exit(1);
    }
    strcpy(str_time, str_folder);
    strcat(str_time, "/");
    strcat(str_time, filename);
#ifdef DEBUG
    printf("%s\n", str_json);
#endif
    FILE * fpout = fopen(str_time, "w");
    if(NULL == fpout){
        printf("open file error!\n");
        exit(EXIT_FAILURE);
    }

    fprintf(fpout, "%f s", time);

    fclose(fpout);
    fpout = NULL;
    free(str_time);
    str_time = NULL;
}

void output(const char * str_folder, char * json, long memoryUsed, double time){
    output_json(str_folder, "haplotype_loci.json", json);
    //output_memory(str_folder, "memory", memoryUsed);
    output_time(str_folder, "time", time);
}