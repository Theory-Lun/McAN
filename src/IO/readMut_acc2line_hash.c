//
// Created by lun on 5/18/21.
//
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <hashmap.h>
#include <string.h>
#include "Hashmap_acc2line.h"
bool chomp( char *str ) ;
char* strip(char* str);
typedef int(*Callback)(const char*, size_t, size_t, size_t, void*);
int vcf2mut(const char* vcfName, const char* mutName, Callback lineHandler, void* arg, const char* missingValue);

void readMutFromLine_acc2line_hash(char* line, struct Hashmap_acc2line *hash) {
        int columnIndex = 0; //column index starting from 0
        char *r = line;
        char *lineCopy = (char *)malloc(sizeof(char) * (strlen(line) + 1));
        strncpy(lineCopy, line, strlen(line) + 1);
        for (char *token1 = strsep(&r, "\t"); token1 != NULL; token1 = strsep(&r, "\t")) {
            token1 = strip(token1);
            if(1 == columnIndex){//2ed column, accession ID
                char *token1Copy = (char *)malloc(sizeof(char) * (strlen(token1) + 1));
                strncpy(token1Copy, token1, strlen(token1) + 1);
                hashmap_put(&(hash->acc2line), token1Copy, lineCopy);
            }
            columnIndex ++;
        }
}

void readMut_acc2line_hash(char *str_input, struct Hashmap_acc2line *hash){
    hashmap_init(&(hash->acc2line), hashmap_hash_string, strcmp);
    FILE *fpin = fopen(str_input, "r");
    if (NULL == fpin) {
        exit(1);
    }
    char *line = NULL;
    size_t len = 0;
    ssize_t read;

    int indexOfline = 0;//from 0
    while ((read = getline(&line, &len, fpin)) != -1) {//read line by line
        chomp(line);
        readMutFromLine_acc2line_hash(line, hash);
        indexOfline ++;
    }
    //close
    fclose(fpin);
    fpin = NULL;

}

int lineHandler_acc2line_hash(const char* line, size_t len, size_t nonuse1, size_t nonuse2, void* arg) {
    if (!line || !arg)
        return -1;

    struct Hashmap_acc2line* hash = (struct Hashmap_acc2line*)arg;
    char* line_cp = malloc((len + 1) * sizeof(char));
    memcpy(line_cp, line, len + 1);
    readMutFromLine_acc2line_hash(line_cp, hash);
    free(line_cp);

    return 0;
}

void readVCF_acc2line_hash(const char* str_vcf, const char* str_mut, struct Hashmap_acc2line *hash) {
    hashmap_init(&(hash->acc2line), hashmap_hash_string, strcmp);
    if (vcf2mut(str_vcf, str_mut, lineHandler_acc2line_hash, (void*)hash, "*") < 0)
        exit(EXIT_FAILURE);
}

