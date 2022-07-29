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
void extractDateFromMetaLine(char *line, struct Date *pDate);

void readMeta_acc2line_hash(char *str_input, struct Hashmap_acc2line *hash){
    hashmap_init(&(hash->acc2line), hashmap_hash_string, strcmp);
    FILE *fpin = fopen(str_input, "r");
    if(NULL == fpin){
        exit(1);
    }
    char *line = NULL;
    size_t len = 0;
    ssize_t read;

    int indexOfline = 0;//from 0
    while ((read = getline(&line, &len, fpin)) != -1) {//read line by line
        chomp(line);
        int columnIndex = 0; //column index starting from 0
        char *r = line;
        char *lineCopy = (char *)malloc(sizeof(char) * (strlen(line) + 1));
        strncpy(lineCopy, line, strlen(line) + 1);
        for (char *token1 = strsep(&r, "\t"); token1 != NULL; token1 = strsep(&r, "\t")) {
            token1 = strip(token1);
            if (1 == columnIndex) { //2ed column, acc
                char *token1Copy = (char *)malloc(sizeof(char) * (strlen(token1) + 1));
                strncpy(token1Copy, token1, strlen(token1) + 1);
                hashmap_put(&(hash->acc2line), token1Copy, lineCopy);
            }
            columnIndex ++;
        }
        indexOfline ++;
    }
    //close
    fclose(fpin);
    fpin = NULL;
}


void readMeta_acc2Date_hash(char *str_input, struct Hashmap_acc2Date *hash){
    hashmap_init(&(hash->acc2date), hashmap_hash_string, strcmp);
    struct Hashmap_acc2line acc2Meta;
    readMeta_acc2line_hash(str_input, &acc2Meta);
    struct Date StructDate;
    const char *acc;
    char *line;
    hashmap_foreach(acc, line, &(acc2Meta.acc2line)) {
        extractDateFromMetaLine(line, &StructDate);
        char *accCp = (char *)malloc(sizeof(char)*(strlen(acc)+1));
        strcpy(accCp, acc);
        struct Date *DateCp = (struct Date *)malloc(sizeof(struct Date));
        *DateCp = StructDate;
        hashmap_put(&(hash->acc2date), accCp, DateCp);
    }
    //free acc2Meta
    hashmap_foreach(acc, line, &(acc2Meta.acc2line)) {
        free((char*)acc);
        acc = NULL;
        free(line);
        line = NULL;
    }
    hashmap_cleanup(&(acc2Meta.acc2line));
}

