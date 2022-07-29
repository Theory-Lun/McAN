//
// Created by lun on 3/4/21.
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "virus.h"
#include <stdbool.h>
#include <hashmap.h>

bool chomp( char *str );
char* strip(char* str);
void readDate2Virus(char * str_date, struct Virus * pVirus);
void show_date(struct Date * date);
void showVirus(struct Virus * viruses, int numOfViruses);
void show_location(struct Location * pLocation);
void location_init(struct Location * pLocation);
void date_init(struct Date * date);
bool emptyStr(char * str);

//#define DEBUG
void readMeta_hash(const char * fin, int numOfViruses, struct Virus * viruses){
    //construct hasmap acc2virus
    HASHMAP(char, struct Virus) map;
    const char *key;
    struct Virus *value;
    int r_hashmap;
    hashmap_init(&map, hashmap_hash_string, strcmp);
    for (int i = 0; i < numOfViruses; ++i) {
        if(NULL != &(viruses[i])){
            r_hashmap = hashmap_put(&map, (viruses[i]).acc, &(viruses[i]));
        }
        if(r_hashmap < 0){
            fprintf(stderr, "same acc in different viruses\n");
            exit(1);
        }
    }

    if(NULL == fin){
        return;
    }
    FILE * fpin = fopen(fin, "r");
    if (NULL == fpin){
        printf("could not open file:\n%s\n", fin);
        exit(EXIT_FAILURE);
    }

    char *line = NULL;
    char* virusName = NULL;
    size_t len = 0;
    ssize_t read;

    int indexOfline = 0;//from 0
    int numOfMatchedViruses = 0;
    while ((read = getline(&line, &len, fpin)) != -1) {//read line by line

        chomp(line);
        //char *ptr1 = NULL;//for strtok_r

        int columnIndex = 0; //column index starting from 0
        struct Virus * pVirus = NULL;
        char * r = line;

        for (char * token1 =strsep(&r,"\t"); token1!=NULL ; token1=strsep(&r,"\t")) {
            //if all viruses are matched, break
            if(numOfMatchedViruses > numOfViruses){
                //break;
            }

            token1 = strip(token1);
            if (token1[0] == '*' && token1[1] == 0)
                token1[0] = 0;

            if (0 == columnIndex) {//1st column, virus strain name
                virusName = (char *)realloc(virusName, sizeof(char) * (1 + strlen(token1)));
                if(NULL == virusName)
                    exit(EXIT_FAILURE);
                strcpy(virusName, token1);
            }
            else if(1 == columnIndex){ //2ed column, acc
                //seach each virus, find associated virus
                //pVirus = findVirusOfAcc(token1, numOfViruses, viruses);
                pVirus = hashmap_get(&map, token1);
                if(NULL == pVirus){
                    break;//next metadata (row)
                }
                else{
#ifdef DEBUG
                    //printf("\n%s\n", line);
                    //printf("%d\n", indexOfline);
#endif
                    pVirus->location = (struct Location *)malloc(sizeof(struct Location));
                    location_init(pVirus->location);
                    pVirus->date = (struct Date *)malloc(sizeof(struct Date));
                    date_init(pVirus->date);
                    numOfMatchedViruses ++;

                    if (!pVirus->name || !pVirus->name[0]) {
                        char* tmp = pVirus->name;
                        pVirus->name = virusName;
                        virusName = tmp;
                    }
                }
            }
            //else if(7 == columnIndex){ //8th column, date
            else if (2 == columnIndex) { //date

#ifdef DEBUG
                printf("%s\n", token1);
#endif

                if(emptyStr(token1)){
                    columnIndex ++;
                    continue;
                }
                readDate2Virus(token1, pVirus);
            }
            //else if(12 == columnIndex){//country
            else if (3 == columnIndex) {//country
                if(emptyStr(token1)){
                    columnIndex ++;
                    continue;
                }
                if(NULL == pVirus->location->country){
                    pVirus->location->country = (char *)malloc(sizeof(char) * (1+strlen(token1)));
                }
                strcpy(pVirus->location->country, token1);
            }
            //else if(13 == columnIndex){//state
            else if (4 == columnIndex) {//state
                if(emptyStr(token1)){
                    columnIndex ++;
                    continue;
                }
                if(NULL == pVirus->location->state){
                    pVirus->location->state = (char *)malloc(sizeof(char) * (1+strlen(token1)));
                }
                strcpy(pVirus->location->state, token1);
            }
            //else if(14 == columnIndex){//city
            else if (5 == columnIndex) {//city
                if(emptyStr(token1)){
                    columnIndex ++;
                    continue;
                }
                if(NULL == pVirus->location->city){
                    pVirus->location->city = (char *)malloc(sizeof(char) * (1+strlen(token1)));
                }
                strcpy(pVirus->location->city, token1);
            }
            columnIndex ++;
        }
        indexOfline ++;
    }

#ifdef DEBUG
    /*
    struct Virus * asdf = hashmap_get(&map, "EPI_ISL_857152");
    if(asdf == NULL){
        printf("asdf == NULL\n");
    }
    else{
        printf("asdf != NULL\n");
        if(NULL == asdf->location){
            printf("NULL == asdf->location\n");
        }
        else{
            printf("%s\n", asdf->location->country);
        }

    }
    exit(1);
     */
    //printf("numOfMatchedViruses = %d\n", numOfMatchedViruses);
#endif
    //close
    fclose(fpin);
    fpin = NULL;
    free(line);
    line= NULL;
    free(virusName);
    virusName = NULL;
    hashmap_cleanup(&map);
}
