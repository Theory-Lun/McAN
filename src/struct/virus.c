//
// Created by lun on 3/5/21.
//

#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#include "virus.h"

int find_string(char * str, char * substr);

//#define DEBUG
void virus_init(struct Virus * viruses, int numOfViruses){
    for (int i = 0; i < numOfViruses; ++i) {
        viruses[i].mutations = NULL;
        viruses[i].date = NULL;
        viruses[i].numOfMut = -1;
        viruses[i].name = NULL;
        viruses[i].acc = NULL;
        viruses[i].location = NULL;
        viruses[i].indexOfHap = -1;
    }
}

void show_date(struct Date * date){
    printf("%04d-%02d-%02d\n", date->year, date->month, date->day);
}

void date_init(struct Date * date){
    date->year = -1;
    date->month = -1;
    date->day = -1;
}

void readDate2Date(char * str_date, struct Date * pDate){
    if(NULL == pDate){
        pDate = (struct Date *)malloc(sizeof(struct Date));
    }
    date_init(pDate);

    int numOfSeparator = find_string(str_date, "-");
#ifdef DEBUG
    printf("%s\n%d\n", str_date, numOfSeparator);
#endif

    if(2 == numOfSeparator){
        sscanf(str_date, "%hd-%hd-%hd", &(pDate->year), &(pDate->month), &(pDate->day));
        return;
    }
    else if(1 == numOfSeparator){
        sscanf(str_date, "%hd-%hd", &(pDate->year), &(pDate->month));
        return;
    }
    else if(0 == numOfSeparator){
        sscanf(str_date, "%hd", &(pDate->year));
        return;
    }
    printf("date format error\n");
    exit(1);
}

void readDate2Virus(char * str_date, struct Virus * pVirus){
    if(NULL == pVirus){
        printf("struct Virus * pVirus = NULL\n");
        exit(1);
    }
    if(NULL == pVirus->date){
        pVirus->date = (struct Date *)malloc(sizeof(struct Date));
    }
    date_init(pVirus->date);
#ifdef DEBUG
    //printf("%lu\t%c\t%c\n", strlen(str_date), str_date[4], str_date[7]);
#endif
    readDate2Date(str_date, pVirus->date);
}

void show_location(struct Location * pLocation){
    if(NULL == pLocation){
        printf("location = NULL\n");
    }
    else{
        if(NULL != pLocation->country){
            printf("country = %s\n", pLocation->country);
        }
        if(NULL != pLocation->state){
            printf("state = %s\n", pLocation->state);
        }
        if(NULL != pLocation->city){
            printf("city = %s\n", pLocation->city);
        }
    }
}

void location_init(struct Location * pLocation){
    if(NULL == pLocation){
        printf("location = NULL\n");
    }
    pLocation->country = NULL;
    pLocation->state = NULL;
    pLocation->city = NULL;
}

void date2str(struct Date * date, char ** str){
    if(NULL == date){
        return;
    }
    if(NULL == str){
        *str = (char *)malloc(sizeof(char) * 11);//2020-01-02
    }
    if(date->year >= 0){
        int j = sprintf(*str, "%04d", date->year);
        if(date->month >= 1 && date->month <= 12){
            j += sprintf((*str)+j, "-%02d", date->month);
            if(date->day >= 1 && date->day <= 31){
                sprintf((*str)+j, "-%02d", date->day);
            }
        }
    }
}

void free_viruses(struct Virus ** viruses, int numOfViruses){
    for (int i = 0; i < numOfViruses; ++i) {
        struct Virus * virus = &((*viruses)[i]);
        if(NULL != virus->acc){
            free(virus->acc);
            virus->acc = NULL;
        }
        if(NULL != virus->name){
            free(virus->name);
            virus->name = NULL;
        }
        if(NULL != virus->location){
            //free location
            if(NULL != virus->location->country){
                free(virus->location->country);
                virus->location->country = NULL;
            }
            if(NULL != virus->location->state){
                free(virus->location->state);
                virus->location->state = NULL;
            }
            if(NULL != virus->location->city){
                free(virus->location->city);
                virus->location->city = NULL;
            }
            free(virus->location);
            virus->location = NULL;
        }

        if(NULL != virus->date){
            //free date
            free(virus->date);
            virus->date = NULL;
        }

        if(NULL != virus->mutations){
            //free mutations
            if(NULL != virus->mutations->ref){
                free(virus->mutations->ref);
                virus->mutations->ref = NULL;
            }
            if(NULL != virus->mutations->alt){
                free(virus->mutations->alt);
                virus->mutations->alt = NULL;
            }
            free(virus->mutations);
            virus->mutations = NULL;
        }
    }
    free(* viruses);
    *viruses = NULL;
}
