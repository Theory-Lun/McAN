//
// Created by lun on 3/4/21.
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#include "virus.h"
#include "indel.h"

bool chomp( char *str );
char* strip(char* str);
int find_string(char * str, char * substr);
void showMut(struct Mutation * mut, int numOfMut);
void virus_init(struct Virus * viruses, int numOfviruses);
typedef int(*Callback)(const char*, size_t, size_t, size_t, void*);
int vcf2mut(const char* vcfName, const char* mutName, Callback lineHandler, void* arg, const char* missingValue);

//#define DEBUG
void readNumOfVirusesFromMut(const char * fin, int * numOfViruses){
    //char * fin: mutation file
    //int * numOfViruses: number of viruses in mutation file

    FILE * fpin = fopen(fin, "r");
    if (NULL == fpin){
        printf("could not open file:\n%s\n", fin);
        exit(EXIT_FAILURE);
    }

    char *line = NULL;
    size_t len = 0;
    ssize_t read;
    int indexOfline = 0;//from 0

    while ((read = getline(&line, &len, fpin)) != -1) {//read line by line
        chomp(line);
        indexOfline ++;
    }
    *numOfViruses = indexOfline;

    //close
    fclose(fpin);
    fpin = NULL;

    free(line);
    line = NULL;

#ifdef DEBUG
    printf("\n###READ NUMBER OF VIRUSES START###\n");
    printf("numOfViruses = %d\n", *numOfViruses);
    printf("###READ NUMBER OF VIRUSES END###\n");
#endif
}

void readMutaionsFromLine(char* line, int indexOfline, bool * mutKeep, struct Virus * viruses, int* lenthOfTmpMut, struct Mutation * tmpMut, char *str_ref, char *str_alt) {
        char *ptr1 = NULL;//for strtok_r
        char *ptr2 = NULL;//for strtok_r
        int columnIndex = 0; //column index starting from 0
        for( char * token1 =strtok_r(line,"\t",&ptr1); token1!=NULL ; token1=strtok_r(NULL,"\t",&ptr1) ) {
#ifdef DEBUG
            printf("%d, %s\n", columnIndex, token1);
#endif
            token1 = strip(token1);
            if (token1[0] == '*' && token1[1] == 0)
                token1[0] = 0;

            if(0 == columnIndex){//1st column, virus strain name
                viruses[indexOfline].name = (char *)malloc(sizeof(char) * (1+ strlen(token1)));
                strcpy(viruses[indexOfline].name, token1);
            }
            else if(1 == columnIndex){//2ed column, accession ID
                viruses[indexOfline].acc = (char *)malloc(sizeof(char) * (1+ strlen(token1)));
                strcpy(viruses[indexOfline].acc, token1);
            }
            else if(2 == columnIndex){//3ed column, mutations
                int numOfSnps = find_string(token1, "(SNP:");
                if(numOfSnps > 30000){
                    printf("number of SNPs > 30000\n");
                    exit(1);
                }
                if(numOfSnps > *lenthOfTmpMut){
                    free(tmpMut);
                    *lenthOfTmpMut = numOfSnps;
                    tmpMut = (struct Mutation *)malloc(sizeof(struct Mutation)*(*lenthOfTmpMut));
                }
#ifdef DEBUG
                //printf("num of all SNPs = %d\n", numOfSnps);
#endif
                //read each mutation
                int numOfFeasibleMut = 0;
                for( char * token2 =strtok_r(token1,";",&ptr2) ; token2!=NULL ; token2=strtok_r(NULL,";",&ptr2) ){
#ifdef DEBUG
                    //printf("%s\n", token2);
#endif
                    int pos;
                    char str_type[10];//SNP, Deletion, Insertion, Indel
                    sscanf(token2, "%d(%9[^:]%*s", &pos, str_type);//example: 25452(SNP:C->T)
#ifdef DEBUG
                    //printf("pos = %d\n", pos);
                    //printf("type = %s\n", str_type);
#endif
                    //ignore pos not in mutKeep
                    if(false == mutKeep[pos-1]){
#ifdef DEBUG
                        //printf("ignore. pose\n");
#endif
                        continue;
                    }

                    sscanf(token2, "%*[0-9](%*[^:]:%[^-]->%[^)])", str_ref, str_alt);
#ifdef DEBUG
                    //printf("\n%s\n%s\n%s\n", token2, str_ref, str_alt);
#endif

#ifdef TEST_SPLIT_DELETION
                    if(0 == strcmp(str_type, "Deletion")){
                        if(str_ref[0] == '-'){
                            tmpMut[numOfFeasibleMut].ref = (char *)malloc(sizeof(char) * 2);
                            tmpMut[numOfFeasibleMut].ref[0] = str_ref[0];
                            tmpMut[numOfFeasibleMut].ref[1] = '\0';

                            tmpMut[numOfFeasibleMut].alt = (char *)malloc(sizeof(char) * 2);
                            tmpMut[numOfFeasibleMut].alt[0] = '-';
                            tmpMut[numOfFeasibleMut].alt[1] = '\0';

                            tmpMut[numOfFeasibleMut].pos = (short )pos;
                            tmpMut[numOfFeasibleMut].mutationType = INDEL_DELETION;
                            numOfFeasibleMut ++;
                        }
                        for (int i = 1; i < strlen(str_ref); ++i) {
                            tmpMut[numOfFeasibleMut].ref = (char *)malloc(sizeof(char) * 2);
                            tmpMut[numOfFeasibleMut].ref[0] = str_ref[i];
                            tmpMut[numOfFeasibleMut].ref[1] = '\0';

                            tmpMut[numOfFeasibleMut].alt = (char *)malloc(sizeof(char) * 2);
                            tmpMut[numOfFeasibleMut].alt[0] = '-';
                            tmpMut[numOfFeasibleMut].alt[1] = '\0';

                            tmpMut[numOfFeasibleMut].pos = (short )(pos+i);
                            tmpMut[numOfFeasibleMut].mutationType = INDEL_DELETION;
                            numOfFeasibleMut ++;
                        }

                    }
#else
                    if(0 == strcmp(str_type, "Deletion")){
                        tmpMut[numOfFeasibleMut].ref = (char *)malloc(sizeof(char)*(1+strlen(str_ref)));
                        strncpy(tmpMut[numOfFeasibleMut].ref, str_ref, 1+strlen(str_ref));
                        tmpMut[numOfFeasibleMut].alt = (char *)malloc(sizeof(char)*(1+strlen(str_alt)));
                        strncpy(tmpMut[numOfFeasibleMut].alt, str_alt, 1+strlen(str_alt));
                        tmpMut[numOfFeasibleMut].pos = (short )pos;
                        tmpMut[numOfFeasibleMut].mutationType = INDEL_DELETION;
                        numOfFeasibleMut ++;
                    }
#endif
                    else if( 0 == strcmp(str_type, "Insertion") ){
                        tmpMut[numOfFeasibleMut].ref = (char *)malloc(sizeof(char)*(1+strlen(str_ref)));
                        strncpy(tmpMut[numOfFeasibleMut].ref, str_ref, 1+strlen(str_ref));
                        tmpMut[numOfFeasibleMut].alt = (char *)malloc(sizeof(char)*(1+strlen(str_alt)));
                        strncpy(tmpMut[numOfFeasibleMut].alt, str_alt, 1+strlen(str_alt));
                        tmpMut[numOfFeasibleMut].pos = (short )pos;
                        tmpMut[numOfFeasibleMut].mutationType = INDEL_INSERTION;
                        numOfFeasibleMut ++;
                    }
                    else if( 0 == strcmp(str_type, "Indel") ){
                        tmpMut[numOfFeasibleMut].ref = (char *)malloc(sizeof(char)*(1+strlen(str_ref)));
                        strncpy(tmpMut[numOfFeasibleMut].ref, str_ref, 1+strlen(str_ref));
                        tmpMut[numOfFeasibleMut].alt = (char *)malloc(sizeof(char)*(1+strlen(str_alt)));
                        strncpy(tmpMut[numOfFeasibleMut].alt, str_alt, 1+strlen(str_alt));
                        tmpMut[numOfFeasibleMut].pos = (short )pos;
                        tmpMut[numOfFeasibleMut].mutationType = INDEL_INDEL;
                        numOfFeasibleMut ++;
                    }
                    else if( 0 == strcmp(str_type, "SNP") ){
                        tmpMut[numOfFeasibleMut].ref = (char *)malloc(sizeof(char)*(1+strlen(str_ref)));
                        strncpy(tmpMut[numOfFeasibleMut].ref, str_ref, 1+strlen(str_ref));
                        tmpMut[numOfFeasibleMut].alt = (char *)malloc(sizeof(char)*(1+strlen(str_alt)));
                        strncpy(tmpMut[numOfFeasibleMut].alt, str_alt, 1+strlen(str_alt));
                        tmpMut[numOfFeasibleMut].pos = (short )pos;
                        tmpMut[numOfFeasibleMut].mutationType = INDEL_SNP;
                        numOfFeasibleMut ++;
                    }
                }

                //creat mutations
                viruses[indexOfline].numOfMut = numOfFeasibleMut;
                if( 0 == numOfFeasibleMut ){
                    viruses[indexOfline].mutations = NULL;
                }
                else{
                    viruses[indexOfline].mutations = (struct Mutation *)malloc(sizeof(struct Mutation) * numOfFeasibleMut);
                    if(NULL == viruses[indexOfline].mutations){
                        fprintf(stderr,"malloc failed\n");
                        exit(1);
                    }
                }

                for (int i = 0; i < numOfFeasibleMut; ++i) { //i: index of mutations
                    (viruses[indexOfline].mutations)[i].pos = tmpMut[i].pos;

                    (viruses[indexOfline].mutations)[i].ref = (char *)malloc(sizeof(char) *(1+strlen(tmpMut[i].ref)));
                    strncpy((viruses[indexOfline].mutations)[i].ref, tmpMut[i].ref, 1+strlen(tmpMut[i].ref));
                    free(tmpMut[i].ref);
                    tmpMut[i].ref = NULL;

                    (viruses[indexOfline].mutations)[i].alt = (char *)malloc(sizeof(char) *(1+strlen(tmpMut[i].alt)));
                    strncpy((viruses[indexOfline].mutations)[i].alt, tmpMut[i].alt, 1+strlen(tmpMut[i].alt));
                    free(tmpMut[i].alt);
                    tmpMut[i].alt = NULL;

                    (viruses[indexOfline].mutations)[i].mutationType = tmpMut[i].mutationType;
                }


#ifdef DEBUG
                //printf("num of feasible SNPs = %d\n", numOfFeasibleSnps);
                //showMut( viruses[indexOfline].mutations, viruses[indexOfline].numOfMut);
#endif
            }
            else{
                fprintf(stderr, "format of mutation file error!\n%s\ncolumnIndex = %d\nindexOfline = %d\nline = %s\n", token1, columnIndex, indexOfline, line);
                exit(1);
            }
            columnIndex ++;
        }
}

void readMutaions(const char * fin, bool * mutKeep, const int maxLengthOfSeq, int numOfViruses, struct Virus * viruses){
    //int lenthOfTmpMut = 1000;
    int lenthOfTmpMut = maxLengthOfSeq;
    struct Mutation * tmpMut = (struct Mutation *)malloc(sizeof(struct Mutation)*lenthOfTmpMut);
    char *str_ref = (char *)malloc(sizeof(char)*maxLengthOfSeq);
    char *str_alt = (char *)malloc(sizeof(char)*maxLengthOfSeq);

    FILE * fpin = fopen(fin, "r");
    if (NULL == fpin){
        printf("could not open file:\n%s\n", fin);
        exit(EXIT_FAILURE);
    }

    char *line = NULL;
    size_t len = 0;
    ssize_t read;

    int indexOfline = 0;//from 0
    while ((read = getline(&line, &len, fpin)) != -1) {//read line by line
        if(indexOfline >= numOfViruses){
            break;
        }
        chomp(line);
#ifdef DEBUG
        //printf("\n%zu\n", len);
        //printf("\n%s\n", line);
#endif
        readMutaionsFromLine(line, indexOfline, mutKeep, viruses, &lenthOfTmpMut, tmpMut, str_ref, str_alt);
        indexOfline ++;
    }


    for (int i = 0; i < numOfViruses; ++i) {
        if(NULL == viruses[i].mutations){
            viruses[i].numOfMut = 0;
        }
    }
    //close
    fclose(fpin);
    fpin = NULL;

    //free
    free(tmpMut);
    tmpMut = NULL;

    free(str_ref);
    str_ref = NULL;

    free(str_alt);
    str_alt = NULL;

    free(line);
    line = NULL;

#ifdef DEBUG
    printf("\n###READ MUTATIONS START###\n");
    printf("i\tname\tacc\tnumOfSnps\n");
    for (int i = 0; i < numOfViruses; ++i) {
        printf("%d\t%s\t%s\t%d\n", i, viruses[i].name, viruses[i].acc, viruses[i].numOfMut);
        showMut( viruses[i].mutations, viruses[i].numOfMut);
    }
    printf("###READ MUTATIONS END###\n");
#endif
}

typedef struct {
    int numOfViruses;
    bool* mutKeep;
    struct Virus* viruses;
    int* lenthOfTmpMut;
    struct Mutation* tmpMut;
    char* str_ref;
    char* str_alt;
} Args;

int lineHandler(const char* line, size_t len, size_t indexOfline, size_t numOfViruses, void* arg) {
    if (!line || !arg)
        return -1;
    Args* args = (Args*)arg;

    if (indexOfline == 0) {
        if (args->numOfViruses > numOfViruses)
            args->numOfViruses = numOfViruses;

        if (args->numOfViruses <= 0 || args->viruses)
            return -1;

        args->viruses = (struct Virus *)malloc(sizeof(struct Virus)*(args->numOfViruses));
        if (args->viruses == NULL) {
            fprintf(stderr, "Memory allocation failed for viruses\n");
            return -1;
        }
        virus_init(args->viruses, args->numOfViruses);
    }

    if (indexOfline >= args->numOfViruses)
        return 1;

    char* line_cp = malloc((len + 1) * sizeof(char));
    memcpy(line_cp, line, len + 1);
    readMutaionsFromLine(line_cp, indexOfline, args->mutKeep, args->viruses, args->lenthOfTmpMut, args->tmpMut, args->str_ref, args->str_alt);
    free(line_cp);

    return 0;
}

void readVCF(const char* str_vcf, const char* str_mut, bool * mutKeep, const int maxLengthOfSeq, int* numOfViruses, struct Virus ** viruses) {
    //int lenthOfTmpMut = 1000;
    int lenthOfTmpMut = maxLengthOfSeq;
    struct Mutation * tmpMut = (struct Mutation *)malloc(sizeof(struct Mutation)*lenthOfTmpMut);
    char *str_ref = (char *)malloc(sizeof(char)*maxLengthOfSeq);
    char *str_alt = (char *)malloc(sizeof(char)*maxLengthOfSeq);

    Args args = { *numOfViruses, mutKeep, NULL, &lenthOfTmpMut, tmpMut, str_ref, str_alt };
    if (vcf2mut(str_vcf, str_mut, lineHandler, (void*)&args, "*") < 0)
        exit(EXIT_FAILURE);
    *numOfViruses = args.numOfViruses;
    *viruses = args.viruses;

    for (int i = 0; i < *numOfViruses; ++i) {
        if(NULL == (*viruses)[i].mutations){
            (*viruses)[i].numOfMut = 0;
        }
    }

    //free
    free(tmpMut);
    tmpMut = NULL;

    free(str_ref);
    str_ref = NULL;

    free(str_alt);
    str_alt = NULL;
}






