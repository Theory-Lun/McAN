//
// Created by lun on 4/13/21.
//

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
//#include <unistd.h>
#include <time.h>
//#include <sys/sysinfo.h>

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
void date2str(struct Date * date, char ** str);

bool checkUsageAll(const char * str_vcf, const char * str_mut, const char * str_sta_anno, const char * str_meta, const char * str_minFreq, const char * str_maxvirus, const char * str_folder,  double * minFreq, int * maxvirus );
#ifdef MUTILTHREADING
void calcDistDirected_sorthap(struct Hap * haps, int numOfHaps, struct Edge * edges, int * numOfEdges, int nThread);
#else
void calcDistDirected_sorthap(struct Hap * haps, int numOfHaps, struct Edge * edges, int * numOfEdges );
#endif
void calcDisAll(int *network, int numOfEdgesInNet, struct Edge * edges, int *disAll);
void calcSrcSizeAve(int *network, int numOfEdgesInNet, struct Edge *edges, double *aveOfSourceSize);
void calcNumOfEdgesTimeGood(int *network, int numOfEdgesInNet, struct Edge *edges, int *numOfEdgesTimeGood);
void findHap_hashmap(struct Virus * viruses, int numOfViruses, struct Hap ** haps, int * numOfHaps);
void free_haps2(struct Hap ** pHaps, int numOfHaps);
typedef int(*Callback)(const char*, size_t, size_t, size_t, void*);
int vcf2mut(const char* vcfName, const char* mutName, Callback lineHandler, void* arg, const char* missingValue);
void readVCF(const char* str_vcf, const char* str_mut, bool * mutKeep, const int maxLengthOfSeq, int* numOfViruses, struct Virus ** viruses);
int outputGraphML(const char* str_folder, struct Virus* viruses, int numOfViruses, struct Hap* haps, int numOfHaps, struct Edge* edges, int* network, int numOfEdgesInNet);

#ifdef MUTILTHREADING
void hapNetDirected(const char * str_vcf, const char * str_mut, const char * str_meta, const char * str_minFreq, const char * str_maxvirus, const char * str_folder, char *str_sitemask, int* outFrmt, int nThread) {
#else
void hapNetDirected(const char * str_vcf, const char * str_mut, const char * str_meta, const char * str_minFreq, const char * str_maxvirus, const char * str_folder, char *str_sitemask, int* outFrmt) {
#endif
    // const char * str_mut, input file mutations
    // const char * str_sta_anno, input file states [optional]
    // const char * str_meta, input file metadata
    // const char * str_minFreq, filtered out SNPs with frequency < minFreq [optional]
    // const char * str_maxvirus, maximum number of virus wanted to read [optional]
    // const char * str_folder, output folder

    //for memory
    //struct sysinfo s_info;
    //long freeram_max;
    //long freeram_min;
    //long freeswap_max;
    //long freeswap_min;
    //long memoryUsed;

    //if(sysinfo(&s_info)==0){
    //    freeram_max = s_info.freeram;
    //    freeswap_max = s_info.freeswap;
    //}

    clock_t start_all, end_all;
    double duration_all;
    start_all = clock();

#ifdef DEBUG_TIME
    clock_t start, end;
    double duration;
#endif

    //some const here
    const int maxLengthOfSeq = 30000;

    //check usage and read numbers
    double minFreq = 0;
    int maxvirus = 0;
    if(false == checkUsageAll(str_vcf, str_mut, NULL, str_meta, str_minFreq, str_maxvirus, str_folder, &minFreq, &maxvirus ) ){
        //help();
        exit(1);
    }

    char* str_mut_out = NULL;
    if (outFrmt[3]) {
        size_t l = strlen(str_folder) + strlen("mut") + 2;//2 for '\0' and '/'
        str_mut_out = (char *)malloc(sizeof(char) * l);
        if (NULL == str_mut_out) {
            fprintf(stderr, "Memory allocation failed for str_mut_out\n");
            exit(1);
        }
        strcpy(str_mut_out, str_folder);
        strcat(str_mut_out, "/");
        strcat(str_mut_out, "mut");
    }

    if (str_vcf && outFrmt[0] + outFrmt[1] + outFrmt[2] == 0) {
        if (-1 == vcf2mut(str_vcf, str_mut_out, NULL, NULL, "*"))
            exit(1);
        return;
    }

    int numOfViruses = 0;

    if (str_vcf) {
        numOfViruses = INT_MAX;
    }

    if (str_mut) {
        //read numOfViruses
        readNumOfVirusesFromMut(str_mut, &numOfViruses);
    }

    //mutation filter
    //mutKeep[i]: 0 (false), ignore pos i+1; 1 (true), keep pos i+1
    bool * mutKeep = (bool *)malloc(sizeof(int)*maxLengthOfSeq);
    if(mutKeep == NULL){
        fprintf(stderr, "Memory allocation failed for mutKeep\n");
        return;
    }
#ifdef DEBUG_TIME
    start = clock();
#endif
    //read site kept
    {
        for (int i = 0; i < maxLengthOfSeq; ++i) {
            mutKeep[i] = false;
        }
        FILE *fpin = fopen(str_sitemask, "r");
        if(NULL == fpin){
            fprintf(stderr, "\nError: failed to open file: %s\n", str_sitemask);
            exit(1);
        }
        int num=1;
        while(num > 0){
            int tmp;
            num = fscanf(fpin, "%d\n", &tmp);
            if(num > 0){
                if(tmp <= maxLengthOfSeq && tmp >= 1){
                    mutKeep[tmp-1] = true;
                }
            }
        }
        fclose(fpin);
        fpin = NULL;
    }

    //mutationFilter(freqEachPos, maxLengthOfSeq, minFreq, numOfViruses, mutKeep);
#ifdef DEBUG_TIME
    end = clock();
    duration = -((double)(start - end)) / CLOCKS_PER_SEC;
    printf("%f s for constructing mutation filter\n", duration);
#endif
    //reduce numOfViruses for test
    if(NULL != str_maxvirus && maxvirus > 0 && maxvirus < numOfViruses){
        numOfViruses = maxvirus;
    }
#ifdef DEBUG
    printf("numOfViruses = %d\n", numOfViruses);
#endif

    struct Virus * viruses = NULL;

    if (str_vcf) {
#ifdef DEBUG_TIME
        start = clock();
#endif
        readVCF(str_vcf, str_mut_out, mutKeep, maxLengthOfSeq, &numOfViruses, &viruses);
#ifdef DEBUG_TIME
        end = clock();
        duration = -((double)(start - end)) / CLOCKS_PER_SEC;
        printf("%f s for reading mutations from vcf\n", duration);
#endif
    }

    if (str_mut) {
        //read mutations
        //only keep Snps and consider mutKeep
        viruses = (struct Virus *)malloc(sizeof(struct Virus)*numOfViruses);
        if (viruses == NULL) {
            fprintf(stderr, "Memory allocation failed for viruses\n");
            return;
        }
        virus_init(viruses, numOfViruses);
#ifdef DEBUG_TIME
        start = clock();
#endif
        readMutaions(str_mut, mutKeep, maxLengthOfSeq, numOfViruses, viruses);
#ifdef DEBUG_TIME
        end = clock();
        duration = -((double)(start - end)) / CLOCKS_PER_SEC;
        printf("%f s for reading mutations\n", duration);
#endif
    }

    //read metadata
#ifdef DEBUG_TIME
    start = clock();
#endif
    readMeta_hash(str_meta, numOfViruses, viruses);
    //readMeta_3col_hash(str_meta, numOfViruses, viruses);
    //readMetaLine(str_meta, numOfViruses, viruses);
#ifdef DEBUG_TIME
    end = clock();
    duration = -((double)(start - end)) / CLOCKS_PER_SEC;
    printf("%f s for reading metadata and matching with mutations\n", duration);
#endif

#ifdef DEBUG
    showVirus(viruses, numOfViruses);
    //getchar();
#endif

    //sort time, needed for calc time in each haplotype
#ifdef DEBUG_TIME
    start = clock();
#endif
    qsort(viruses, numOfViruses, sizeof(struct Virus), cmp_time_virus);
#ifdef DEBUG_TIME
    end = clock();

    duration = -((double)(start - end)) / CLOCKS_PER_SEC;
    printf("%f s for sort all viruses (time)\n", duration);
#endif


#ifdef DEBUG
    for (int i = 0; i < numOfViruses; ++i) {
        printf("%04d %02d %02d\n", viruses[i].date->year, viruses[i].date->month,viruses[i].date->day);
    }
#endif
    int numOfVirusesALL = numOfViruses;

    //find haplotype
    struct Hap * haps = NULL;
    int numOfHaps = -1;
#ifdef DEBUG_TIME
    start = clock();
#endif
    findHap_hashmap(viruses, numOfViruses, &haps, &numOfHaps);

#ifdef DEBUG_TIME
    end = clock();
    duration = -((double)(start - end)) / CLOCKS_PER_SEC;
    printf("%f s for finding haplotypes\n", duration);
    printf("virus = %d, hap = %d\n", numOfViruses, numOfHaps);
#endif
#ifdef DEBUG
    hap_show(haps, numOfHaps);
#endif

    //calculate distance (find ancestor)
    int numOfEdges = numOfHaps;
    struct Edge * edges = (struct Edge *)malloc(sizeof(struct Edge) * numOfEdges);
    if(edges == NULL){
        fprintf(stderr, "Memory allocation failed for edges\n");
        return;
    }
#ifdef DEBUG_TIME
    start = clock();
#endif
    //calcDistDirected(haps, numOfHaps, edges, &numOfEdges );
#ifdef MUTILTHREADING
    calcDistDirected_sorthap(haps, numOfHaps, edges, &numOfEdges, nThread);
#else
    calcDistDirected_sorthap(haps, numOfHaps, edges, &numOfEdges );
#endif
    //calcDistDirectedNew(haps, numOfHaps, edges, numOfEdges);
#ifdef DEBUG_TIME
    end = clock();
    duration = -((double)(start - end)) / CLOCKS_PER_SEC;
    printf("%f s for calculating distances\n", duration);
#endif

#ifdef DEBUG
    edges_show(edges, numOfEdges);
    //getchar();
#endif

//construct hap net
#ifdef DEBUG_TIME
    start = clock();
#endif
    int numOfEdgesInNet = numOfHaps - 1;
    int * network = (int *)malloc(sizeof(int) * numOfEdgesInNet);
    if(network == NULL){
        fprintf(stderr, "Memory allocation failed for network\n");
        return;
    }
    int count = 0;
    for (int i = 0; i < numOfEdges; ++i) {
        if( 0 != edges[i].pHap2->numOfMut ){
            network[count] = i; //index of useful edges
            count ++;
        }
    }
    if(count != numOfEdges-1){
        fprintf(stderr, "Error: count = %d, numOfHaps = %d", count, numOfHaps);
        exit(1);
    }
#ifdef DEBUG_TIME
    end = clock();
    duration = -((double)(start - end)) / CLOCKS_PER_SEC;
    printf("%f s for constructing hap net\n", duration);
#endif

    //calculate country percentage for each haplotype
#ifdef DEBUG_TIME
    start = clock();
#endif
    countryDistribution(haps, numOfHaps, viruses, numOfViruses);
#ifdef DEBUG_TIME
    end = clock();
    duration = -((double)(start - end)) / CLOCKS_PER_SEC;
    printf("%f s for calculating country percentage for each haplotype\n", duration);
#endif

    //analysis
    start = clock();
    int heightOfTree;
    treeHeight(haps, numOfHaps, &heightOfTree);
    printf("heightOfTree = %d\n", heightOfTree);
    //calculate cost of the network (sum of distance, sum of number of samples in source of each edge, number of good edge according to time criterion)


    int disAll;
    calcDisAll(network, numOfEdgesInNet, edges, &disAll);
    printf("sumOfDist = %d\n", disAll);

    double aveOfSourceSize;
    calcSrcSizeAve(network, numOfEdgesInNet, edges, &aveOfSourceSize);
    printf("aveSourceSize = %f\n", aveOfSourceSize);

    int numOfEdgesTimeGood;
    calcNumOfEdgesTimeGood(network, numOfEdgesInNet, edges, &numOfEdgesTimeGood);
    printf("numOfEdgesTimeGood = %d\n", numOfEdgesTimeGood);

    end = clock();
    duration = -((double)(start - end)) / CLOCKS_PER_SEC;
    printf("%f s for test\n", duration);
    fflush(stdout);

    //creat json
    char * json = NULL;
    if (outFrmt[0]) {
#ifdef DEBUG_TIME
        start = clock();
#endif
        json = creat_network(viruses, numOfViruses, haps, numOfHaps, edges, numOfEdges, network, numOfEdgesInNet, heightOfTree);

#ifdef DEBUG_TIME
        end = clock();
        duration = -((double)(start - end)) / CLOCKS_PER_SEC;
        printf("%f s for creating json\n", duration);
#endif
#ifdef DEBUG
        printf("###network###:\n%s\n", json);
#endif
    }

    //for memory
    //if(sysinfo(&s_info)==0){
    //    freeram_min = s_info.freeram;
    //    freeswap_min = s_info.freeswap;
    //}
    //memoryUsed = freeram_max +  freeswap_max - freeram_min -  freeswap_min;// Byte

#ifdef DEBUG
    system("sh /home/lun/work/experiments/test-mem/mem-my.sh HapNet_v0_2_5");
#endif
    //all time
    end_all = clock();
    duration_all = -((double)(start_all - end_all)) / CLOCKS_PER_SEC;

    //output json
    if (outFrmt[0]) {
#ifdef DEBUG_TIME
        start = clock();
#endif
        output(str_folder, json, /*memoryUsed*/0, duration_all);
        free(json);
        json = NULL;
#ifdef DEBUG_TIME
        end = clock();
        duration = -((double)(start - end)) / CLOCKS_PER_SEC;
        printf("%f s for writing results in JSON format\n", duration);
#endif
    }

    //output tsv
    if (outFrmt[1]) {
#ifdef DEBUG_TIME
        start = clock();
#endif

    {
        //input str_folder, duration_all, edges, haplotype, network
        //edges
        char *str_edges = NULL;
        size_t l = strlen(str_folder) + strlen("Edges_out.tsv") + 2;//2 for '\0' and '/'
        str_edges = (char *) malloc(sizeof(char) * l);
        if (NULL == str_edges) {
            fprintf(stderr, "Memory allocation failed for str_edges\n");
            exit(1);
        }
        strcpy(str_edges, str_folder);
        strcat(str_edges, "/");
        strcat(str_edges, "Edges_out.tsv");
        FILE *fpout = fopen(str_edges, "w");
        if (NULL == fpout) {
            printf("open file error!\n");
            exit(EXIT_FAILURE);
        }
        for (int i = 0; i < numOfEdgesInNet; ++i) {
            fprintf(fpout, "Node_%d\t", edges[network[i]].source);
            fprintf(fpout, "Node_%d\t", edges[network[i]].target);
            fprintf(fpout, "%d\n", edges[network[i]].distance);
        }
        fclose(fpout);
        fpout = NULL;
        free(str_edges);
        str_edges = NULL;
    }
    {
        //haplotype
        char *str_haps = NULL;
        size_t l = strlen(str_folder) + strlen("Hap_out.tsv") + 2;//2 for '\0' and '/'
        str_haps = (char *) malloc(sizeof(char) * l);
        if (NULL == str_haps) {
            fprintf(stderr, "Memory allocation failed for str_haps\n");
            exit(1);
        }
        strcpy(str_haps, str_folder);
        strcat(str_haps, "/");
        strcat(str_haps, "Hap_out.tsv");
        FILE *fpout = fopen(str_haps, "w");
        if (NULL == fpout) {
            printf("open file error!\n");
            exit(EXIT_FAILURE);
        }
        for (int i = 0; i < numOfHaps; ++i) {
            fprintf(fpout, "Node_%d\t", i);//node ID
            for (int j = 0; j < haps[i].numOfVirus; ++j) {
                if (j!=0) {
                    fprintf(fpout, ";");
                }
                fprintf(fpout, "%s",viruses[haps[i].allViruses[j]].acc);
            }
            fprintf(fpout, "\t");
            char tmp_mutType[20];
            if (haps[i].numOfMut==0) {
                fprintf(fpout, "N/A");
            }
            for (int j = 0; j < haps[i].numOfMut; ++j) {
                if(MUTATION_INDEL == haps[i].mutations[j].mutationType ){
                    strcpy(tmp_mutType, "Indel");
                }
                else if(MUTATION_DELETION == haps[i].mutations[j].mutationType ){
                    strcpy(tmp_mutType, "Deletion");
                }
                else if(MUTATION_INSERTION == haps[i].mutations[j].mutationType ){
                    strcpy(tmp_mutType, "Insertion");
                }
                else if(MUTATION_SNP == haps[i].mutations[j].mutationType ){
                    strcpy(tmp_mutType, "SNP");
                }
                if (j!=0) {
                    fprintf(fpout, ";");
                }
                fprintf(fpout, "%d(%s:%s->%s)", haps[i].mutations[j].pos, tmp_mutType, haps[i].mutations[j].ref, haps[i].mutations[j].alt);
            }
            fprintf(fpout, "\n");
        }
        fclose(fpout);
        fpout = NULL;
        free(str_haps);
        str_haps = NULL;
    }

#ifdef DEBUG_TIME
        end = clock();
        duration = -((double)(start - end)) / CLOCKS_PER_SEC;
        printf("%f s for writing results in TSV format\n", duration);
#endif
    }

    //output GraphML
    if (outFrmt[2]) {
#ifdef DEBUG_TIME
        start = clock();
#endif

        if (outputGraphML(str_folder, viruses, numOfViruses, haps, numOfHaps, edges, network, numOfEdgesInNet))
            exit(EXIT_FAILURE);

#ifdef DEBUG_TIME
        end = clock();
        duration = -((double)(start - end)) / CLOCKS_PER_SEC;
        printf("%f s for writing results in GraphML format\n", duration);
#endif
    }

    start = clock();
    //free, need to free hap.*, ...,

    free(mutKeep);
    mutKeep = NULL;
    free(edges);
    edges = NULL;
    free(str_mut_out);
    str_mut_out = NULL;


    free_haps2(&haps, numOfHaps);


    free_viruses(&viruses, numOfVirusesALL);
    //free(json);
    //free(json);//!!!!!!!!!!!!!malloc_consolidate(): invalid chunk size when numOfSample is smaller than 100
    //json = NULL;
    end = clock();
    duration = -((double)(start - end)) / CLOCKS_PER_SEC;
    printf("%f s for free\n", duration);
}

int outputGraphML(const char* str_folder, struct Virus* viruses, int numOfViruses, struct Hap* haps, int numOfHaps, struct Edge* edges, int* network, int numOfEdgesInNet) {
    char *str_graphml = NULL;
    size_t l = strlen(str_folder) + strlen("haplotype_loci.graphml") + 2;//2 for '\0' and '/'
    str_graphml = (char *)malloc(sizeof(char) * l);
    if (NULL == str_graphml) {
        fprintf(stderr, "Memory allocation failed for str_graphml\n");
        return -1;
    }
    strcpy(str_graphml, str_folder);
    strcat(str_graphml, "/");
    strcat(str_graphml, "haplotype_loci.graphml");
    FILE *fpout = fopen(str_graphml, "w");
    if (NULL == fpout) {
        printf("open file error!\n");
        return -1;
    }

    fprintf(fpout, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
    fprintf(fpout, "<!-- Created by McAN -->\n");
    fprintf(fpout, "<graphml xmlns=\"http://graphml.graphdrawing.org/xmlns\"\n");
    fprintf(fpout, "         xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n");
    fprintf(fpout, "         xsi:schemaLocation=\"http://graphml.graphdrawing.org/xmlns\n");
    fprintf(fpout, "                             http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd\">\n");

    fprintf(fpout, "  <key id=\"g_numOfViruses\" for=\"graph\" attr.name=\"numOfViruses\" attr.type=\"int\"/>\n");
    fprintf(fpout, "  <key id=\"g_numOfEdges\" for=\"graph\" attr.name=\"numOfEdges\" attr.type=\"int\"/>\n");
    fprintf(fpout, "  <key id=\"g_numOfHaps\" for=\"graph\" attr.name=\"numOfHaps\" attr.type=\"int\"/>\n");
    fprintf(fpout, "  <key id=\"v_date\" for=\"node\" attr.name=\"date\" attr.type=\"string\"/>\n");
    fprintf(fpout, "  <key id=\"v_numOfMut\" for=\"node\" attr.name=\"numOfMut\" attr.type=\"int\"/>\n");
    fprintf(fpout, "  <key id=\"v_mutation\" for=\"node\" attr.name=\"mutation\" attr.type=\"string\"/>\n");
    fprintf(fpout, "  <key id=\"v_numOfVirus\" for=\"node\" attr.name=\"numOfVirus\" attr.type=\"int\"/>\n");
    fprintf(fpout, "  <key id=\"v_virus_acc\" for=\"node\" attr.name=\"virus_acc\" attr.type=\"string\"/>\n");
    fprintf(fpout, "  <key id=\"v_virus_date\" for=\"node\" attr.name=\"virus_date\" attr.type=\"string\"/>\n");
    fprintf(fpout, "  <key id=\"v_virus_loc_1\" for=\"node\" attr.name=\"virus_loc_1\" attr.type=\"string\"/>\n");
    fprintf(fpout, "  <key id=\"v_virus_loc_2\" for=\"node\" attr.name=\"virus_loc_2\" attr.type=\"string\"/>\n");
    fprintf(fpout, "  <key id=\"v_virus_loc_3\" for=\"node\" attr.name=\"virus_loc_3\" attr.type=\"string\"/>\n");
    fprintf(fpout, "  <key id=\"v_virus_name\" for=\"node\" attr.name=\"virus_name\" attr.type=\"string\"/>\n");
    fprintf(fpout, "  <key id=\"e_distance\" for=\"edge\" attr.name=\"distance\" attr.type=\"double\"/>\n");

    fprintf(fpout, "  <graph id=\"G\" edgedefault=\"directed\"\n");
    fprintf(fpout, "                parse.nodes=\"%d\" parse.edges=\"%d\"\n", numOfHaps, numOfEdgesInNet);
    fprintf(fpout, "                parse.nodeids=\"canonical\" parse.edgeids=\"canonical\"\n");
    fprintf(fpout, "                parse.order=\"nodesfirst\">\n");

    fprintf(fpout, "    <data key=\"g_numOfViruses\">%d</data>\n", numOfViruses);
    fprintf(fpout, "    <data key=\"g_numOfEdges\">%d</data>\n", numOfEdgesInNet);
    fprintf(fpout, "    <data key=\"g_numOfHaps\">%d</data>\n", numOfHaps);

    char tmpstr_date[11];
    char* ptr = tmpstr_date;
    tmpstr_date[10] = 0;

    for (int i = 0; i < numOfHaps; ++i) {
        fprintf(fpout, "    <node id=\"n%d\">\n", i);

        fprintf(fpout, "      <data key=\"v_date\">");
        if (haps[i].date) {
            date2str(haps[i].date, &ptr);
            fprintf(fpout, "%s", tmpstr_date);
        }
        else
            fprintf(fpout, "N/A");
        fprintf(fpout, "</data>\n");

        fprintf(fpout, "      <data key=\"v_numOfMut\">%d</data>\n", haps[i].numOfMut);

        fprintf(fpout, "      <data key=\"v_mutation\">");
        static const char* mutTypeStr[] = { "Indel", "Deletion", "Insertion", "SNP" };
        if (haps[i].numOfMut == 0)
            fprintf(fpout, "N/A");
        for (int j = 0; j < haps[i].numOfMut; ++j) {
            if (j != 0)
                fprintf(fpout, ";");
            fprintf(fpout, "%d(%s:%s-&gt;%s)", haps[i].mutations[j].pos, mutTypeStr[haps[i].mutations[j].mutationType], haps[i].mutations[j].ref, haps[i].mutations[j].alt);
        }
        fprintf(fpout, "</data>\n");

        fprintf(fpout, "      <data key=\"v_numOfVirus\">%d</data>\n", haps[i].numOfVirus);

        fprintf(fpout, "      <data key=\"v_virus_acc\">");
        for (int j = 0; j < haps[i].numOfVirus; ++j) {
            if (j != 0)
                fprintf(fpout, ";");
            fprintf(fpout, "%s", viruses[haps[i].allViruses[j]].acc);
        }
        fprintf(fpout, "</data>\n");

        fprintf(fpout, "      <data key=\"v_virus_date\">");
        for (int j = 0; j < haps[i].numOfVirus; ++j) {
            if (j != 0)
                fprintf(fpout, ";");
            if (viruses[haps[i].allViruses[j]].date) {
                date2str(viruses[haps[i].allViruses[j]].date, &ptr);
                fprintf(fpout, "%s", tmpstr_date);
            }
            else
                fprintf(fpout, "N/A");
        }
        fprintf(fpout, "</data>\n");

        fprintf(fpout, "      <data key=\"v_virus_loc_1\">");
        for (int j = 0; j < haps[i].numOfVirus; ++j) {
            if (j != 0)
                fprintf(fpout, ";");
            if (viruses[haps[i].allViruses[j]].location->country)
                fprintf(fpout, "%s", viruses[haps[i].allViruses[j]].location->country);
            else
                fprintf(fpout, "N/A");
        }
        fprintf(fpout, "</data>\n");

        fprintf(fpout, "      <data key=\"v_virus_loc_2\">");
        for (int j = 0; j < haps[i].numOfVirus; ++j) {
            if (j != 0)
                fprintf(fpout, ";");
            if (viruses[haps[i].allViruses[j]].location->state)
                fprintf(fpout, "%s", viruses[haps[i].allViruses[j]].location->state);
            else
                fprintf(fpout, "N/A");
        }
        fprintf(fpout, "</data>\n");

        fprintf(fpout, "      <data key=\"v_virus_loc_3\">");
        for (int j = 0; j < haps[i].numOfVirus; ++j) {
            if (j != 0)
                fprintf(fpout, ";");
            if (viruses[haps[i].allViruses[j]].location->city)
                fprintf(fpout, "%s", viruses[haps[i].allViruses[j]].location->city);
            else
                fprintf(fpout, "N/A");
        }
        fprintf(fpout, "</data>\n");

        fprintf(fpout, "      <data key=\"v_virus_name\">");
        for (int j = 0; j < haps[i].numOfVirus; ++j) {
            if (j != 0)
                fprintf(fpout, ";");
            if (viruses[haps[i].allViruses[j]].name && viruses[haps[i].allViruses[j]].name[0])
                fprintf(fpout, "%s", viruses[haps[i].allViruses[j]].name);
            else
                fprintf(fpout, "N/A");
        }
        fprintf(fpout, "</data>\n");

        fprintf(fpout, "    </node>\n");
    }

    for (int i = 0; i < numOfEdgesInNet; ++i) {
        fprintf(fpout, "    <edge id=\"e%d\" source=\"n%d\" target=\"n%d\">\n", i, edges[network[i]].source, edges[network[i]].target);
        fprintf(fpout, "      <data key=\"e_distance\">%d</data>\n", edges[network[i]].distance);
        fprintf(fpout, "    </edge>\n");
    }

    fprintf(fpout, "  </graph>\n");
    fprintf(fpout, "</graphml>\n");

    fclose(fpout);
    fpout = NULL;
    free(str_graphml);
    str_graphml = NULL;

    return 0;
}
