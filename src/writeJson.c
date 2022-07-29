//
// Created by lun on 2/9/21.
//

//#define DEBUG

#define OLDVERSION //old node_0; new 0
#define NEWITERMS //jump
//#define JUMP_PIE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cJSON.h>
#include <stdbool.h>
#include <hashmap.h>
#include <time.h>

#include "virus.h"
#include "haplotype.h"
#include "edges.h"

void date2str(struct Date * date, char ** str);

//create a monitor with a list of supported resolutions
//NOTE: Returns a heap allocated string, you are required to free it after use.
int NumOfDigits(int n){
    int count = 0;

    if(n == 0)
    {
        return 0;
    }
    else
    {
        while(n)
        {
            n=n/10;
            count++;
        }
    }
    return count;
}

char * creat_network(struct Virus * viruses, int numOfViruses, struct Hap * haps, int numOfHaps, struct Edge * edges, int numOfEdges, int * network, int numOfEdgesInNet, int heightOfTree ){
#ifdef DEBUG
    printf("start creat network\n");
#endif
    char *string = NULL;
    cJSON *nodes = NULL;
    cJSON *node = NULL;
    cJSON *id = NULL;
    cJSON *radius = NULL;
    cJSON *jump = NULL;
    cJSON *SNPs = NULL;
    cJSON *SNP = NULL;

    cJSON *pieChart = NULL;
    cJSON *pie_cJSON = NULL;
    cJSON *color_cJSON = NULL;
    cJSON *percent_cJSON = NULL;
    cJSON *nodeDate_cJSON = NULL;

    cJSON *Viruses = NULL;
    cJSON *Virus = NULL;
    cJSON *acc = NULL;
    cJSON *date = NULL;
    cJSON *loci = NULL;
    cJSON *name = NULL;
    //infor
    cJSON *infor = NULL;
    //links
    cJSON *links = NULL;
    cJSON *link = NULL;
    cJSON *source = NULL;
    cJSON *target = NULL;
    cJSON *distance = NULL;
    cJSON *subset = NULL;
    cJSON *mintime_cJSON = NULL;

    //netinfor
    cJSON *netinfor = NULL;
    cJSON *netinfor_sub = NULL;

    //network
    cJSON *network_cJSON = cJSON_CreateObject();
    char * nodeIDstring = NULL;
    char * tmpstr_date =NULL;
    char * str_SNP = NULL;
    if (network_cJSON == NULL)
    {
        fprintf(stderr, "network_cJSON == NULL");
        goto end;
    }

    nodeIDstring = (char *)malloc(sizeof(char) * (1 + strlen("Node_") + NumOfDigits(numOfEdges)) );
    //nodes
    nodes = cJSON_CreateArray();
    if (nodes == NULL)
    {
        fprintf(stderr, "nodes == NULL");
        goto end;
    }
    cJSON_AddItemToObject(network_cJSON, "nodes", nodes);
    tmpstr_date = (char *)malloc(sizeof(char)*11);
    //str_SNP = (char *)malloc(sizeof(char) * (1+strlen("30000(Insertion:A->T)")));
    str_SNP = (char *)malloc(sizeof(char) * (1 + 30000 ));
    for (int i = 0; i < numOfHaps; ++i) {
        //nodes -- node
        node = cJSON_CreateObject();
        if (node == NULL)
        {
            fprintf(stderr, "node == NULL");
            goto end;
        }
        cJSON_AddItemToArray(nodes, node);

        //nodes -- node -- date//2021/4/8
        date2str(haps[i].date, &tmpstr_date);
        nodeDate_cJSON = cJSON_CreateString(tmpstr_date);
        if (nodeDate_cJSON == NULL)
        {
            fprintf(stderr, "fail to creat cjson object!\n");
            goto end;
        }
        cJSON_AddItemToObject(node, "date", nodeDate_cJSON);


        //nodes -- node -- id
#ifdef OLDVERSION
        sprintf(nodeIDstring, "Node_%d", i);
        id = cJSON_CreateString(nodeIDstring);
#else
        id = cJSON_CreateNumber(i);
#endif
        if (id == NULL)
        {
            fprintf(stderr, "id == NULL");
            goto end;
        }
        cJSON_AddItemToObject(node, "id", id);

#ifdef NEWITERMS
        //jump
        jump = cJSON_CreateNumber(haps[i].jump);
        if (jump == NULL)
        {
            fprintf(stderr, "jump == NULL");
            goto end;
        }
        cJSON_AddItemToObject(node, "jump", jump);
        //SNP
        SNPs = cJSON_CreateArray();
        //SNPs = cJSON_CreateObject();
        if (SNPs == NULL)
        {
            fprintf(stderr, "SNPs == NULL");
            goto end;
        }
        cJSON_AddItemToObject(node, "SNPs", SNPs);

        char tmp_mutType[20];
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
            sprintf(str_SNP, "%d(%s:%s->%s)", haps[i].mutations[j].pos, tmp_mutType, haps[i].mutations[j].ref, haps[i].mutations[j].alt);

            //sprintf(str_SNP, "%d(SNP:%s->%s)", haps[i].mutations[j].pos, haps[i].mutations[j].ref, haps[i].mutations[j].alt);
            SNP = cJSON_CreateString(str_SNP);
            if (SNP == NULL)
            {
                fprintf(stderr, "SNP == NULL");
                goto end;
            }
            cJSON_AddItemToArray(SNPs, SNP);
        }

#endif

        //nodes -- node -- radius
        radius = cJSON_CreateNumber(haps[i].numOfVirus);
        if (radius == NULL)
        {
            fprintf(stderr, "radius == NULL");
            goto end;
        }
        cJSON_AddItemToObject(node, "radius", radius);

        //nodes -- node -- pieChart
        pieChart = cJSON_CreateArray();
        if (pieChart == NULL)
        {
            fprintf(stderr, "fail to creat cjson object!\n");
            goto end;
        }
        cJSON_AddItemToObject(node, "pieChart", pieChart);
#ifdef JUMP_PIE
        sprintf(nodeIDstring, "%d", haps[i].jump);
        pie_cJSON = cJSON_CreateObject();
        if (pie_cJSON == NULL)
        {
            fprintf(stderr, "fail to creat cjson object!\n");
            goto end;
        }
        cJSON_AddItemToArray(pieChart, pie_cJSON);
        color_cJSON = cJSON_CreateString(nodeIDstring);
        if (color_cJSON == NULL)
        {
            fprintf(stderr, "fail to creat cjson object!\n");
            goto end;
        }
        cJSON_AddItemToObject(pie_cJSON, "color",color_cJSON);
        percent_cJSON = cJSON_CreateNumber(1);
        if (percent_cJSON == NULL)
        {
            fprintf(stderr, "fail to creat cjson object!\n");
            goto end;
        }
        cJSON_AddItemToObject(pie_cJSON, "percent",percent_cJSON);
#else
        const char *tmpcountry;
        int * value;
        hashmap_foreach(tmpcountry, value, &haps[i].country) {
        //for (int j = 0; j < 2; ++j) {
            pie_cJSON = cJSON_CreateObject();
            if (pie_cJSON == NULL)
            {
                fprintf(stderr, "fail to creat cjson object!\n");
                goto end;
            }
            cJSON_AddItemToArray(pieChart, pie_cJSON);
            //color_cJSON = cJSON_CreateString("color");
            color_cJSON = cJSON_CreateString(tmpcountry);
            if (color_cJSON == NULL)
            {
                fprintf(stderr, "fail to creat cjson object!\n");
                goto end;
            }
            cJSON_AddItemToObject(pie_cJSON, "color",color_cJSON);
            percent_cJSON = cJSON_CreateNumber((double)(*value)/(double)(haps[i].numOfVirus)*100.0);
            if (percent_cJSON == NULL)
            {
                fprintf(stderr, "fail to creat cjson object!\n");
                goto end;
            }
            cJSON_AddItemToObject(pie_cJSON, "percent",percent_cJSON);
        }
#endif

        //nodes -- node -- Viruses
        Viruses = cJSON_CreateArray();
        if (Viruses == NULL)
        {
            fprintf(stderr, "fail to creat cjson object!\n");
            goto end;
        }

        cJSON_AddItemToObject(node, "Virus", Viruses);

        for (int j = 0; j < haps[i].numOfVirus; ++j) {
            //nodes -- node -- Viruses -- Virus
            Virus = cJSON_CreateObject();
            if (Virus == NULL)
            {
                fprintf(stderr, "fail to creat cjson object!\n");
                goto end;
            }
            cJSON_AddItemToArray(Viruses, Virus);

            //nodes -- node -- Viruses -- Virus -- acc    (string)
            //acc = cJSON_CreateString("EPI_ISL_412972");

            acc = cJSON_CreateString(viruses[haps[i].allViruses[j]].acc);
            if (acc == NULL)
            {
                fprintf(stderr, "fail to creat cjson object!\n");
                goto end;
            }
            cJSON_AddItemToObject(Virus, "acc", acc);

            //nodes -- node -- Viruses -- Virus -- date   (string)
            //date = cJSON_CreateString("date");
            date2str(viruses[haps[i].allViruses[j]].date, &tmpstr_date);
            date = cJSON_CreateString(tmpstr_date);

            if (date == NULL)
            {
                fprintf(stderr, "fail to creat cjson object!\n");
                goto end;
            }
            cJSON_AddItemToObject(Virus, "date", date);

            //nodes -- node -- Viruses -- Virus -- loci   (string)
            //loci = cJSON_CreateString("Mexico");
            if(NULL != viruses[haps[i].allViruses[j]].location){
                if(NULL != viruses[haps[i].allViruses[j]].location->country){
                    loci = cJSON_CreateString(viruses[haps[i].allViruses[j]].location->country);
                    if (loci == NULL)
                    {
                        fprintf(stderr, "fail to creat cjson object!\n");
                        goto end;
                    }
                    cJSON_AddItemToObject(Virus, "loci", loci);
                }
            }

            //nodes -- node -- Viruses -- Virus -- name   (string)
            //name = cJSON_CreateString("BetaCoV/Mexico/CDMX/InDRE_01/2020");
            name = cJSON_CreateString(viruses[haps[i].allViruses[j]].name);
            if (name == NULL)
            {
                fprintf(stderr, "fail to creat cjson object!\n");
                goto end;
            }
            cJSON_AddItemToObject(Virus, "name", name);
        }

    }

    //infor
    char * str_infor = (char *)malloc(sizeof(char) * 256);
    time_t currentTime;
    struct tm * currentTM;

    time(&currentTime);
    currentTM = localtime(&currentTime);
    sprintf(str_infor, "As of %04d-%02d-%02d %02d:%02d:%02d, %d strains of virus have been sampled, with a total of %d haplotypes. heightOfTree = %d", currentTM->tm_year+1900, currentTM->tm_mon+1, currentTM->tm_mday, currentTM->tm_hour, currentTM->tm_min, currentTM->tm_sec, numOfViruses, numOfHaps, heightOfTree);
    infor = cJSON_CreateString(str_infor);
    //infor = cJSON_CreateString("As of 2020-10-15, 3237 strains of virus have been sampled, with a total of 207 haplotypes.");
    free(str_infor);
    str_infor = NULL;
    if (infor == NULL)
    {
        fprintf(stderr, "fail to creat cjson object!\n");
        goto end;
    }
    cJSON_AddItemToObject(network_cJSON, "infor", infor);

//links
    links = cJSON_CreateArray();
    if (links == NULL)
    {
        fprintf(stderr, "fail to creat cjson object!\n");
        goto end;
    }
    cJSON_AddItemToObject(network_cJSON, "links", links);

    for (int i = 0; i < numOfEdgesInNet; ++i) {
        //links -- link
        link = cJSON_CreateObject();
        if (link == NULL)
        {
            fprintf(stderr, "fail to creat cjson object!\n");
            goto end;
        }
        cJSON_AddItemToArray(links, link);

        //links -- link -- source
#ifdef OLDVERSION
        sprintf(nodeIDstring, "Node_%d", edges[network[i]].source);
        source = cJSON_CreateString(nodeIDstring);
#else
        source = cJSON_CreateNumber(edges[network[i]].source);
#endif
        if (source == NULL)
        {
            fprintf(stderr, "fail to creat cjson object!\n");
            goto end;
        }
        cJSON_AddItemToObject(link, "source", source);

        //links -- link -- target
#ifdef OLDVERSION
        sprintf(nodeIDstring, "Node_%d", edges[network[i]].target);
        target = cJSON_CreateString(nodeIDstring);
#else
        target = cJSON_CreateNumber(edges[network[i]].target);
#endif
        if (target == NULL)
        {
            fprintf(stderr, "fail to creat cjson object!\n");
            goto end;
        }
        cJSON_AddItemToObject(link, "target", target);

        //links -- link -- distance
        distance = cJSON_CreateNumber(edges[network[i]].distance);
        if (distance == NULL)
        {
            fprintf(stderr, "fail to creat cjson object!\n");
            goto end;
        }
        cJSON_AddItemToObject(link, "distance", distance);

        //subset
        //links -- link -- subset
        subset = cJSON_CreateNumber(edges[network[i]].subset);
        if (subset == NULL)
        {
            fprintf(stderr, "fail to creat cjson object!\n");
            goto end;
        }
        cJSON_AddItemToObject(link, "subset", subset);

        //links -- link -- timemin
        mintime_cJSON = cJSON_CreateNumber(edges[network[i]].mintime);
        if (mintime_cJSON == NULL)
        {
            fprintf(stderr, "fail to creat cjson object!\n");
            goto end;
        }
        cJSON_AddItemToObject(link, "minTime", mintime_cJSON);

    }

    //net infor
    netinfor = cJSON_CreateObject();
    if (nodes == NULL)
    {
        fprintf(stderr, "fail to creat cjson object!\n");
        goto end;
    }
    cJSON_AddItemToObject(network_cJSON, "netInfor", netinfor);
    //numOfViruses
    netinfor_sub = cJSON_CreateNumber(numOfViruses);
    cJSON_AddItemToObject(netinfor, "numOfViruses", netinfor_sub);
    //numOfEdgesInNet
    netinfor_sub = cJSON_CreateNumber(numOfEdgesInNet);
    cJSON_AddItemToObject(netinfor, "numOfEdges", netinfor_sub);
    //numOfHaps
    netinfor_sub = cJSON_CreateNumber(numOfHaps);
    cJSON_AddItemToObject(netinfor, "numOfHaps", netinfor_sub);
    //heightOfTree
    netinfor_sub = cJSON_CreateNumber(heightOfTree);
    cJSON_AddItemToObject(netinfor, "heightOfTree", netinfor_sub);

    //print
    string = cJSON_Print(network_cJSON);
    if (string == NULL)
    {
        fprintf(stderr, "Failed to print network.\n");
    }
#ifdef DEBUG
    printf("%s\n", string);
#endif
    //end
    end:
    if(NULL != network_cJSON) {
        cJSON_Delete(network_cJSON);
    }
    else {
        exit(1);
    }

    //free
    free(nodeIDstring);
    nodeIDstring = NULL;
    free(tmpstr_date);
    tmpstr_date = NULL;
    free(str_SNP);
    str_SNP = NULL;
    return string;
}