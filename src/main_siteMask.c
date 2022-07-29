//
// Created by lun on 2021/5/27.
//


#include <stdio.h>
#include <stdlib.h>
//#include <getopt.h>
#include <string.h>
#include <stdbool.h>
#include "Hashmap_acc2line.h"
#include "virus.h"
#include "optionparser.h"

//#define DEBUG
void NextMonth_struct(struct Date *structDate);
void readMut_acc2line_hash(char *str_input, struct Hashmap_acc2line *hash);
void readMeta_acc2line_hash(char *str_input, struct Hashmap_acc2line *hash);
void extractDateFromMetaLine(char *line, struct Date *pDate);
void readMeta_acc2Date_hash(char *str_input, struct Hashmap_acc2Date *hash);
int cmp_time_date(const void* _a, const void* _b);
void show_date(struct Date * date);
bool chomp( char *str );
char* strip(char* str);
void readVCF_acc2line_hash(const char* str_vcf, const char* str_mut, struct Hashmap_acc2line *hash);
void version();

void calcMask(int * freqEachPos, int maxLengthOfSeq, double minFreq, bool * mutKeep, int numOfSamples){
    double minNum = minFreq * numOfSamples;
    for (int i = 0; i < maxLengthOfSeq; ++i) {
        if(0 == freqEachPos[i]){ // freq = 0 or UTR
            mutKeep[i] = false;
        }
        else if( minNum > freqEachPos[i]){
            mutKeep[i] = false;
        }
        else{
            mutKeep[i] = true;
        }
    }
}

void siteMask(char *str_inputVCF, char *str_inputMut, char *str_inputMeta, char *str_output, char *str_minFreq, bool mutOutput){
    //minfreq
    double minFreq;

    if(NULL != str_minFreq){
        char *endptr = NULL; // for strtol
        minFreq = strtod(str_minFreq, &endptr);
        if('\0' != endptr[0] || '\0' == str_minFreq[0]){
            printf("error: minfreq is not a real number:\n%s\n", str_minFreq);
            exit(1);
        }
        endptr = NULL;

        if(minFreq < 0 || minFreq >= 1.0){
            printf("error: minFreq = %f\n", minFreq);
            printf("minFreq must belong to [0,1)\n");
            exit(1);
        }
    }
    else{
        fprintf(stderr, "no minimum frequency\n");
        exit(1);
    }

    printf("minimum frequency = %f\n", minFreq);


    char* str_mut_out = NULL;
    if (mutOutput) {
        char* ptr = strrchr(str_output, '/');
        size_t n = (ptr) ? ptr - str_output + 1 : 0;
        size_t l = n + strlen("mut") + 1;//1 for '\0'
        str_mut_out = (char *)malloc(sizeof(char) * l);
        if (NULL == str_mut_out) {
            exit(1);
        }
        strncpy(str_mut_out, str_output, n);
        strcat(str_mut_out, "mut");
    }


    //read Mut and Meta
    printf("read Mutation and Metadata ...\n");
    struct Hashmap_acc2line Mut;
    if (str_inputVCF) {
        readVCF_acc2line_hash(str_inputVCF, str_mut_out, &Mut);
    }
    if (str_inputMut) {
        readMut_acc2line_hash(str_inputMut, &Mut);
    }
    struct Hashmap_acc2Date acc2date;
    readMeta_acc2Date_hash(str_inputMeta, &acc2date);
    //readMeta_3col_acc2Date_hash(str_inputMeta, &acc2date);
    printf("done\n");

    const char *acc;
    char *line;
    struct Date *structDate;
    //min max
    struct Date minStructDate, maxStructDate;
    int count = 0;
    hashmap_foreach(acc, structDate, &(acc2date.acc2date)) {
        structDate->day = -1;

        if(0 == count){
            minStructDate = *structDate;
            maxStructDate = *structDate;
        }
        else{
            if( cmp_time_date(&minStructDate, structDate) > 0 ){
                minStructDate = *structDate;
            }
            if( cmp_time_date(&maxStructDate, structDate) < 0 ){
                maxStructDate = *structDate;
            }
        }
        count ++;
    }

    int minSite = 265;
    int maxSite = 29674;
    int lengthOfMask = 30001;
    bool *mask = (bool *)malloc(sizeof(bool)*lengthOfMask);
    for (int i = 0; i < lengthOfMask; ++i) {
        mask[i] = false;
    }

    bool *maskTmp = (bool *)malloc(sizeof(bool)*lengthOfMask);
    for (int i = 0; i < lengthOfMask; ++i) {
        maskTmp[i] = false;
    }
    printf("min month:\t");
    show_date(&minStructDate);
    printf("max month:\t");
    show_date(&maxStructDate);
    int *numOfMut_site = (int *)malloc(sizeof(int)*lengthOfMask);
    for (struct Date structDate1 = minStructDate; cmp_time_date(&structDate1, &maxStructDate) <= 0; NextMonth_struct(&structDate1)) {
        show_date(&structDate1);
        //init
        for(int i = 0; i < lengthOfMask; ++i){
            numOfMut_site[i] = 0;
        }
        count = 0;
        hashmap_foreach(acc, line, &(Mut.acc2line)) {
            structDate = hashmap_get(&acc2date.acc2date, acc);
            if(cmp_time_date(structDate, &structDate1) != 0){
                continue;
            }
            else{
                count ++;
            }
            char *lineCp = (char *)malloc(sizeof(char) * (strlen(line)+1));
            strcpy(lineCp, line);

            //calc freq
            if( 0 == cmp_time_date(structDate, &structDate1) ){
                chomp(lineCp);
                char *ptr1 = NULL;//for strtok_r
                char *ptr2 = NULL;//for strtok_r
                int columnIndex = 0; //column index starting from 0
                for (char *token1 = strtok_r(lineCp, "\t", &ptr1); token1 != NULL; token1 = strtok_r(NULL, "\t", &ptr1)) {
                    token1 = strip(token1);
                    if (2 == columnIndex) {//3ed column, mutations
                        for( char * token2 =strtok_r(token1,";",&ptr2) ; token2!=NULL ; token2=strtok_r(NULL,";",&ptr2) ){
                            //printf("%s\n", token2);
                            int pos;
                            sscanf(token2, "%d(%*9[^:]%*s", &pos);//example: 25452(SNP:C->T);
                            if(pos >= minSite && pos <= maxSite){
                                numOfMut_site[pos] += 1;
                                //printf("%d\n", pos);
                            }
                        }
                    }
                    columnIndex ++;
                }
            }
            free(lineCp);
            lineCp = NULL;
        }

        //keep or not
        printf("numOfSamples = %d\n", count);
        calcMask(numOfMut_site, lengthOfMask, minFreq, maskTmp, count);

        for (int i = 0; i < lengthOfMask; ++i) {
            if(mask[i] == false && maskTmp[i] == true){
                mask[i] = true;
            }
        }

    }

    FILE *fpout = fopen(str_output, "w");
    if(NULL == fpout){
        fprintf(stderr, "\nError: failed to create file: %s\n", str_output);
        exit(1);
    }
    for (int i = 0; i < lengthOfMask; ++i) {
        if(mask[i] == true){
            fprintf(fpout, "%d\n", i);
        }
    }

    //close
    fclose(fpout);
    fpout = NULL;
    //free
    hashmap_foreach(acc, line, &(Mut.acc2line)) {
        free((char*)acc);
        acc = NULL;
        free(line);
        line = NULL;
    }
    hashmap_cleanup(&(Mut.acc2line));

    hashmap_foreach(acc, structDate, &(acc2date.acc2date)) {
        //printf("%s\n%d\n%d\n%d\n\n", acc, structDate->year, structDate->month, structDate->day);
        free((char*)acc);
        acc = NULL;
        free(structDate);
        structDate = NULL;
    }
    hashmap_cleanup(&(acc2date.acc2date));
    free(maskTmp);
    maskTmp =NULL;
    free(mask);
    mask =NULL;
    free(str_mut_out);
    str_mut_out =NULL;
}

//void version() {
//    printf("McAN_siteMask v0_3\n");
//}

void help_sm() {
    printf("Usage: McAN siteMask [option]...\n\n");

    printf("Options:\n");
    printf("  -f, --vcf <file>       input vcf file\n");
    printf("  -u, --mutation <file>  input mutation file\n");
    printf("  -m, --meta <file>      input meta file\n");
    printf("  -q, --minfreq <int>    filtered out variants with frequency < minfreq\n");
    printf("  -o, --out <file>       output sitemask file\n");
    printf("  -M, --oMutation        convert vcf into mutation format and output\n\n");

    printf("  -h, --help             display this help and exit\n");
    printf("  -v, --version          output version information and exit\n");
}

enum EOption_sm {
    Evcf_sm = 0,
    Emutation_sm,
    Emeta_sm,
    Eminfreq_sm,
    Eout_sm,
    EoMutation_sm,
    Eversion_sm,
    Ehelp_sm,
    ELast_sm
};

typedef struct {
    const char* vcf;
    const char* mutation;
    const char* meta;
    const char* minfreq;
    const char* out;
    bool oMutation;
    bool version;
    bool help;
    bool used[ELast_sm];
} Argument_sm;

static const Option optionList_sm[] = {
    { "--vcf", 1, Evcf_sm },
    { "--mutation", 1, Emutation_sm },
    { "--meta", 1, Emeta_sm },
    { "--minfreq", 1, Eminfreq_sm },
    { "--out", 1, Eout_sm },
    { "--oMutation", 0, EoMutation_sm },
    { "--version", 0, Eversion_sm },
    { "--help", 0, Ehelp_sm },
    { "-f", 1, Evcf_sm },
    { "-u", 1, Emutation_sm },
    { "-m", 1, Emeta_sm },
    { "-q", 1, Eminfreq_sm },
    { "-o", 1, Eout_sm },
    { "-M", 0, EoMutation_sm },
    { "-v", 0, Eversion_sm },
    { "-h", 0, Ehelp_sm }
};

bool getOptions_sm(int argc, char* argv[], Argument_sm* arg) {
    if (argc <= 1) {
        help_sm();
        return false;
    }

    Opt opt;
    int index = 1;
    bool stopParsing = false;
    bool ret = true;
    int optionListSize = sizeof(optionList_sm) / sizeof(*optionList_sm);

    while (parseCommandLine(argc, argv, optionList_sm, optionListSize, &opt, &index, 1, 1, 0)) {
        int i = optionList_sm[opt.index].val;
        const char* name = optionList_sm[opt.index].name;
        const char* val = opt.val;

        if (arg->used[i]) {
            if (name[0] == 0)
                fprintf(stderr, "Error: repeated argument '%s'\n", val);
            else
                fprintf(stderr, "Error: repeated option '%s'\n", name);
            ret = false;
            break;
        }
        arg->used[i] = true;

        if (optionList_sm[opt.index].hasArg && (val == NULL || val[0] == 0)) {
            fprintf(stderr, "Error: option '%s' requires an argument.\n", name);
            ret = false;
            break;
        }

        char* pEnd;
        double dnum = 0;
        bool invalid = false;

        if (i == Eminfreq_sm) {
            dnum = strtod(val, &pEnd);
            if (*pEnd != 0)
                invalid = true;
        }

        switch (i)
        {
        case Evcf_sm:
            arg->vcf = val;
            break;
        case Emutation_sm:
            arg->mutation = val;
            break;
        case Emeta_sm:
            arg->meta = val;
            break;
        case Eminfreq_sm:
            arg->minfreq = val;
            break;
        case Eout_sm:
            arg->out = val;
            break;
        case EoMutation_sm:
            arg->oMutation = true;
            break;
        case Eversion_sm:
            version();
            stopParsing = true;
            break;
        case Ehelp_sm:
            help_sm();
            stopParsing = true;
            break;
        default:
            invalid = true;
        }

        if (invalid) {
            if (name[0] == 0)
                fprintf(stderr, "Error: invalid argument '%s'\n", val);
            else
                fprintf(stderr, "Error: invalid '%s' argument '%s'\n", name, val);
            ret = false;
            break;
        }

        if (stopParsing)
            break;
    }

    if (ret && !stopParsing && index < argc)
        ret = false;

    return ret;
}

void opt_siteMask(int argc, char **argv){
    char * str_inputVCF = NULL;
    char * str_inputMut = NULL;
    char * str_inputMeta = NULL;
    char * str_output = NULL;
    char * str_minFreq = NULL;
    bool mutOutput = false;
    //bool bool_help = false;
    //bool bool_version = false;
    //int c;

    Argument_sm arg;
    memset(&arg, 0, sizeof(arg));
    if (!getOptions_sm(argc, argv, &arg))
        exit(EXIT_FAILURE);
    else if (arg.used[Ehelp_sm] || arg.used[Eversion_sm])
        exit(EXIT_SUCCESS);

    if (arg.vcf) {
        str_inputVCF = (char *)malloc(sizeof(char) * (1 + strlen(arg.vcf)));
        strcpy(str_inputVCF, arg.vcf);
    }

    if (arg.mutation) {
        str_inputMut = (char *)malloc(sizeof(char) * (1 + strlen(arg.mutation)));
        strcpy(str_inputMut, arg.mutation);
    }

    if (arg.meta) {
        str_inputMeta = (char *)malloc(sizeof(char) * (1 + strlen(arg.meta)));
        strcpy(str_inputMeta, arg.meta);
    }

    if (arg.out) {
        str_output = (char *)malloc(sizeof(char) * (1 + strlen(arg.out)));
        strcpy(str_output, arg.out);
    }

    mutOutput = arg.oMutation;

    if (arg.minfreq) {
        str_minFreq = (char *)malloc(sizeof(char) * (1 + strlen(arg.minfreq)));
        strcpy(str_minFreq, arg.minfreq);
    }

    //while (1) {
    //    int this_option_optind = optind ? optind : 1;
    //    int option_index = 0;
    //    static struct option long_options[] = {
    //            {"vcf",     required_argument, 0, 0},
    //            {"mutation",     required_argument, 0, 0},
    //            {"meta",     required_argument, 0, 0},
    //            {"out",    required_argument, 0, 0},
    //            {"oMutation",    no_argument, 0, 0},
    //            {"minfreq", required_argument, 0, 0},
    //            {"help", no_argument, 0, 0},
    //            {"version", no_argument, 0, 0},
    //            {0, 0,                           0, 0}
    //    };
    //    c = getopt_long(argc, argv, "",
    //                    long_options, &option_index);
    //    if (c == -1)
    //        break;
    //    switch (c) {
    //        case 0:
    //            if (0 == strcmp(long_options[option_index].name, "vcf")) {
    //                str_inputVCF = (char *) malloc(sizeof(char) * (1 + strlen(optarg)));
    //                strcpy(str_inputVCF, optarg);
    //            }
    //            if (0 == strcmp(long_options[option_index].name, "mutation")) {
    //                str_inputMut = (char *) malloc(sizeof(char) * (1 + strlen(optarg)));
    //                strcpy(str_inputMut, optarg);
    //            }
    //            else if (0 == strcmp(long_options[option_index].name, "meta")) {
    //                str_inputMeta = (char *) malloc(sizeof(char) * (1 + strlen(optarg)));
    //                strcpy(str_inputMeta, optarg);
    //            }
    //            else if (0 == strcmp(long_options[option_index].name, "out")) {
    //                str_output = (char *) malloc(sizeof(char) * (1 + strlen(optarg)));
    //                strcpy(str_output, optarg);
    //            }
    //            else if( 0 == strcmp(long_options[option_index].name, "oMutation") ){
    //                mutOutput = true;
    //            }
    //            else if( 0 == strcmp(long_options[option_index].name, "minfreq") ){
    //                str_minFreq = (char *)malloc(sizeof(char) * (1+strlen(optarg)));
    //                strcpy(str_minFreq, optarg);
    //            }
    //            else if( 0 == strcmp(long_options[option_index].name, "help") ){
    //                bool_help = true;
    //            }
    //            else if( 0 == strcmp(long_options[option_index].name, "version") ){
    //                bool_version = true;
    //            }
    //            break;
    //        case '?':
    //            break;

    //        default:
    //            printf("?? getopt returned character code 0%o ??\n", c);
    //    }
    //}

    //if (optind < argc) {
    //    printf("non-option ARGV-elements: ");
    //    while (optind < argc)
    //        printf("%s ", argv[optind++]);
    //    printf("\n");
    //}
    //if(true == bool_help){
    //    help_sm();
    //    exit(0);
    //}
    //else if(true == bool_version){
    //    version();
    //    exit(0);
    //}
    /*else*/ if(!str_inputVCF && !str_inputMut) {
        fprintf(stderr, "No input file.\n");
        exit(-1);
    }
    else if(str_inputVCF && str_inputMut) {
        fprintf(stderr, "Duplicate input file.\n");
        exit(-1);
    }
    else if(!str_inputMeta) {
        fprintf(stderr, "No input meta file.\n");
        exit(-1);
    }
    else if(!str_output) {
        fprintf(stderr, "No output file.\n");
        exit(-1);
    }
    else{
        if (str_inputMut && mutOutput) {
            fprintf(stderr, "Take mutation format as input file. Option '--oMutation' is ignored.\n");
            mutOutput = false;
        }
        siteMask(str_inputVCF, str_inputMut, str_inputMeta, str_output, str_minFreq, mutOutput);
        //run
    }

    //free
    free(str_inputVCF);
    str_inputVCF = NULL;
    free(str_inputMut);
    str_inputMut = NULL;
    free(str_inputMeta);
    str_inputMeta = NULL;
    free(str_output);
    str_output = NULL;
    free(str_minFreq);
    str_minFreq = NULL;
}


int main_siteMask(int argc, char* argv[]) {
    opt_siteMask(argc, argv);
    printf("\nEND\n");
    return 0;
}