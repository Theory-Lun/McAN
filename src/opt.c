//
// Created by lun on 3/3/21.
//
//#define DEBUG

#include <stdio.h>
#include <stdlib.h>
//#include <getopt.h>
#include <string.h>
#include <stdbool.h>
#include "optionparser.h"

//#ifdef MUTILTHREADING
//void hapNetDirected_inputMask_sysMem( const char * str_mut, const char * str_meta, const char * str_minFreq, const char * str_maxvirus, const char * str_folder, char *str_sitemask, int nThread);
//#else
//void hapNetDirected_inputMask_sysMem( const char * str_mut, const char * str_meta, const char * str_minFreq, const char * str_maxvirus, const char * str_folder, char *str_sitemask);
//#endif
#ifdef MUTILTHREADING
void hapNetDirected(const char * str_vcf, const char * str_mut, const char * str_meta, const char * str_minFreq, const char * str_maxvirus, const char * str_folder, char *str_sitemask, int* outFrmt, int nThread);
#else
void hapNetDirected(const char * str_vcf, const char * str_mut, const char * str_meta, const char * str_minFreq, const char * str_maxvirus, const char * str_folder, char *str_sitemask, int* outFrmt);
#endif

void opt__printStr(char * arg, char * opt){
    printf("%s = ", opt);
    if( NULL ==  arg){
        printf("NULL\n");
    }
    else{
        printf("%s\n", arg);
    }
}

void printAllArg( char * str_mut, char * str_sta_anno, char * str_meta, char * str_minFreq, char * str_maxvirus, char * str_folder, bool bool_help ){
    //head
    printf("\n###ALL ARGUMENTS START###\n");
    //contents
    opt__printStr(str_mut, "mutation");
    opt__printStr(str_sta_anno, "staanno");
    opt__printStr(str_meta, "meta");
    opt__printStr(str_minFreq, "minfreq");
    opt__printStr(str_maxvirus, "maxvirus");
    opt__printStr(str_folder, "folder");
    printf("%s = %d\n", "help", bool_help);
    //end
    printf("###ALL ARGUMENTS ###\n");
}

void version() {
    printf("McAN v1.2\n");
}

void help(){
    printf("Usage: McAN [option]...\n");
    printf("       McAN siteMask [option]...\n\n");

    printf("Options:\n");
    printf("  -f, --vcf <file>        input vcf file\n");
    printf("  -u, --mutation <file>   input mutation file\n");
    printf("  -m, --meta <file>       input meta file\n");
    printf("  -s, --sitemask <file>   input sitemask file [optional]\n");
    printf("  -x, --maxnsample <int>  maximum number of samples [default: Inf]\n");
#ifdef MUTILTHREADING
    printf("  -t, --nthread <int>     number of worker threads [default: 1]\n");
#endif
    printf("  -o, --outDir <dir>      directory for output\n");
    printf("  -J, --oJSON             output network with JSON format\n");
    printf("  -T, --oTSV              output network with TSV format\n");
    printf("  -G, --oGraphML          output network with GraphML format\n");
    printf("  -M, --oMutation         convert vcf into mutation format and output\n\n");

    printf("  -h, --help              display this help and exit\n");
    printf("  -v, --version           output version information and exit\n");
}

enum EOption {
    Evcf = 0,
    Emutation,
    Emeta,
    Esitemask,
    Emaxnsample,
#ifdef MUTILTHREADING
    Enthread,
#endif
    EoutDir,
    EoJSON,
    EoTSV,
    EoGraphML,
    EoMutation,
    Eversion,
    Ehelp,
    ELast
};

typedef struct {
    const char* vcf;
    const char* mutation;
    const char* meta;
    const char* sitemask;
    const char* maxnsample;
#ifdef MUTILTHREADING
    int nthread;
#endif
    const char* outDir;
    int outFrmt[4];
    bool version;
    bool help;
    bool used[ELast];
} Argument;

static const Option optionList[] = {
    { "--vcf", 1, Evcf },
    { "--mutation", 1, Emutation },
    { "--meta", 1, Emeta },
    { "--sitemask", 1, Esitemask },
    { "--maxnsample", 1, Emaxnsample },
#ifdef MUTILTHREADING
    { "--nthread", 1, Enthread },
#endif
    { "--outDir", 1, EoutDir },
    { "--oJSON", 0, EoJSON },
    { "--oTSV", 0, EoTSV },
    { "--oGraphML", 0, EoGraphML },
    { "--oMutation", 0, EoMutation },
    { "--version", 0, Eversion },
    { "--help", 0, Ehelp },
    { "-f", 1, Evcf },
    { "-u", 1, Emutation },
    { "-m", 1, Emeta },
    { "-s", 1, Esitemask },
    { "-x", 1, Emaxnsample },
#ifdef MUTILTHREADING
    { "-t", 1, Enthread },
#endif
    { "-o", 1, EoutDir },
    { "-J", 0, EoJSON },
    { "-T", 0, EoTSV },
    { "-G", 0, EoGraphML },
    { "-M", 0, EoMutation },
    { "-v", 0, Eversion },
    { "-h", 0, Ehelp }
};

bool getOptions(int argc, char* argv[], Argument* arg) {
    if (argc <= 1) {
        help();
        return false;
    }

    Opt opt;
    int index = 1;
    bool stopParsing = false;
    bool ret = true;
    int optionListSize = sizeof(optionList) / sizeof(*optionList);

    while (parseCommandLine(argc, argv, optionList, optionListSize, &opt, &index, 1, 1, 0)) {
        int i = optionList[opt.index].val;
        const char* name = optionList[opt.index].name;
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

        if (optionList[opt.index].hasArg && (val == NULL || val[0] == 0)) {
            fprintf(stderr, "Error: option '%s' requires an argument.\n", name);
            ret = false;
            break;
        }

        char* pEnd;
        int inum = 0;
        bool invalid = false;

#ifdef MUTILTHREADING
        if (i == Emaxnsample || i == Enthread) {
#else
        if (i == Emaxnsample) {
#endif
            inum = (int)strtol(val, &pEnd, 10);
            if (*pEnd != 0)
                invalid = true;
        }

        switch (i)
        {
        case Evcf:
            arg->vcf = val;
            break;
        case Emutation:
            arg->mutation = val;
            break;
        case Emeta:
            arg->meta = val;
            break;
        case Esitemask:
            arg->sitemask = val;
            break;
        case Emaxnsample:
            arg->maxnsample = val;
            break;
#ifdef MUTILTHREADING
        case Enthread:
            arg->nthread = inum;
            break;
#endif
        case EoutDir:
            arg->outDir = val;
            break;
        case EoJSON:
            arg->outFrmt[0] = 1;
            break;
        case EoTSV:
            arg->outFrmt[1] = 1;
            break;
        case EoGraphML:
            arg->outFrmt[2] = 1;
            break;
        case EoMutation:
            arg->outFrmt[3] = 1;
            break;
        case Eversion:
            version();
            stopParsing = true;
            break;
        case Ehelp:
            help();
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

void opt_hapNet(int argc, char **argv){
    //all argv
    char * str_vcf = NULL;
    char * str_mut = NULL;
    //char * str_sta_anno = NULL;
    char * str_meta = NULL;
    char * str_minFreq = NULL;
    char * str_maxvirus = NULL;
    char * str_folder = NULL;
    char * str_sitemask = NULL;
#ifdef MUTILTHREADING
    int nThread = 1;
#endif
    int outFrmt[4];
    outFrmt[0] = outFrmt[1] = outFrmt[2] = outFrmt[3] = 0;
    //bool bool_help = false;
    //bool bool_version = false;

    Argument arg;
    memset(&arg, 0, sizeof(arg));
#ifdef MUTILTHREADING
    arg.nthread = 1;
#endif
    if (!getOptions(argc, argv, &arg))
        exit(EXIT_FAILURE);
    else if (arg.used[Ehelp] || arg.used[Eversion])
        exit(EXIT_SUCCESS);

    if (arg.vcf) {
        str_vcf = (char *)malloc(sizeof(char) * (1 + strlen(arg.vcf))); // the last char store '\0'
        strcpy(str_vcf, arg.vcf);
    }

    if (arg.mutation) {
        str_mut = (char *)malloc(sizeof(char) * (1 + strlen(arg.mutation))); // the last char store '\0'
        strcpy(str_mut, arg.mutation);
    }

    if (arg.meta) {
        str_meta = (char *)malloc(sizeof(char) * (1 + strlen(arg.meta)));
        strcpy(str_meta, arg.meta);
    }

    if (arg.maxnsample) {
        str_maxvirus = (char *)malloc(sizeof(char) * (1 + strlen(arg.maxnsample)));
        strcpy(str_maxvirus, arg.maxnsample);
    }

    if (arg.outDir) {
        str_folder = (char *)malloc(sizeof(char) * (1 + strlen(arg.outDir)));
        strcpy(str_folder, arg.outDir);
    }

    if (arg.sitemask) {
        str_sitemask = (char *)malloc(sizeof(char) * (1 + strlen(arg.sitemask)));
        strcpy(str_sitemask, arg.sitemask);
    }

#ifdef MUTILTHREADING
    nThread = arg.nthread;
#endif

    memcpy(outFrmt, arg.outFrmt, sizeof(outFrmt));

    //read all arg
//    int c;
//    int digit_optind = 0;
//    while (1) {
//        int this_option_optind = optind ? optind : 1;
//        int option_index = 0;
//        static struct option long_options[] = {
//                {"vcf", required_argument, 0, 0},
//                {"mutation", required_argument, 0, 0},
//                //{"staanno", required_argument, 0, 0},
//                {"meta", required_argument, 0, 0},
//                {"minfreq", required_argument, 0, 0},
//                {"maxnsample", required_argument, 0, 0},
//                {"outDir", required_argument, 0, 0},
//                {"sitemask", required_argument, 0, 0},
//#ifdef MUTILTHREADING
//                {"nthread", required_argument, 0, 0},
//#endif
//                {"oJSON", no_argument, 0, 0},
//                {"oTSV", no_argument, 0, 0},
//                {"oGraphML", no_argument, 0, 0},
//                {"oMutation", no_argument, 0, 0},
//                {"help", no_argument, 0, 0},
//                {"version", no_argument, 0, 0},
//                {0,         0,                 0,  0 }
//        };
//        c = getopt_long(argc, argv, "",
//                        long_options, &option_index);
//        if (c == -1)
//            break;
//
//        switch (c) {
//            case 0:
//                if(0 == strcmp(long_options[option_index].name, "vcf")){
//                    str_vcf = (char *)malloc(sizeof(char) * (1+strlen(optarg))); // the last char store '\0'
//                    strcpy(str_vcf, optarg);
//                }
//                else if(0 == strcmp(long_options[option_index].name, "mutation")){
//                    str_mut = (char *)malloc(sizeof(char) * (1+strlen(optarg))); // the last char store '\0'
//                    strcpy(str_mut, optarg);
//                }
//                //else if( 0 == strcmp(long_options[option_index].name, "staanno") ){
//                //    str_sta_anno = (char *)malloc(sizeof(char) * (1+strlen(optarg)));
//                //    strcpy(str_sta_anno, optarg);
//                //}
//                else if(0 == strcmp(long_options[option_index].name, "meta")){
//                    str_meta = (char *)malloc(sizeof(char) * (1+strlen(optarg)));
//                    strcpy(str_meta, optarg);
//                }
//                else if( 0 == strcmp(long_options[option_index].name, "minfreq") ){
//                    str_minFreq = (char *)malloc(sizeof(char) * (1+strlen(optarg)));
//                    strcpy(str_minFreq, optarg);
//                }
//                else if( 0 == strcmp(long_options[option_index].name, "maxnsample") ){
//                    str_maxvirus = (char *)malloc(sizeof(char) * (1+strlen(optarg)));
//                    strcpy(str_maxvirus, optarg);
//                }
//                else if (0 == strcmp(long_options[option_index].name, "outDir")){
//                    str_folder = (char *)malloc(sizeof(char) * (1+strlen(optarg)));
//                    strcpy(str_folder, optarg);
//                }
//                else if (0 == strcmp(long_options[option_index].name, "sitemask")){
//                    str_sitemask = (char *)malloc(sizeof(char) * (1+strlen(optarg)));
//                    strcpy(str_sitemask, optarg);
//                }
//#ifdef MUTILTHREADING
//                else if (0 == strcmp(long_options[option_index].name, "nthread")){
//                    nThread = atoi(optarg);
//                    if (nThread < 1) {
//                        fprintf(stderr, "Number of threads should not be less than 1.\n");
//                        exit(-1);
//                    }
//                }
//#endif
//                else if (0 == strcmp(long_options[option_index].name, "oJSON")){
//                    outFrmt[0] = 1;
//                }
//                else if (0 == strcmp(long_options[option_index].name, "oTSV")){
//                    outFrmt[1] = 1;
//                }
//                 else if (0 == strcmp(long_options[option_index].name, "oGraphML")){
//                    outFrmt[2] = 1;
//                }
//                 else if (0 == strcmp(long_options[option_index].name, "oMutation")){
//                    outFrmt[3] = 1;
//                }
//                else if (0 == strcmp(long_options[option_index].name, "help")){
//                    bool_help = true;
//                }
//                else if (0 == strcmp(long_options[option_index].name, "version")){
//                    bool_version = true;
//                }
//                break;
//            case '?':
//                break;
//
//            default:
//                printf("?? getopt returned character code 0%o ??\n", c);
//        }
//    }
//    if (optind < argc) {
//        printf("non-option ARGV-elements: ");
//        while (optind < argc)
//            printf("%s ", argv[optind++]);
//        printf("\n");
//    }
    //run
    //if(true == bool_help){
    //    help();
    //    exit(0);
    //}
    //else if(true == bool_version){
    //    version();
    //    exit(0);
    //}
    /*else*/ if(!str_vcf && !str_mut) {
        fprintf(stderr, "No input file.\n");
        exit(-1);
    }
    else if(str_vcf && str_mut) {
        fprintf(stderr, "Duplicate input file.\n");
        exit(-1);
    }
    else if(!str_meta) {
        fprintf(stderr, "No input meta file.\n");
        exit(-1);
    }
    else if(!str_folder) {
        fprintf(stderr, "No output directory.\n");
        exit(-1);
    }
    else if(outFrmt[0] + outFrmt[1] + outFrmt[2] + outFrmt[3] == 0) {
        fprintf(stderr, "At least one output file format is required.\n");
        exit(-1);
    }
#ifdef MUTILTHREADING
    else if (nThread < 1) {
        fprintf(stderr, "Number of threads should not be less than 1.\n");
        exit(-1);
    }
#endif
    else{//run all
        if (str_mut && outFrmt[3]) {
            fprintf(stderr, "Take mutation format as input file. Option '--oMutation' is ignored.\n");
            outFrmt[3] = 0;
        }
#ifdef DEBUG
        printAllArg( str_mut, str_sta_anno, str_meta, str_minFreq, str_maxvirus, str_folder, bool_help );
#endif
        //hapNetAll( str_mut, str_sta_anno, str_meta, str_minFreq, str_maxvirus, str_folder );
        //hapNetDirected_inputMask( str_mut, str_meta, str_minFreq, str_maxvirus, str_folder, str_sitemask );
//#ifdef MUTILTHREADING
//        hapNetDirected_inputMask_sysMem(str_mut, str_meta, str_minFreq, str_maxvirus, str_folder, str_sitemask, nThread);
//#else
//        hapNetDirected_inputMask_sysMem(str_mut, str_meta, str_minFreq, str_maxvirus, str_folder, str_sitemask);
//#endif
#ifdef MUTILTHREADING
        hapNetDirected(str_vcf, str_mut, str_meta, str_minFreq, str_maxvirus, str_folder, str_sitemask, outFrmt, nThread);
#else
        hapNetDirected(str_vcf, str_mut, str_meta, str_minFreq, str_maxvirus, str_folder, str_sitemask, outFrmt);
#endif
    }
    //free
    free(str_vcf);
    str_mut = NULL;
    free(str_mut);
    str_mut = NULL;
    free(str_folder);
    str_folder = NULL;
    //free(str_sta_anno);
    //str_sta_anno = NULL;
    free(str_minFreq);
    str_minFreq = NULL;
    free(str_maxvirus);
    str_maxvirus = NULL;
    free(str_meta);
    str_meta = NULL;
    free(str_sitemask);
    str_sitemask = NULL;
}
