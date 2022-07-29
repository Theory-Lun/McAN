#include <stdio.h>
#include <time.h>
#include <string.h>

#define DEBUG

void opt_hapNet(int argc, char **argv);
int main_siteMask(int argc, char* argv[]);

int main_McAN(int argc, char* argv[]) {
#ifdef DEBUG
    clock_t start, end;
    double duration;
    start = clock();
#endif

    opt_hapNet(argc, argv);

#ifdef DEBUG
    end = clock();
    duration = -((double)(start - end)) / CLOCKS_PER_SEC;
    if(duration < 0){
        duration = 0;
    }
    printf("\n###RUNNING TIME START###\n");
    printf("running time = %f s\n", duration);
    printf("###RUNNING TIME END###\n");
    printf("\nEND\n");
#endif
    return 0;
}

int main(int argc, char* argv[]) {
    if (argc >= 2 && argv[1] && strcmp(argv[1], "siteMask") == 0)
        return main_siteMask(argc - 1, argv + 1);
    else
        return main_McAN(argc, argv);
}
