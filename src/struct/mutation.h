//
// Created by lun on 3/4/21.
//

#ifndef HAPNET_V0_2_1_SNP_H
#define HAPNET_V0_2_1_SNP_H

#define MUTATION_INDEL 0
#define MUTATION_DELETION 1
#define MUTATION_INSERTION 2
#define MUTATION_SNP 3

struct Mutation{
    short int pos;//position. range: 1, 2, ..., const int maxLengthOfSeq
    char *ref;//
    char *alt;//
    short int substitutionType;//1: Transition (A<->G, T<->C); 0: Transversion (other AGTC); -1: don't know (others)
    int mutationType;
};

#endif //HAPNET_V0_2_1_SNP_H
