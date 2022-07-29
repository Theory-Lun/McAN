//
// Created by lun on 3/5/21.
//

#include <string.h>

int find_string(char * str, char * substr){
    //char * str: string
    //char * substr: substring
    //return number of substring in string

    int count = 0,i,j,check;
    int len = (int)strlen(str);
    int sublen = (int)strlen(substr);
    for(i = 0; i < len; i++){
        check = 1;
        for(j = 0; j + i < len && j < sublen; j++){
            if(str[i+j] != substr[j]){
                check = 0;
                break;
            }
        }
        if(check == 1){
            count++;
            i = i + sublen;
        }
    }
    return count;
}