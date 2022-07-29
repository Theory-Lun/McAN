//
// Created by lun on 3/8/21.
//

#include <string.h>
#include "virus.h"
#include <stdbool.h>


bool emptyStr(char * str){
    //true: all empty

    if(0 == strlen(str)){
        return true;
    }

    for (int i = 0; i < strlen(str); ++i) {
        if(NULL == strchr(" \t\r\n\f\v", str[i])){//str[i] is not empty
            return false;
        }
    }

    return true;
}