//
// Created by lun on 1/28/21.
//

#include <ctype.h>
#include <stdbool.h>
#include <string.h>

bool chomp( char *str ) {
    size_t len = strlen(str);

    /* Empty string */
    if( len == 0 ) {
        return false;
    }

    size_t last_idx = len - 1;
    if( str[last_idx] == '\n' ) {
        str[last_idx] = '\0';
        return true;
    }
    else {
        return false;
    }
}

char* strip(char* str) {
    if (str) {
        char* ptr = str + strspn(str, " \t\r\n\f\v");
        size_t size = strlen(ptr);
        while (size-- > 0)
            if (!isspace(ptr[size]))
                break;
        ptr[size + 1] = 0;
        memmove(str, ptr, size + 1);
    }
    return str;
}
