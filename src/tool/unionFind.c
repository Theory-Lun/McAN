//
// Created by lun on 3/10/21.
//

int findset(int r, int *pr) { //slow
    if(pr[r] == r) return r; // x is the root
    return findset(pr[r], pr);
}