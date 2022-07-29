//
// Created by lun on 3/18/21.
//

#ifndef HAPNET_V0_2_4_QUEUE_H
#define HAPNET_V0_2_4_QUEUE_H


// A linked list (LL) node to store a queue entry
struct QNode {
    int key;
    struct QNode* next;
};

// The queue, front stores the front node of LL and rear stores the
// last node of LL
struct Queue {
    struct QNode *front, *rear;
};

#endif //HAPNET_V0_2_4_QUEUE_H
