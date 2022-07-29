//---------------------------------------------------------------------------
// Copyright (C) 2021 Bo Xu <xubo123@big.ac.cn>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU Affero General Public License as
// published by the Free Software Foundation, either version 3 of the
// License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Affero General Public License for more details.
//
// You should have received a copy of the GNU Affero General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
// Mutil-threading support for calculating distance between two haplotypes
//---------------------------------------------------------------------------

#include <atomic>
#include <cstdio>
#include <thread>
#include <vector>
#include "edges.h"
#include "haplotype.h"

extern "C" {
    void calcSubset_linear_forsorted(struct Hap * pHap_source, struct Hap * pHap_target, struct Edge * edge);
    void calcMintime(struct Hap * pHap1, struct Hap * pHap2, struct Edge * edge);
    void calcDist_thread(struct Hap* haps, int numOfHaps, struct Edge* edges, int nTask, int nThread);
}

void thread_start(Hap* haps, int numOfHaps, Edge* edges, int nTask, std::atomic<int>* firstWaitingTaskID, std::atomic<int>* finishedTasks) {
    while (true) {
        int i = (*firstWaitingTaskID)++;
        if (i >= nTask) break;

        struct Edge *tmpE = NULL;
        tmpE = &(edges[i]);
        tmpE->target = i;
        tmpE->pHap2 = &(haps[tmpE->target]);
        for (int j = i + 1; j < numOfHaps; ++j) {//source, determine if j is i's ancestor
            tmpE->source = j;
            tmpE->pHap1 = &(haps[tmpE->source]);
            calcSubset_linear_forsorted(tmpE->pHap1, tmpE->pHap2, tmpE);//2021.3.28
            if (NOSUBSET == tmpE->subset) {//only keep subset
                continue;
            }
            else {
                tmpE->distance = tmpE->pHap2->numOfMut - tmpE->pHap1->numOfMut;
                calcMintime(tmpE->pHap1, tmpE->pHap2, tmpE);
                break;//find ancestor for next hap
            }
        }

        (*finishedTasks)++;
    }
}

void calcDist_thread(struct Hap* haps, int numOfHaps, struct Edge* edges, int nTask, int nThread){
        printf("  Number of worker threads: %d\n", nThread);
        printf("  Total tasks: %d\n", nTask);
        fflush(stdout);

        std::atomic<int> firstWaitingTaskID(0);
        std::atomic<int> finishedTasks(0);
        std::vector<std::thread> threads;

        for (int i = 0; i < nThread; i++)
            threads.emplace_back(thread_start, haps, numOfHaps, edges, nTask, &firstWaitingTaskID, &finishedTasks);

        while (true) {
            printf("  Percent completed: %5.1f%%\r", finishedTasks * 100 / (double)nTask);
            fflush(stdout);

            if (finishedTasks == nTask) {
                printf("\n");
                fflush(stdout);
                break;
            }

            std::this_thread::sleep_for(std::chrono::milliseconds(100));
        }

        for (auto& elem : threads)
            elem.join();
}
