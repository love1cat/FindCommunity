#include <cstdio>
#include "cluster.h"

int Cluster::id_count = 0;

Cluster Cluster::Merge(const Cluster& cls2){
    ids_.insert(cls2.GetIDs().begin(), cls2.GetIDs().end());
    return *this;
}

void Cluster::Print(){
    printf("****************\n");
    printf("Cluster content:\n");
    std::set<int>::iterator it;
    for(it=ids_.begin();it!=ids_.end();++it){
        printf("%d ", (*it));
    }
    printf("\n****************\n\n");
}
