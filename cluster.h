/* 
 * File:   cluster.h
 * Author: andy
 *
 * Created on February 11, 2013, 12:40 AM
 */

#ifndef CLUSTER_H
#define	CLUSTER_H

#include <set>
#include <boost/shared_ptr.hpp>

class Cluster;
typedef boost::shared_ptr<Cluster> ClusterPtr;

class Cluster{
public:
    Cluster Merge(const Cluster& cls2);
    void Print();
    inline void Insert(int x) {ids_.insert(x);}
    inline int GetSize() const{return ids_.size();}
    inline const std::set<int>& GetIDs() const{return ids_;}
private:
    std::set<int> ids_;   
};

#endif	/* CLUSTER_H */

