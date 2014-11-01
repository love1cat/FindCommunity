/*
 * File:   cluster.h
 * Author: andy
 *
 * Created on February 11, 2013, 12:40 AM
 */

#ifndef CLUSTER_H
#define	CLUSTER_H

#include <set>
#include <boost/unordered_set.hpp>
#include <boost/shared_ptr.hpp>

class Cluster;
typedef boost::shared_ptr<Cluster> ClusterPtr;

class Cluster{
public:
  Cluster() {id_ = ++id_count;}
  Cluster Merge(const Cluster& cls2);
  void Print() const;
  void WriteToFile(FILE *fp) const;
  inline void Insert(int x) {ids_.insert(x);}
  inline int GetSize() const {return ids_.size();}
  inline const boost::unordered_set<int>& GetIDs() const{return ids_;}
  inline int id() const {return id_;}
  inline void set_id(int id) {id_ = id;}
private:
  int id_;
  static int id_count;
  boost::unordered_set<int> ids_;
};

#endif	/* CLUSTER_H */

