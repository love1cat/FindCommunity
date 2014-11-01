/*
 * File:   findcluster.h
 * Author: andy
 *
 * Created on February 11, 2013, 6:48 PM
 */

#ifndef FINDCLUSTER_H
#define	FINDCLUSTER_H

#include <boost/unordered_map.hpp>

#include "cluster.h"

class Input;
class FindCluster{
public:
  typedef boost::unordered_map<int, ClusterPtr> ClusterMap_t;
  inline FindCluster(Input* ip) : ip_(ip) {}
  void run(const int threshold) const;
  void run2(const int threshold) const;
  void PrintClusters(const ClusterMap_t& clmap) const;
private:
  Input* ip_;
};

#endif	/* FINDCLUSTER_H */

