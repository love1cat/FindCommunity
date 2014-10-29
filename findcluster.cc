#include <set>
#include "input.h"
#include "cluster.h"
#include "findcluster.h"

void FindCluster::run(const int THRESHOLD) const {
  std::cout << "Starting finding community..." << std::endl;
  int n = ip_->GetNodeCount();
  std::set<ClusterPtr> clset;
  for (int i = 0; i < n; ++i) {
    ClusterPtr cp(new Cluster());
    cp->Insert(i);
    clset.insert(cp);
  }
  
  // begin dividing
  std::cout << "Threadhold = " << THRESHOLD << std::endl;
  while (clset.size() > THRESHOLD) {
    std::cout << "Size of clusters = " << clset.size() << std::endl;
    std::cout << "Running new merging round..." << std::endl;
    // find two cluster with best similarity (max in pearson similarity case)
    double max = -std::numeric_limits<double>::max();
    std::set<ClusterPtr>::iterator it1, it2, maxit1, maxit2;
    for (it1 = clset.begin(); it1 != clset.end(); ++it1)
      for (it2 = clset.begin(); it2 != clset.end(); ++it2) if (it1 != it2) {
        Cluster& c1 = *((*it1).get());
        Cluster& c2 = *((*it2).get());
        double s = ip_->ComputeSimilarity(c1, c2);
        if (s > max) {
          max = s;
          maxit1 = it1;
          maxit2 = it2;
        }
      }
    
    // merge the two cluster with max similarity and delete second cluster
    Cluster& c1 = *((*maxit1).get());
    Cluster& c2 = *((*maxit2).get());
    //        printf("Two clusters with max similarity found are:\n");
    //        printf("Cluster 1 is:\n");
    //        c1.Print();
    //        printf("Cluster 2 is:\n");
    //        c2.Print();
    std::cout << "Merging Cluster " << c1.id() << " and Cluster " << c2.id() << "..." << std::endl;
    c1.Merge(c2);
    clset.erase(maxit2);
  }
  
  // print clusters
  std::cout << "Printing community..." << std::endl;
  std::set<ClusterPtr>::iterator it;
  for(it=clset.begin();it!=clset.end();++it){
    (*it)->Print();
  }
  std::cout << "Done finding community..." << std::endl;
}
