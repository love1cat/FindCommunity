#include <set>
#include "input.h"
#include "cluster.h"
#include "findcluster.h"

namespace {
  const std::string CLUSTER_FILE = "clusters.txt";
}

void FindCluster::run(const int THRESHOLD) const {
  std::cout << "Starting finding community..." << std::endl;
  int n = ip_->GetNodeCount();
  std::set<ClusterPtr> clset;
  for (int i = 0; i < n; ++i) {
    ClusterPtr cp(new Cluster());
    if (!(ip_->IsIsolatedNode(i))) {
      cp->Insert(i);
      clset.insert(cp);
    }
  }
  
  // begin dividing
  std::cout << "Threadhold = " << THRESHOLD << std::endl;
  int count = 0;
  while (clset.size() > THRESHOLD) {
    std::cout << "Size of clusters = " << clset.size() << std::endl;
    std::cout << "Running new merging round " << ++count << std::endl;
    // find two cluster with best similarity (max in pearson similarity case)
    double max = -std::numeric_limits<double>::max();
    std::set<ClusterPtr>::iterator it1, it2, maxit1, maxit2;
    for (it1 = clset.begin(); it1 != clset.end(); ++it1) {
      Cluster& c1 = *((*it1).get());
      std::cout << "Computing similarity -- Cluster " << c1.id() << " vs others..." << std::endl;
      for (it2 = clset.begin(); it2 != clset.end(); ++it2) if (it1->get()->id() != it2->get()->id()) {
        Cluster& c2 = *((*it2).get());
        double s = ip_->ComputeSimilarity(c1, c2);
        if (s > max) {
          max = s;
          maxit1 = it1;
          maxit2 = it2;
        }
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
  std::cout << "Printing and saving community..." << std::endl;
  FILE *fp = fopen(CLUSTER_FILE.c_str(), "w");
  assert(fp && "Open cluser file failed.");
  std::set<ClusterPtr>::iterator it;
  for(it=clset.begin();it!=clset.end();++it){
    (*it)->Print();
    (*it)->WriteToFile(fp);
  }
  fclose(fp);
  std::cout << "Done finding community..." << std::endl;
}
