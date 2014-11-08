#include <set>

#include "input.h"
#include "findcluster.h"

namespace {
  
  const std::string CLUSTER_FILE = "clusters.txt";
  
  typedef std::pair<int, int> Pair;
  typedef boost::unordered_map<Pair, double> ClusterSimilarity_t;
  typedef boost::unordered_map<int, boost::unordered_set<int> > ClusterNonZeroSimPairs_t;
  
  ClusterSimilarity_t cls_sim;
  ClusterNonZeroSimPairs_t clssim_pr;
  
  void AddSimilarityPair(int id1, int id2) {
    ClusterNonZeroSimPairs_t::iterator it = clssim_pr.find(id1);
    if (it == clssim_pr.end()) {
      boost::unordered_set<int> prset1;
      prset1.insert(id2);
      clssim_pr.insert(ClusterNonZeroSimPairs_t::value_type(id1, prset1));
      return;
    }
    
    it->second.insert(id2);
  }
  
  void RemoveSimilarities(int id) {
    // Remove id related sim entry from hash table.
    ClusterNonZeroSimPairs_t::iterator it = clssim_pr.find(id);
    if (it != clssim_pr.end()) {
      boost::unordered_set<int> &ps = it->second;
      for (boost::unordered_set<int>::iterator sit = ps.begin(); sit != ps.end(); ++sit) {
        int id2 = *sit;
        cls_sim.erase(Pair(id, id2));
        cls_sim.erase(Pair(id2, id));
        
        ClusterNonZeroSimPairs_t::iterator it2 = clssim_pr.find(id2);
        it2->second.erase(id);
      }
      clssim_pr.erase(id);
    }
  }
  
  bool AddSimilarity(int id, FindCluster::ClusterMap_t& clmap, const Input* ip) {
    FindCluster::ClusterMap_t::iterator it1 = clmap.find(id);
    assert(it1 != clmap.end());
    for (FindCluster::ClusterMap_t::iterator it2 = clmap.begin(); it2 != clmap.end(); ++it2) {
      int id2 = it2->first;
      if (id == id2) continue;
      ClusterSimilarity_t::iterator sim_it = cls_sim.find(Pair(id2, id));
      if (sim_it != cls_sim.end()) {
        continue;
      }
      const Cluster &c1 = *(it1->second.get());
      const Cluster &c2 = *(it2->second.get());
      double sim = ip->ComputeSimilarity(c1, c2);
      if (sim != 0) {
        cls_sim.insert(ClusterSimilarity_t::value_type(Pair(id, id2), sim));
        if (cls_sim.size() > MEMORY_LIMIT) {
          std::cout << "Too many similarities found. The size of the hash table exceeds the predefined limit." << std::endl;
          return false;
        }
        
        AddSimilarityPair(id, id2);
        AddSimilarityPair(id2, id);
      }
    }
    return true;
  }
  
}

void FindCluster::PrintClusters(const ClusterMap_t &clmap) const {
  std::cout << "Printing and saving community..." << std::endl;
  FILE *fp = fopen(CLUSTER_FILE.c_str(), "w");
  assert(fp && "Open cluser file failed.");
  ClusterMap_t::const_iterator it;
  for(it = clmap.begin();it != clmap.end();++it){
    it->second->Print();
    it->second->WriteToFile(fp);
  }
  fclose(fp);
}

void FindCluster::run(const int THRESHOLD) const {
  std::cout << "Starting finding community..." << std::endl;
  int n = ip_->GetNodeCount();
  
  ClusterMap_t clmap;
  for (int i = 0; i < n; ++i) {
    ClusterPtr cp(new Cluster());
    cp->Insert(i);
    cp->set_id(i);
    clmap.insert(ClusterMap_t::value_type(i, cp));
  }
  
  // Initialize similarity
  std::cout << "Initializing cluster similarities... " << std::endl;
  ClusterMap_t::iterator it1, it2;
  bool memory_limit_reached = false;
  
  long long scount = 0;
  long long raw_scount = 0;
  for (int id1 = 0; id1 < n - 1; ++id1) {
    if (memory_limit_reached) break;
    for (int id2 = id1 + 1; id2 < n; ++id2) {
      ++raw_scount;
      ClusterSimilarity_t::iterator sim_it = cls_sim.find(Pair(id2, id1));
      if (sim_it != cls_sim.end()) {
        continue;
      }
      
      /************** HACK *****************/
      // Precompute all pearson similarity here, not in input.cc
      // we then just need to compute all similarities once instead of twice
      double sim = ip_->GetPearsonSimilarity(id1, id2);
      if (sim != 0) {
        ++scount;
        if (scount % 2000 == 0) {
          // Decrease output frequency.
          std::cout << "Found 2000 non-zero similarities. Count = " << scount << std::endl;
          std::cout << "Checked " << raw_scount << " pairs." << std::endl;
        }
        cls_sim.insert(ClusterSimilarity_t::value_type(Pair(id1, id2), sim));
        if (cls_sim.size() > MEMORY_LIMIT) {
          std::cout << "Too many similarities found. The size of the hash table exceeds the predefined limit." << std::endl;
          PrintClusters(clmap);
          memory_limit_reached = true;
          break;
        }
        
        AddSimilarityPair(id1, id2);
        AddSimilarityPair(id2, id1);
      }
    }
  }
  std::cout << "----------------------------------------------------------------" << std::endl;
  std::cout << "All similarities pre-computed as cluster are singltons initially." << std::endl;
  std::cout << "The similarity hash table contains " << ip_->GetSimilaritySize() << " entries." << std::endl;
  std::cout << "----------------------------------------------------------------" << std::endl;
  std::cout << "The cluster similarity hash table contains " << cls_sim.size() << " entries." << std::endl;
  
  // Begin aggregating clusters
  std::cout << "Threadhold = " << THRESHOLD << std::endl;
  std::cout << "Starting while loop..." << std::endl;
  unsigned int count = 0;
  while (clmap.size() > THRESHOLD && !cls_sim.empty()) {
    double max = 0;
    Pair max_pair;
    
    std::cout << std::endl << "***********************" << std::endl;
    std::cout << "Round " << ++count << std::endl;
    std::cout << "Find the cluster pair with maximum similarity..." << std::endl;
    for (ClusterSimilarity_t::const_iterator it = cls_sim.begin();
         it != cls_sim.end(); ++it) {
      if (it->second >= max) {
        max = it->second;
        max_pair = it->first;
      }
    }
    
    int id1 = max_pair.first;
    int id2 = max_pair.second;
    std::cout << "Merge the pair: " << id1 << ", " << id2 << std::endl;
    ClusterMap_t::iterator it = clmap.find(id1);
    assert(it != clmap.end());
    Cluster &c1 = *(it->second.get());
    it = clmap.find(id2);
    assert(it != clmap.end());
    Cluster &c2 = *(it->second.get());
    c1.Merge(c2);
    clmap.erase(id2);
    
    std::cout << "Remove these two cluster related entries from data structures..." << std::endl;
    RemoveSimilarities(id1);
    RemoveSimilarities(id2);
    
    std::cout << "Adding new similarities from new cluster..." << std::endl;
    if(!AddSimilarity(id1, clmap, ip_)) {
      std::cout << "Memory limit reached. Go print existing clusters." << std::endl;
      break;
    }
    
    std::cout << "The cluster similarity hash table is updated to " << cls_sim.size() << " entries." << std::endl;
  }
  
  // print clusters
  PrintClusters(clmap);
  std::cout << "Done finding community..." << std::endl;
}

void FindCluster::run2(const int THRESHOLD) const {
  std::cout << "Starting finding community..." << std::endl;
  int n = ip_->GetNodeCount();
  std::set<ClusterPtr> clset;
  for (int i = 0; i < n; ++i) {
    ClusterPtr cp(new Cluster());
    //if (!(ip_->IsIsolatedNode(i))) {
    cp->Insert(i);
    clset.insert(cp);
    //}
  }
  
  // begin dividing
  std::cout << "Threadhold = " << THRESHOLD << std::endl;
  int count = 0;
  int subcount = 0;
  while (clset.size() > THRESHOLD) {
    subcount = 0;
    std::cout << "Size of clusters = " << clset.size() << std::endl;
    std::cout << "Running new merging round " << ++count << std::endl;
    // find two cluster with best similarity (max in pearson similarity case)
    double max = -std::numeric_limits<double>::max();
    std::set<ClusterPtr>::iterator it1, it2, maxit1, maxit2;
    for (it1 = clset.begin(); it1 != clset.end(); ++it1) {
      Cluster& c1 = *((*it1).get());
      if (++subcount % 100 == 0) {
        std::cout << subcount << " cluster vs others completed..." << std::endl;
      }
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
