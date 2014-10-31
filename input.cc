#include <limits>
#include <cmath>
#include <queue>
#include <iostream>
#include <fstream>
#include <sstream>

#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>

#include "input.h"
#include "cluster.h"

// if using orginal graph, comment out unify part below, see comments.
namespace {
  
  const char COMMENT_CHAR = '#';
  
  bool is_comment(const char *line)
  {
    // find first charactor
    while (isspace(*line)) {
      ++line;
    }
    return *line == COMMENT_CHAR;
  }
  
  // Assistant node vector
  
  struct Node {
    double sigma;
    double miu;
    boost::unordered_set<int> nb;
  };
  
  std::vector<Node> ns;
  
  // Weight hash table instead of array to save memory
  typedef std::pair<int, int> Pair;
  typedef boost::unordered_map<Pair, int> Weight_t;
  Weight_t w;
  
  double get_weight(int x1, int x2) {
    double weight = 0.0;
    Weight_t::const_iterator cit = w.find(std::make_pair(x1, x2));
    if (cit != w.end()) {
      weight = cit->second;
    }
    return weight;
  }
  
  // Caches
  unsigned long long SIM_CACHE_CAPACITY = 500000000;
  unsigned long long CLUSTER_SIM_CACHE_CAPACITY = 300000000;
  typedef boost::unordered_map<Pair, double> Similarity_t;
  Similarity_t sim_cache;
  
  typedef std::pair<const Cluster*, int> ClusterInfo_t;
  typedef std::pair<ClusterInfo_t, ClusterInfo_t> ClusterPair_t;
  typedef boost::unordered_map<ClusterPair_t, double> ClusterSimilarity_t;
  ClusterSimilarity_t cluster_sim_cache;
  
}

Input* Input::inp_ptr_ = NULL;

Input* Input::inst(const char * INPUT_FILE)
{
  if (!inp_ptr_) {
    inp_ptr_ = new Input(INPUT_FILE);
  }
  return inp_ptr_;
}

bool Input::IsIsolatedNode(int x1) const {
  return ns[x1].nb.empty();
}

double Input::GetPearsonSimilarity(int x1, int x2) const
{
  Similarity_t::const_iterator it = sim_cache.find(std::make_pair(x1, x2));
  if (it != sim_cache.end()) {
    return it->second;
  }
  
  double sigma1 = ns[x1].sigma;
  double sigma2 = ns[x2].sigma;
  if (sigma1 == 0 || sigma2 == 0) return std::numeric_limits<double>::min();
  
  double divid = sigma1 * sigma2 * n_;
  
  double ret = 0.0;
  double miu1 = ns[x1].miu;
  double miu2 = ns[x2].miu;
  
  // Neighbor part
  boost::unordered_set<int> nbset;
  nbset.insert(ns[x1].nb.begin(), ns[x1].nb.end());
  nbset.insert(ns[x2].nb.begin(), ns[x2].nb.end());
  if (nbset.size() == ns[x1].nb.size() + ns[x2].nb.size()) {
    // Not sharing any neighbor, return minimum value 0.0
    return 0.0;
  }
  
  for (boost::unordered_set<int>::const_iterator it = nbset.begin(); it != nbset.end(); ++it) {
    int nb = *it;
    double weight1 = get_weight(x1, nb);
    double weight2 = get_weight(x2, nb);
    ret += (weight1 - miu1) * (weight2 - miu2) / divid;
  }
  
  // Non-neighbor part
  ret += miu1 * miu2 * (n_ - nbset.size()) / divid;
  
  if (sim_cache.size() >= SIM_CACHE_CAPACITY) {
    sim_cache.erase(sim_cache.begin());
  }
  
  sim_cache.insert(Similarity_t::value_type(std::make_pair(x1, x2), ret));
  
  return abs(ret);
}

Input::Input(const char * INPUT_FILE)
{
  std::ifstream infile(INPUT_FILE);
  if (!infile) {
    throw "Cannot open input file.";
  }
  
  // ID hash table for ID check
  typedef boost::unordered_set<int> ID_t;
  ID_t ids;
  unsigned int idcount = 0;
  unsigned int maxid = 0;
  unsigned int minid = std::numeric_limits<unsigned int>::max();
  
  std::string line;
  while (std::getline(infile, line)) {
    if (is_comment(line.c_str())) {
      continue;
    }
    std::istringstream ssline(line);
    int id1, id2;
    ssline >> id1 >> id2;
    
    assert(id1 >= 0 && id2 >= 0);
    
    if (id1 == id2) {
      continue;
    }
    
    std::pair < Weight_t::iterator, bool> pr =
    w.insert(Weight_t::value_type(std::make_pair(id1, id2), 1));
    if (!(pr.second)) {
      // Already in the hash table
      // Increase the count
      Weight_t::iterator &wit = pr.first;
      ++(wit->second);
    }
    
    // ID check
    std::pair < ID_t::iterator, bool> pr2 =
    ids.insert(ID_t::value_type(id1));
    if (pr2.second) {
      ++idcount;
    }
    
    pr2 = ids.insert(ID_t::value_type(id2));
    if (pr2.second) {
      ++idcount;
    }
    
    if (id1 > maxid) {
      maxid = id1;
    }
    
    if (id1 < minid) {
      minid = id1;
    }
    
    if (id2 > maxid) {
      maxid = id2;
    }
    
    if (id2 < minid) {
      minid = id2;
    }
  }
  
  std::cout << "Done reading file..." <<std::endl;
  
  // ID validity check -- IDs must be 0-indexed
  std::cout << "Max ID: " << maxid << ", Min ID: " << minid << ", idcount: "<< idcount << std::endl;
  //assert(maxid == idcount - 1 && idcount > 0);
  assert(minid == 0 && idcount > 0);
  
  n_ = maxid + 1;
  // Go over hash table and obtain node neighbor,
  // miu and sigma
  
  
  // Init node vector
  ns.resize(n_);
  for (int i = 0; i < n_; ++i) {
    ns[i].miu = 0.0;
    ns[i].sigma = 0.0;
    ns[i].nb.clear();
  }
  
  // Compute mius
  std::cout << "Computing MIUs..." <<std::endl;
  for(Weight_t::const_iterator it = w.begin(); it != w.end(); ++it) {
    const Pair &p = it->first;
    int weight = it->second;
    int id1 = p.first;
    int id2 = p.second;
    // IDs must be 0-indexed
    ns[id1].miu += (double)weight / (double)n_;
    ns[id1].nb.insert(id2);
  }
  
  // Compute sigmas
  // First step, compute partial value from all neighbors
  std::cout << "Computing Sigmas, first part..." <<std::endl;
  for(Weight_t::const_iterator it = w.begin(); it != w.end(); ++it) {
    const Pair &p = it->first;
    int weight = it->second;
    int id1 = p.first;
    // IDs must be 0-indexed
    double miu = ns[id1].miu;
    ns[id1].sigma += (double)(weight - miu) / (double)n_ * (double)(weight - miu);
  }
  
  // Second step, compute partial value from all non-neighbors
  std::cout << "Computing Sigmas, second part..." <<std::endl;
  for (int i = 0; i < n_; ++i) {
    //    std::cout << "i = " << i <<std::endl;
    double &miu = ns[i].miu;
    if (miu == 0.0) continue;
    ns[i].sigma += miu / (double)n_ * miu * (n_ - ns[i].nb.size());
  }
  
  std::cout << "Done processing input..." <<std::endl;
}



double Input::ComputeSimilarity(const Cluster& cls1, const Cluster& cls2) const
{
  ClusterPair_t cls_pr(ClusterInfo_t(&cls1, cls1.GetSize()), ClusterInfo_t(&cls2, cls2.GetSize()));
  ClusterSimilarity_t::const_iterator it = cluster_sim_cache.find(cls_pr);
  if (it != cluster_sim_cache.end()) {
    return it->second;
  }
  
  // cls1 should not be equal to cls2
  if (cls1.GetIDs().size() == 0 || cls2.GetIDs().size() == 0) throw ("cluster cannot be size 0!");
  boost::unordered_set<int>::iterator it1, it2;
  
  double sum = 0;
  if (cls1.GetIDs().size() == 1 && cls2.GetIDs().size() == 1) {
    // If both are singleton clusters, only compute similarity if they are nbs.
    const int id1 = *(cls1.GetIDs().begin());
    const int id2 = *(cls2.GetIDs().begin());
    boost::unordered_set<int>::iterator it = ns[id1].nb.find(id1);
    if (it == ns[id1].nb.end()) {
      return 0.0;
    }
    
    return GetPearsonSimilarity(*it1, *it2);
  }
  
  for (it1 = cls1.GetIDs().begin(); it1 != cls1.GetIDs().end(); ++it1)
    for (it2 = cls2.GetIDs().begin(); it2 != cls2.GetIDs().end(); ++it2) {
      //sum += sim_[*it1][*it2];
      sum += GetPearsonSimilarity(*it1, *it2);
    }
  
  double ret = sum / (double) (cls1.GetSize() * cls2.GetSize());
  
  if (cluster_sim_cache.size() >= CLUSTER_SIM_CACHE_CAPACITY) {
    cluster_sim_cache.erase(cluster_sim_cache.begin());
  }
  cluster_sim_cache.insert(ClusterSimilarity_t::value_type(cls_pr, ret));
  
  return ret;
}

void Input::ProcessTopCommunities(const char *TOP_COMMUNITY_FILE, const char *EDGE_FILE, const char *OUTFILE, const char *CONV_TOP_CMTY_FILE) {
  std::ifstream infile(TOP_COMMUNITY_FILE);
  if (!infile) {
    throw "Cannot open input file.";
  }
  
  // ID hash table for ID check
  typedef boost::unordered_map<int, int> IDMap_t;
  IDMap_t idmap;
  unsigned int idcount = 0;
  
  std::string line;
  while (std::getline(infile, line)) {
    if (is_comment(line.c_str())) {
      continue;
    }
    std::istringstream ssline(line);
    int id;
    while (ssline >> id) {
      std::pair<IDMap_t::iterator, bool> pr = idmap.insert(IDMap_t::value_type(id, idcount));
      if (pr.second) {
        ++idcount;
      }
    }
  }
  
  std::cout << "ID count = " << idcount << std::endl;
  
  infile.close();
  
  // Read edge file. Remove edges not in top communities.
  infile.open(EDGE_FILE);
  std::ofstream outfile(OUTFILE, std::ios_base::trunc);
  while (std::getline(infile, line)) {
    if (is_comment(line.c_str())) {
      continue;
    }
    std::istringstream ssline(line);
    int id1, id2;
    ssline >> id1 >> id2;
    IDMap_t::const_iterator cit1 = idmap.find(id1);
    IDMap_t::const_iterator cit2 = idmap.find(id2);
    if (cit1 != idmap.end() && cit2 != idmap.end()) {
      outfile << cit1->second << "\t" << cit2->second << "\n";
    }
  }
  
  infile.close();
  outfile.close();
  
  // Convert top communities file to 0-indexed ids
  infile.open(TOP_COMMUNITY_FILE);
  outfile.open(CONV_TOP_CMTY_FILE, std::ios_base::trunc);
  while (std::getline(infile, line)) {
    if (is_comment(line.c_str())) {
      continue;
    }
    std::istringstream ssline(line);
    int id;
    while(ssline >> id) {
      IDMap_t::const_iterator it = idmap.find(id);
      if (it == idmap.end()) {
        throw "Error: id in community is not in idmap.";
      }
      outfile << it->second << " ";
    }
    outfile << "\n";
  }
}
