#include <limits>
#include <cmath>
#include <queue>

#include <boost/scoped_ptr.hpp>
#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>

#include "input.h"
#include "cluster.h"

// if using orginal graph, comment out unify part below, see comments.
namespace {

const std::string INPUT_FILE = "nwdata.txt";
const char COMMENT_CHAR = '#';

bool is_comment(const char *line)
{
  // find first charactor
  while (isspace(*line)) {
    ++line;
  }
  return *line == COMMENT_CHAR;
}

const int MAXN = 1000;

// Assistant node vector 

struct Node {
  double sigma;
  double miu;
  std::set<int> nb;
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

typedef boost::unordered_map<Pair, double> Similarity_t;
Similarity_t sim_cache;
unsigned long long CACHE_CAPACITY = 500000000;
}

Input* Input::inp_ptr_ = NULL;

Input* Input::inst(bool is_walker_graph)
{
  if (!inp_ptr_) {
    inp_ptr_ = new Input(is_walker_graph);
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
  std::vector<int> nbvec;
  nbvec.insert(nbvec.end(), ns[x1].nb.begin(), ns[x1].nb.end());
  nbvec.insert(nbvec.end(), ns[x2].nb.begin(), ns[x2].nb.end());
  for (int i = 0; i < nbvec.size(); ++i) {
    int nb = nbvec[i];
    double weight1 = get_weight(x1, nb);
    double weight2 = get_weight(x2, nb);
    ret += (weight1 - miu1) * (weight2 - miu2) / divid;
  }
  
  // Non-neighbor part
  ret += miu1 * miu2 * (n_ - nbvec.size()) / divid;
  
  if (sim_cache.size() >= CACHE_CAPACITY) {
    sim_cache.erase(sim_cache.begin());
  }
  
  sim_cache.insert(Similarity_t::value_type(std::make_pair(x1, x2), ret));

  return ret;
}

Input::Input(bool is_walker_graph)
{
  FILE *fp = fopen(INPUT_FILE.c_str(), "r");
  if (!fp) {
    throw "Cannot open input file.";
  }

  char* line = new char[MAXN];
  // ID hash table for ID check
  typedef boost::unordered_set<int> ID_t;
  ID_t ids;
  unsigned int idcount = 0;
  unsigned int maxid = 0;
  unsigned int minid = std::numeric_limits<unsigned int>::max();

  while (fgets(line, MAXN, fp) != NULL) {
    if (is_comment(line)) {
      continue;
    }
    int id1, id2;
    sscanf(line, "%d%d", &id1, &id2);

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
  
  fclose(fp);
  
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

//Input::Input(bool is_walker_graph = true)
//{
//  int is_weighted;
//  FILE *fp = fopen(inputfile_.c_str(), "r");
//  fscanf(fp, "%d", &n_);
//  fscanf(fp, "%d", &is_weighted);
//  int* cnt = new int[n_];
//  int* tw = new int[n_]; // total weight for weighted graph, degree for unweighted graph
//  double **w = new double*[n_];
//  for (int i = 0; i < n_; ++i) w[i] = new double[n_];
//  double* sigma = new double[n_];
//  double* miu = new double[n_];
//  for (int i = 0; i < n_; ++i) {
//    fscanf(fp, "%d", &cnt[i]);
//    fscanf(fp, "%d", &tw[i]);
//    //memset(w[i], 0, sizeof (int) *n_);
//    for (int j = 0; j < cnt[i]; ++j) w[i][j] = 0;
//    for (int j = 0; j < cnt[i]; ++j) {
//      int x;
//      fscanf(fp, "%d", &x);
//      ++w[i][x];
//    }
//  }
//  fclose(fp);
//
//  if (is_walker_graph) {
//    for (int i = 0; i < n_; ++i)
//      for (int j = 0; j < n_; ++j)
//        w[i][j] /= (double) tw[j];
//  }
//
//  // compute miu and sigma of w[][]
//  for (int i = 0; i < n_; ++i) {
//    double miusum = 0.0;
//    for (int j = 0; j < n_; ++j) miusum += w[i][j];
//    miu[i] = miusum / (double) n_;
//    double sigmasum = 0.0;
//    for (int j = 0; j < n_; ++j) sigmasum += (w[i][j] - miu[i])*(w[i][j] - miu[i]);
//    sigma[i] = sqrt(sigmasum / (double) n_);
//  }
//
//  sim_ = new double*[n_];
//  for (int i = 0; i < n_; ++i) sim_[i] = new double[n_];
//
//  //#ifdef DEBUG
//  //    std::priority_queue<double> testq;
//  //#endif
//
//  printf("****Similarity****\n");
//  for (int i = 0; i < n_; ++i)
//    for (int j = 0; j < n_; ++j) {
//      sim_[i][j] = GetPearsonSimilarity(i, j, w, sigma, miu);
//      //            printf("%d - %d:\t%.4f\n", i + 1, j + 1, sim_[i][j]);
//      //#ifdef DEBUG
//      //            testq.push(sim_[i][j])
//      //#endif
//    }
//  printf("****End of Similarity****\n\n");
//
//  //#ifdef DEBUG
//  //    while(!testq.empty()){
//  //        
//  //    }
//  //#endif
//
//  delete[] cnt;
//  for (int i = 0; i < n_; ++i) delete[] w[i];
//  delete[] w;
//  delete[] tw;
//  delete[] sigma;
//  delete[] miu;
//}

double Input::ComputeSimilarity(const Cluster& cls1, const Cluster& cls2) const
{
  // cls1 should not be equal to cls2
  if (cls1.GetIDs().size() == 0 || cls2.GetIDs().size() == 0) throw ("cluster cannot be size 0!");
  std::set<int>::iterator it1, it2;
  //    double max = -std::numeric_limits<double>::max();
  //    for (it1 = cls1.GetIDs().begin(); it1 != cls1.GetIDs().end(); ++it1)
  //        for (it2 = cls2.GetIDs().begin(); it2 != cls2.GetIDs().end(); ++it2) {
  //            double s = sim_[*it1][*it2];
  //            if (s > max) max = s;
  //        }
  //    return max;

  double sum = 0;
  for (it1 = cls1.GetIDs().begin(); it1 != cls1.GetIDs().end(); ++it1)
    for (it2 = cls2.GetIDs().begin(); it2 != cls2.GetIDs().end(); ++it2) {
      //sum += sim_[*it1][*it2];
      sum += GetPearsonSimilarity(*it1, *it2);
    }
  return sum / (double) (cls1.GetSize() * cls2.GetSize());
}
