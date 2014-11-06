/* 
 * File:   input.h
 * Author: andy
 *
 * Created on February 11, 2013, 12:51 AM
 */

#ifndef INPUT_H
#define	INPUT_H

#include <string>
#include <vector>

const long long MEMORY_LIMIT = 1000000000;

class Cluster;

struct InputFile {
  std::string filename;
  bool is_directed;
  InputFile(const std::string &filenameval, bool is_directedval)
  : filename(filenameval), is_directed(is_directedval) {}
};

class Input {
public:
  double ComputeSimilarity(const Cluster& cls1, const Cluster& cls2) const;
  static Input* inst();
  
  bool IsIsolatedNode(int x1) const;
  
  inline double GetNodeCount() {
  	return n_;
  }
  
  static void SetInputFiles(const std::vector<InputFile> &input_files) { input_files_ = input_files; }
  
  static void ProcessTopCommunities(const char *TOP_COMMUNITY_FILE, const char *EDGE_FILE, const char *OUTFILE, const char *CONV_TOP_CMTY_FILE);
	//inline double GetSimilarity(int i, int j){return sim_[i][j];} // for simplicity, not check invalid index
  void test();
private:
  Input();
  Input operator=(const Input&);
  Input(const Input&);
  
  static Input* inp_ptr_;
  
  static std::vector<InputFile> input_files_;

  int n_; // node number
  
  double GetWeight(int x1, int x2) const;
  void AddWeightPair(int x1, int x2);
  
  double GetPrecomputedSimilarity(int x1, int x2) const;
  
//  double GetPearsonSimilarity(int x1, int x2, const double * const * weight, const double* sigma, const double* miu);
  double GetPearsonSimilarity(int x1, int x2) const;
  
};

#endif	/* INPUT_H */

