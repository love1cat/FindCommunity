/* 
 * File:   input.h
 * Author: andy
 *
 * Created on February 11, 2013, 12:51 AM
 */

#ifndef INPUT_H
#define	INPUT_H

#include <string>

class Cluster;

class Input {
public:
  double ComputeSimilarity(const Cluster& cls1, const Cluster& cls2) const;
  static Input* inst(const char * INPUT_FILE = NULL);
  
  bool IsIsolatedNode(int x1) const;
  
  inline double GetNodeCount() {
  	return n_;
  }
  
  static void ProcessTopCommunities(const char *TOP_COMMUNITY_FILE, const char *EDGE_FILE, const char *OUTFILE, const char *CONV_TOP_CMTY_FILE);
	//inline double GetSimilarity(int i, int j){return sim_[i][j];} // for simplicity, not check invalid index
private:
  Input(const char * INPUT_FILE);
  Input operator=(const Input&);
  Input(const Input&);
  
  static Input* inp_ptr_;

  int n_; // node number
  
//  double GetPearsonSimilarity(int x1, int x2, const double * const * weight, const double* sigma, const double* miu);
  double GetPearsonSimilarity(int x1, int x2) const;
};

#endif	/* INPUT_H */

