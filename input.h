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

class Input{
public:
    double ComputeSimilarity(const Cluster& cls1, const Cluster& cls2);
    static Input* inst();
    inline double GetNodeCount(){return n_;}
    //inline double GetSimilarity(int i, int j){return sim_[i][j];} // for simplicity, not check invalid index
private:
    Input();
    Input operator=(const Input&);
    Input(const Input&);
    
    static Input* inp_ptr_;
    static std::string inputfile_;
    
    int n_; // node number
    double **sim_; // similarity
    
    double GetPearsonSimilarity(int x1, int x2, const double * const * weight, const double* sigma, const double* miu);
};

#endif	/* INPUT_H */

