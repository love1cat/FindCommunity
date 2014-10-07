#include <limits>
#include <cmath>
#include <queue>
#include "input.h"
#include "cluster.h"

#define DEBUG

// if using orginal graph, comment out unify part below, see comments.
std::string Input::inputfile_ = "../nwsim/nwdataorg.txt";
Input* Input::inp_ptr_ = NULL;

Input* Input::inst() {
    if (!inp_ptr_) {
        inp_ptr_ = new Input();
    }
    return inp_ptr_;
}

double Input::GetPearsonSimilarity(int x1, int x2, const double * const * weight, const double* sigma, const double* miu) {
    double sum = 0.0;
    for (int i = 0; i < n_; ++i) sum += (weight[x1][i] - miu[x1])*(weight[x2][i] - miu[x2]);
    if (sigma[x1] == 0 || sigma[x2] == 0) return std::numeric_limits<double>::min();
    return sum / ((double) n_ * sigma[x1] * sigma[x2]);
}

Input::Input() {
    int is_weighted;
    FILE *fp = fopen(inputfile_.c_str(), "r");
    fscanf(fp, "%d", &n_);
    fscanf(fp, "%d", &is_weighted);
    int* cnt = new int[n_];
    int* tw = new int[n_]; // total weight for weighted graph, degree for unweighted graph
    double **w = new double*[n_];
    for (int i = 0; i < n_; ++i) w[i] = new double[n_];
    double* sigma = new double[n_];
    double* miu = new double[n_];
    for (int i = 0; i < n_; ++i) {
        fscanf(fp, "%d", &cnt[i]);
        fscanf(fp, "%d", &tw[i]);
        //memset(w[i], 0, sizeof (int) *n_);
        for (int j = 0; j < cnt[i]; ++j) w[i][j] = 0;
        for (int j = 0; j < cnt[i]; ++j) {
            int x;
            fscanf(fp, "%d", &x);
            ++w[i][x];
        }
    }
    fclose(fp);

    // commment this part out if using original graph
//    for (int i = 0; i < n_; ++i)
//        for (int j = 0; j < n_; ++j)
//            w[i][j] /= (double) tw[j];

    // compute miu and sigma of w[][]
    for (int i = 0; i < n_; ++i) {
        double miusum = 0.0;
        for (int j = 0; j < n_; ++j) miusum += w[i][j];
        miu[i] = miusum / (double) n_;
        double sigmasum = 0.0;
        for (int j = 0; j < n_; ++j) sigmasum += (w[i][j] - miu[i])*(w[i][j] - miu[i]);
        sigma[i] = sqrt(sigmasum / (double) n_);
    }

    sim_ = new double*[n_];
    for (int i = 0; i < n_; ++i) sim_[i] = new double[n_];

    //#ifdef DEBUG
    //    std::priority_queue<double> testq;
    //#endif

    printf("****Similarity****\n");
    for (int i = 0; i < n_; ++i)
        for (int j = 0; j < n_; ++j) {
            sim_[i][j] = GetPearsonSimilarity(i, j, w, sigma, miu);
//            printf("%d - %d:\t%.4f\n", i + 1, j + 1, sim_[i][j]);
            //#ifdef DEBUG
            //            testq.push(sim_[i][j])
            //#endif
        }
    printf("****End of Similarity****\n\n");

    //#ifdef DEBUG
    //    while(!testq.empty()){
    //        
    //    }
    //#endif

    delete[] cnt;
    for (int i = 0; i < n_; ++i) delete[] w[i];
    delete[] w;
    delete[] tw;
    delete[] sigma;
    delete[] miu;
}

double Input::ComputeSimilarity(const Cluster& cls1, const Cluster& cls2) {
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
            sum += sim_[*it1][*it2];
        }
    return sum / (double) (cls1.GetSize() * cls2.GetSize());
}
