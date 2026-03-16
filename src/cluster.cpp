// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::plugins(cpp17)]]

#include <Rcpp.h>
#include <RcppParallel.h>

using namespace Rcpp;
using namespace RcppParallel;

struct SparseGraphBuilder : public Worker {
  
  const RMatrix<double> X;
  const RVector<double> rt;
  double rt_window;
  double alpha;
  int n_samples;
  int n_features;
  
  // Thread-safe containers
  std::vector< std::vector<int> > from_list;
  std::vector< std::vector<int> > to_list;
  std::vector< std::vector<double> > weight_list;
  
  SparseGraphBuilder(NumericMatrix X,
                     NumericVector rt,
                     double rt_window,
                     double alpha)
    : X(X), rt(rt), rt_window(rt_window), alpha(alpha),
      n_samples(X.ncol()), n_features(X.nrow())
  {
    from_list.resize(n_features);
    to_list.resize(n_features);
    weight_list.resize(n_features);
  }
  
  void operator()(std::size_t begin, std::size_t end) {
    
    for (size_t i = begin; i < end; ++i) {
      
      for (size_t j = i + 1; j < n_features; ++j) {
        
        if (std::abs(rt[i] - rt[j]) > rt_window)
          continue;
        
        double num = 0.0;
        double denom_i = 0.0;
        double denom_j = 0.0;
        
        for (int k = 0; k < n_samples; ++k) {
          double xi = X(i, k);
          double xj = X(j, k);
          num += xi * xj;
          denom_i += xi * xi;
          denom_j += xj * xj;
        }
        
        double r = num / std::sqrt(denom_i * denom_j);
        
        if (std::abs(r) < 1e-12)
          continue;
        
        //double t_stat = r * std::sqrt((n_samples - 2.0) / (1.0 - r * r));
        //double pval = 2.0 * R::pt(-std::abs(t_stat), n_samples - 2.0, 1, 0);
        
       //if (pval < alpha) {
         if(r>alpha){ 
          from_list[i].push_back(i + 1);
          to_list[i].push_back(j + 1);
          weight_list[i].push_back(r);
        }
      }
    }
  }
};

// [[Rcpp::export]]
List build_sparse_graph_parallel(NumericMatrix X,
                                 NumericVector rt,
                                 double rt_window,
                                 double alpha) {
  
  SparseGraphBuilder builder(X, rt, rt_window, alpha);
  
  parallelFor(0, X.nrow(), builder);
  
  // Merge thread-safe buffers
  std::vector<int> from;
  std::vector<int> to;
  std::vector<double> weight;
  
  for (size_t i = 0; i < builder.from_list.size(); ++i) {
    from.insert(from.end(),
                builder.from_list[i].begin(),
                builder.from_list[i].end());
    to.insert(to.end(),
              builder.to_list[i].begin(),
              builder.to_list[i].end());
    weight.insert(weight.end(),
                  builder.weight_list[i].begin(),
                  builder.weight_list[i].end());
  }
  
  return List::create(
    Named("from") = from,
    Named("to") = to,
    Named("weight") = weight
  );
}