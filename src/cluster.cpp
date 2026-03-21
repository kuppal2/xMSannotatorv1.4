// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::plugins(cpp17)]]

#include <Rcpp.h>
#include <RcppParallel.h>
#include <algorithm>

using namespace Rcpp;
using namespace RcppParallel;

struct SparseGraphBuilder : public Worker {
  
  const RMatrix<double> X;
  const RVector<double> rt;
  double rt_window;
  double alpha;
  int top_k;
  
  int n_samples;
  int n_features;
  
  std::vector< std::vector<int> > from_list;
  std::vector< std::vector<int> > to_list;
  std::vector< std::vector<double> > weight_list;
  
  SparseGraphBuilder(NumericMatrix X,
                     NumericVector rt,
                     double rt_window,
                     double alpha,
                     int top_k)
    : X(X), rt(rt), rt_window(rt_window), alpha(alpha), top_k(top_k),
      n_samples(X.ncol()), n_features(X.nrow())
  {
    from_list.resize(n_features);
    to_list.resize(n_features);
    weight_list.resize(n_features);
  }
  
  void operator()(std::size_t begin, std::size_t end) {
    
    for (size_t i = begin; i < end; ++i) {
      
      std::vector<std::pair<int, double>> neighbors;
      
      for (size_t j = 0; j < n_features; ++j) {
        if (i == j) continue;
        
        // --- compute correlation ---
        double num = 0.0, denom_i = 0.0, denom_j = 0.0;
        
        for (int k = 0; k < n_samples; ++k) {
          double xi = X(i, k);
          double xj = X(j, k);
          num += xi * xj;
          denom_i += xi * xi;
          denom_j += xj * xj;
        }
        
        if (denom_i == 0 || denom_j == 0) continue;
        
        double r = num / std::sqrt(denom_i * denom_j);
        
        if (std::abs(r) < alpha) continue;
        
        // --- soft RT weighting ---
        double rt_diff = std::abs(rt[i] - rt[j]);
        double rt_weight = std::exp(-rt_diff / rt_window);
        
        double w = std::abs(r) * rt_weight;
        
        neighbors.emplace_back(j, w);
      }
      
      // --- top-K pruning ---
      if (neighbors.size() > (size_t)top_k) {
        std::partial_sort(
          neighbors.begin(),
          neighbors.begin() + top_k,
          neighbors.end(),
          [](const std::pair<int,double>& a, const std::pair<int,double>& b) {
            return a.second > b.second;
          }
        );
        neighbors.resize(top_k);
      }
      
      // --- store edges ---
      for (auto &p : neighbors) {
        from_list[i].push_back(i + 1);
        to_list[i].push_back(p.first + 1);
        weight_list[i].push_back(p.second);
      }
    }
  }
};

// [[Rcpp::export]]
List build_sparse_graph_parallel(
    NumericMatrix X,
    NumericVector rt,
    double rt_window,
    double alpha,
    int top_k = 15
) {
  
  SparseGraphBuilder builder(X, rt, rt_window, alpha, top_k);
  
  parallelFor(0, X.nrow(), builder);
  
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