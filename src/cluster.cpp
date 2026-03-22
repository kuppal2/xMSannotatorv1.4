// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
#include <RcppParallel.h>
#include <algorithm>
#include <cmath>
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

  std::vector<std::vector<int>>    from_list;
  std::vector<std::vector<int>>    to_list;
  std::vector<std::vector<double>> weight_list;

  SparseGraphBuilder(NumericMatrix X_, NumericVector rt_,
                     double rt_window_, double alpha_, int top_k_)
    : X(X_), rt(rt_), rt_window(rt_window_), alpha(alpha_), top_k(top_k_),
      n_samples(X_.ncol()), n_features(X_.nrow())
  {
    from_list.resize(n_features);
    to_list.resize(n_features);
    weight_list.resize(n_features);
  }

  void operator()(std::size_t begin, std::size_t end) {

    for (std::size_t i = begin; i < end; ++i) {

      std::vector<std::pair<double, int>> candidates; // (weight, j)

      for (std::size_t j = 0; j < (std::size_t)n_features; ++j) {
        if (i == j) continue;

        // RT gate first — cheap filter before computing correlation
        double rt_diff = std::abs(rt[i] - rt[j]);
        if (rt_diff > rt_window * 3.0) continue;  // skip distant features early

        // Mean-centered Pearson (X is already mean-centered)
        double num = 0.0, denom_i = 0.0, denom_j = 0.0;
        for (int k = 0; k < n_samples; ++k) {
          double xi = X(i, k);
          double xj = X(j, k);
          num     += xi * xj;
          denom_i += xi * xi;
          denom_j += xj * xj;
        }

        if (denom_i == 0.0 || denom_j == 0.0) continue;

        double r = num / std::sqrt(denom_i * denom_j);

        // Keep only POSITIVE correlations — negative = not co-pathway
        if (r < alpha) continue;

        // Soft-threshold: square the correlation (mimics WGCNA power=2)
        // then decay by RT distance
        double rt_weight = std::exp(-rt_diff / rt_window);
        double w = r * r * rt_weight;

        candidates.emplace_back(w, (int)j);
      }

      // Top-K pruning — keep only the strongest top_k edges per node
      if ((int)candidates.size() > top_k) {
        std::partial_sort(candidates.begin(),
                          candidates.begin() + top_k,
                          candidates.end(),
                          [](const std::pair<double,int>& a,
                             const std::pair<double,int>& b){
                            return a.first > b.first; // descending weight
                          });
        candidates.resize(top_k);
      }

      for (auto& p : candidates) {
        from_list[i].push_back((int)i + 1);
        to_list[i].push_back(p.second + 1);
        weight_list[i].push_back(p.first);
      }
    }
  }
};

// [[Rcpp::export]]
List build_sparse_graph_parallel(
    NumericMatrix X,
    NumericVector rt,
    double rt_window,
    double alpha  = 0.7,
    int    top_k  = 15
) {
  SparseGraphBuilder builder(X, rt, rt_window, alpha, top_k);
  parallelFor(0, X.nrow(), builder);

  std::vector<int>    from, to;
  std::vector<double> weight;

  for (std::size_t i = 0; i < builder.from_list.size(); ++i) {
    from.insert(from.end(),
                builder.from_list[i].begin(), builder.from_list[i].end());
    to.insert(to.end(),
              builder.to_list[i].begin(), builder.to_list[i].end());
    weight.insert(weight.end(),
                  builder.weight_list[i].begin(), builder.weight_list[i].end());
  }

  return List::create(Named("from") = from,
                      Named("to")   = to,
                      Named("weight") = weight);
}