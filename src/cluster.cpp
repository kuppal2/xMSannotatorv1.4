// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
#include <RcppParallel.h>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <vector>
#include <unordered_map>

using namespace Rcpp;
using namespace RcppParallel;

// ============================================================
// STEP 1 — Sparse graph (existing, cleaned up)
// ============================================================

struct SparseGraphBuilder : public Worker {

  const RMatrix<double> X;
  const RVector<double> rt;
  double rt_window;
  double alpha;
  int    top_k;
  int    n_samples;
  int    n_features;

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

      std::vector<std::pair<double,int>> candidates;

      for (std::size_t j = 0; j < (std::size_t)n_features; ++j) {
        if (i == j) continue;

        // RT gate — cheap early exit
        double rt_diff = std::abs(rt[i] - rt[j]);
        if (rt_diff > rt_window * 3.0) continue;

        // Pearson correlation (X already mean-centered)
        double num=0, di=0, dj=0;
        for (int k = 0; k < n_samples; ++k) {
          double xi = X(i,k), xj = X(j,k);
          num += xi*xj; di += xi*xi; dj += xj*xj;
        }
        if (di == 0.0 || dj == 0.0) continue;

        double r = num / std::sqrt(di * dj);
        if (r < alpha) continue;  // positive only

        double rt_weight = std::exp(-rt_diff / rt_window);
        double w = r * r * rt_weight;
        candidates.emplace_back(w, (int)j);
      }

      // Top-K pruning
      if ((int)candidates.size() > top_k) {
        std::partial_sort(
          candidates.begin(),
          candidates.begin() + top_k,
          candidates.end(),
          [](const std::pair<double,int>& a, const std::pair<double,int>& b){
            return a.first > b.first;
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
    double alpha  = 0.70,
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

  return List::create(Named("from")=from, Named("to")=to, Named("weight")=weight);
}


// ============================================================
// STEP 2 — Sparse TOM in C++
// Computes TOM only for edges that exist in the sparse graph.
// Never allocates a dense n×n matrix.
//
// TOM_ij = (adj_ij + shared_ij) / (min(deg_i, deg_j) + 1 - adj_ij)
// where shared_ij = Σ_u adj_iu * adj_uj  (common neighbor strength)
// ============================================================

// [[Rcpp::export]]
NumericVector compute_sparse_tom(
    IntegerVector from,     // 1-based
    IntegerVector to,       // 1-based
    NumericVector weight,   // already soft-thresholded
    int           n         // total number of features
) {
  int n_edges = from.size();

  // ── Build adjacency lookup per node ─────────────────────
  // adj_map[i] = map of {j -> weight} for all edges from i
  // (store both directions)
  std::vector<std::unordered_map<int,double>> adj(n + 1);
  std::vector<double> degree(n + 1, 0.0);

  for (int e = 0; e < n_edges; ++e) {
    int    i = from[e], j = to[e];
    double w = weight[e];
    adj[i][j]  = w;
    adj[j][i]  = w;    // symmetrise
    degree[i] += w;
    degree[j] += w;
  }

  // ── Compute TOM for each edge ────────────────────────────
  NumericVector tom_weight(n_edges);

  for (int e = 0; e < n_edges; ++e) {
    int    i    = from[e];
    int    j    = to[e];
    double a_ij = weight[e];

    // Shared neighbor sum: iterate over the SHORTER neighbor list
    // for efficiency — O(min(deg_i, deg_j))
    double shared = 0.0;

    const auto& ni = adj[i];
    const auto& nj = adj[j];

    if (ni.size() <= nj.size()) {
      for (const auto& kv : ni) {
        int u = kv.first;
        if (u == j) continue;
        auto it = nj.find(u);
        if (it != nj.end()){
          shared += kv.second * it->second;
        }
      }
    } else {
      for (const auto& kv : nj) {
        int u = kv.first;
        if (u == i) continue;
        auto it = ni.find(u);
        if (it != ni.end()){
          shared += kv.second * it->second;
        }
      }
    }

    double denom = std::min(degree[i], degree[j]) - a_ij + 1.0;
    tom_weight[e] = (a_ij + shared) / denom;
  }

  return tom_weight;
}


// ============================================================
// STEP 3 — Infomap clustering on sparse TOM graph (pure C++)
//
// We implement a lightweight version of the map equation
// minimisation via a greedy single-level flow partition.
// For production use, igraph's C-level cluster_infomap is
// superior — this is provided for cases where igraph cannot
// be used (embedded, shinyapps constraints, etc.).
//
// For most users: call igraph::cluster_infomap() in R after
// getting tom_weights from compute_sparse_tom() above.
// That is the recommended path.
// ============================================================

// Lightweight greedy modularity (Louvain-style) as fallback
// Full infomap requires a more complex C++ implementation

struct GreedyCommunity {

  int n;
  std::vector<int> membership;
  std::vector<std::vector<std::pair<int,double>>> adj;
  std::vector<double> vol;   // sum of edge weights per node
  double total_weight;

  GreedyCommunity(int n_, const IntegerVector& from,
                  const IntegerVector& to, const NumericVector& w)
    : n(n_), membership(n_+1), adj(n_+1), vol(n_+1, 0.0), total_weight(0.0)
  {
    std::iota(membership.begin(), membership.end(), 0);

    for (int e = 0; e < from.size(); ++e) {
      int i=from[e], j=to[e]; double wt=w[e];
      adj[i].emplace_back(j, wt);
      adj[j].emplace_back(i, wt);
      vol[i] += wt; vol[j] += wt;
      total_weight += wt;
    }
  }

  // Compute modularity gain of moving node i to community c
  double delta_modularity(int i, int c) {
    double k_i_in = 0.0;
    double k_c    = 0.0;
    double k_i    = vol[i];

    // Sum weights from i to neighbours in community c
    for (auto& kv : adj[i]) {
      if (membership[kv.first] == c) k_i_in += kv.second;
    }
    // Volume of community c
    for (int u = 1; u <= n; ++u) {
      if (membership[u] == c && u != i) k_c += vol[u];
    }

    double m2 = 2.0 * total_weight;
    return 2.0 * k_i_in / m2 - 2.0 * k_i * k_c / (m2 * m2);
  }

  bool one_pass() {
    bool improved = false;
    for (int i = 1; i <= n; ++i) {
      int    best_c     = membership[i];
      double best_delta = 0.0;

      std::unordered_map<int,double> neighbour_communities;
      for (auto& kv : adj[i]) {
        int c = membership[kv.first];
        if (c != membership[i]) neighbour_communities[c] += kv.second;
      }

      for (auto& mc : neighbour_communities) {
        double dq = delta_modularity(i, mc.first);
        if (dq > best_delta) { best_delta = dq; best_c = mc.first; }
      }

      if (best_c != membership[i]) {
        membership[i] = best_c;
        improved = true;
      }
    }
    return improved;
  }

  void run(int max_iter = 20) {
    for (int it = 0; it < max_iter; ++it) {
      if (!one_pass()) break;
    }
    // Relabel to 1..K
    std::unordered_map<int,int> remap;
    int next = 1;
    for (int i = 1; i <= n; ++i) {
      int c = membership[i];
      if (remap.find(c) == remap.end()) remap[c] = next++;
      membership[i] = remap[c];
    }
  }
};

// [[Rcpp::export]]
IntegerVector cluster_sparse_greedy(
    IntegerVector from,
    IntegerVector to,
    NumericVector weight,
    int           n,
    int           max_iter = 20
) {
  GreedyCommunity gc(n, from, to, weight);
  gc.run(max_iter);

  IntegerVector result(n);
  for (int i = 1; i <= n; ++i) result[i-1] = gc.membership[i];
  return result;
}