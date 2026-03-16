// [[Rcpp::depends(RcppParallel)]]

#include <Rcpp.h>
#include <RcppParallel.h>
#include <unordered_map>
#include <algorithm>
#include <set>
#include <limits>

using namespace Rcpp;
using namespace RcppParallel;

struct Stage3Worker : public Worker {
  
  const RVector<double> time;
  const RVector<double> intensity;
  const RVector<double> charge;
  const CharacterVector adduct;
  
  const RVector<int> group_start;
  const RVector<int> group_end;
  
  const std::unordered_map<std::string,double>& weight_map;
  const std::vector<std::string>& primary_adducts;
  
  double max_diff_rt;
  
  RVector<double> score;
  
  Stage3Worker(
    NumericVector time_,
    NumericVector intensity_,
    NumericVector charge_,
    CharacterVector adduct_,
    IntegerVector group_start_,
    IntegerVector group_end_,
    std::unordered_map<std::string,double>& weight_map_,
    std::vector<std::string>& primary_adducts_,
    double max_diff_rt_,
    NumericVector score_)
    :
    time(time_),
    intensity(intensity_),
    charge(charge_),
    adduct(adduct_),
    group_start(group_start_),
    group_end(group_end_),
    weight_map(weight_map_),
    primary_adducts(primary_adducts_),
    max_diff_rt(max_diff_rt_),
    score(score_) {}
  
  Stage3Worker(const Stage3Worker& other, Split)
    :
    time(other.time),
    intensity(other.intensity),
    charge(other.charge),
    adduct(other.adduct),
    group_start(other.group_start),
    group_end(other.group_end),
    weight_map(other.weight_map),
    primary_adducts(other.primary_adducts),
    max_diff_rt(other.max_diff_rt),
    score(other.score) {}
  
  void operator()(std::size_t begin, std::size_t end) {
    
    for (size_t g = begin; g < end; g++) {
      
      int start = group_start[g];
      int stop  = group_end[g];
      
      if (stop <= start) continue;
      
      //-------------------------------------
      // 1. Primary RT anchor
      //-------------------------------------
      
      std::vector<double> primary_times;
      
      for (int i = start; i < stop; i++) {
        
        std::string add = Rcpp::as<std::string>(adduct[i]);
        
        for (auto& p : primary_adducts)
          if (add == p)
            primary_times.push_back(time[i]);
      }
      
      double primary_rt = NA_REAL;
      double primary_rt_mean = NA_REAL;
      
      if (!primary_times.empty()) {
        
        std::sort(primary_times.begin(), primary_times.end());
        
        size_t idx = floor(0.75 * (primary_times.size() - 1));
        primary_rt = primary_times[idx];
        
        // --- NEW: mean RT of primary adduct matches
        double sum_rt = 0.0;
        for (double v : primary_times)
          sum_rt += v;
        primary_rt_mean = sum_rt / primary_times.size();
      }
      
      //-------------------------------------
      // 2. STRICT RT FILTERING
      //-------------------------------------
      
      std::vector<int> keep;
      
      for (int i = start; i < stop; i++) {
        
        bool keep_rt = true;
        
        if (!NumericVector::is_na(primary_rt)) {
          if (std::abs(time[i] - primary_rt) > max_diff_rt)
            keep_rt = false;
        }
        
        // --- NEW: filter against primary mean RT
        if (!NumericVector::is_na(primary_rt_mean)) {
          if (std::abs(time[i] - primary_rt_mean) > max_diff_rt)
            keep_rt = false;
        }
        
        if (!keep_rt)
          continue;
        
        keep.push_back(i);
      }
      
      if (keep.empty()) continue;
      
      //-------------------------------------
      // 3. RT RANGE (IQR)
      //-------------------------------------
      
      std::vector<double> rt_vals;
      
      for (int idx : keep)
        rt_vals.push_back(time[idx]);
      
      std::sort(rt_vals.begin(), rt_vals.end());
      
      int n = rt_vals.size();
      
      double q1 = rt_vals[(int)floor(0.25 * (n - 1))];
      double q3 = rt_vals[(int)floor(0.75 * (n - 1))];
      
      double rt_range = q3 - q1;
      
      double sigma = (rt_range <= max_diff_rt) ? 1.0 : 10.0;
      
      if (rt_range <= max_diff_rt)
        rt_range = 1.0;
      
      //-------------------------------------
      // 4. Adduct diversity
      //-------------------------------------
      
      std::set<std::string> adduct_set;
      
      for (int idx : keep)
        adduct_set.insert(Rcpp::as<std::string>(adduct[idx]));
      
      int n_adducts = adduct_set.size();
      
      //-------------------------------------
      // 5. Adduct weights
      //-------------------------------------
      
      double adduct_sum = 0.0;
      
      for (const auto& a : adduct_set) {
        
        auto it = weight_map.find(a);
        
        if (it != weight_map.end())
          adduct_sum += it->second;
      }
      
      //-------------------------------------
      // 6. Multimer + charge penalties
      //-------------------------------------
      
      double mono_int =
        -std::numeric_limits<double>::infinity();
        
        double multi_int =
        -std::numeric_limits<double>::infinity();
        
        double charge1_int =
        -std::numeric_limits<double>::infinity();
        
        double high_charge_int =
        -std::numeric_limits<double>::infinity();
        
        for (int idx : keep) {
          
          std::string add = Rcpp::as<std::string>(adduct[idx]);
          double inten = intensity[idx];
          
          bool multimer =
            add.rfind("2M",0)==0 ||
            add.rfind("3M",0)==0;
          
          if (!multimer)
            mono_int = std::max(mono_int,inten);
          else
            multi_int = std::max(multi_int,inten);
          
          if (charge[idx] == 1)
            charge1_int = std::max(charge1_int,inten);
          
          if (charge[idx] > 1)
            high_charge_int = std::max(high_charge_int,inten);
        }
        
        double multimer_penalty = 1.0;
        
        if (multi_int > mono_int)
          multimer_penalty = 0.5;
        
        if (high_charge_int > charge1_int)
          multimer_penalty *= 0.5;
        
        //-------------------------------------
        // 7. Final score
        //-------------------------------------
        
        double chemical_score =
          multimer_penalty *
          (1.0 / (rt_range * sigma + 1.0)) *
          n_adducts *
          10.0 *
          adduct_sum;
        
        for (int idx : keep)
          score[idx] = chemical_score;
    }
  }
};


// [[Rcpp::export]]
DataFrame stage3_score_engine_cpp(
    
    DataFrame DT,
    DataFrame adduct_weights,
    CharacterVector primary_adducts,
    double max_diff_rt,
    IntegerVector group_start,
    IntegerVector group_end
) {
  
  NumericVector time = DT["time"];
  NumericVector intensity = DT["mean_int_vec"];
  
  NumericVector charge;
  
  if (DT.containsElementNamed("charge"))
    charge = DT["charge"];
  else
    charge = NumericVector(DT.nrows(),1.0);
  
  CharacterVector adduct = DT["Adduct"];
  
  //-------------------------------------
  // Build weight map
  //-------------------------------------
  
  std::unordered_map<std::string,double> weight_map;
  
  CharacterVector aw_name = adduct_weights[0];
  NumericVector aw_val = adduct_weights[1];
  
  for (int i=0;i<aw_name.size();i++)
    weight_map[Rcpp::as<std::string>(aw_name[i])] = aw_val[i];
  
  std::vector<std::string> prim;
  
  for (int i=0;i<primary_adducts.size();i++)
    prim.push_back(Rcpp::as<std::string>(primary_adducts[i]));
  
  //NumericVector score(DT.nrows());
  NumericVector score(DT.nrows(), NA_REAL);
  
  Stage3Worker worker(
      time,
      intensity,
      charge,
      adduct,
      group_start,
      group_end,
      weight_map,
      prim,
      max_diff_rt,
      score
  );
  
  parallelFor(0, group_start.size(), worker);
  
  DT["score"] = score;
  
  return DT;
}