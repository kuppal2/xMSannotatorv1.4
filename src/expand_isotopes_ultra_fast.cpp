// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::depends(RcppParallel)]]

#include <Rcpp.h>
#include <RcppParallel.h>
#include <cmath>
#include <vector>
#include <string>
#include <unordered_map>

using namespace Rcpp;
using namespace RcppParallel;


// isotope constants
static const double C13  = 1.0033548378;
static const double Cl37 = 1.99704989;
static const double S34  = 1.9957959;
static const double Br81 = 1.9979535;


// extract element count
inline int extract_count(const std::string& formula,
                         const std::string& element)
{
  size_t pos = formula.find(element);
  
  if (pos == std::string::npos)
    return 0;
  
  pos += element.length();
  
  int count = 0;
  
  while (pos < formula.size() && isdigit(formula[pos])) {
    count = count * 10 + (formula[pos] - '0');
    pos++;
  }
  
  return (count == 0 ? 1 : count);
}



// -----------------------------
// Worker
// -----------------------------

struct IsoWorker : public Worker {
  
  const RVector<double> mz;
  const RVector<double> time;
  const RVector<double> mono_mz;
  const RVector<int> module_id;
  
  const RVector<int> group_start;
  const RVector<int> group_end;
  
  const RVector<double> isp_mz;
  const RVector<double> isp_time;
  
  const std::unordered_map<int,std::vector<int>>& isp_index;
  
  const CharacterVector formula;
  
  double max_diff_rt;
  int max_isp;
  double ppm_tol;
  
  std::vector<int> parent;
  std::vector<double> out_mz;
  std::vector<double> out_time;
  std::vector<std::string> out_formula;
  std::vector<std::string> out_adduct;
  
  
  IsoWorker(
    NumericVector mz_,
    NumericVector time_,
    NumericVector mono_mz_,
    IntegerVector module_id_,
    IntegerVector group_start_,
    IntegerVector group_end_,
    NumericVector isp_mz_,
    NumericVector isp_time_,
    const std::unordered_map<int,std::vector<int>>& isp_index_,
    CharacterVector formula_,
    double max_diff_rt_,
    int max_isp_,
    double ppm_tol_)
    :
    mz(mz_),
    time(time_),
    mono_mz(mono_mz_),
    module_id(module_id_),
    group_start(group_start_),
    group_end(group_end_),
    isp_mz(isp_mz_),
    isp_time(isp_time_),
    isp_index(isp_index_),
    formula(formula_),
    max_diff_rt(max_diff_rt_),
    max_isp(max_isp_),
    ppm_tol(ppm_tol_)
  {}
  
  
  
  IsoWorker(const IsoWorker& other, Split)
    :
    mz(other.mz),
    time(other.time),
    mono_mz(other.mono_mz),
    module_id(other.module_id),
    group_start(other.group_start),
    group_end(other.group_end),
    isp_mz(other.isp_mz),
    isp_time(other.isp_time),
    isp_index(other.isp_index),
    formula(other.formula),
    max_diff_rt(other.max_diff_rt),
    max_isp(other.max_isp),
    ppm_tol(other.ppm_tol)
  {}
  
  
  
  void operator()(std::size_t begin,
                std::size_t end)
  {
    
    for (size_t g = begin; g < end; g++) {
      
      int start = group_start[g];
      int stop  = group_end[g];
      
      if (stop <= start)
        continue;
      
      int mod = module_id[start];
      
      auto it = isp_index.find(mod);
      if (it == isp_index.end())
        continue;
      
      const std::vector<int>& isp_rows = it->second;
      
      std::string f =
        Rcpp::as<std::string>(formula[start]);
      
      int nC  = extract_count(f,"C");
      int nCl = extract_count(f,"Cl");
      int nBr = extract_count(f,"Br");
      int nS  = extract_count(f,"S");
      
      for (int i = start; i < stop; i++) {
        
        double mono = mono_mz[i];
        double rt   = time[i];
        
        double ppm_scale =
          mono * ppm_tol * 1e-6;
        
        for (int idx : isp_rows) {
          
          if (std::abs(isp_time[idx] - rt) > max_diff_rt)
            continue;
          
          double shift =
            isp_mz[idx] - mono;
          
          if (shift <= 0)
            continue;
          
          
          int c13 = std::round(shift / C13);
          
          if (c13 > 0 &&
              c13 <= std::min(nC,max_isp) &&
              std::abs(shift - c13*C13) <= ppm_scale)
          {
            
            parent.push_back(i);
            out_mz.push_back(isp_mz[idx]);
            out_time.push_back(isp_time[idx]);
            
            std::string label =
              "M_[+" + std::to_string(c13) + "]";
            
            out_adduct.push_back(label);
            out_formula.push_back(
              f + "_[+" + std::to_string(c13) + "]"
            );
            
            continue;
          }
          
          
          
          int cl = std::round(shift / Cl37);
          
          if (cl > 0 &&
              cl <= std::min(nCl,max_isp) &&
              std::abs(shift - cl*Cl37) <= ppm_scale)
          {
            
            parent.push_back(i);
            out_mz.push_back(isp_mz[idx]);
            out_time.push_back(isp_time[idx]);
            
            std::string label =
              "M_[+" + std::to_string(cl) + "]";
            
            out_adduct.push_back(label);
            out_formula.push_back(
              f + "_[+" + std::to_string(cl) + "]"
            );
            
            continue;
          }
          
          
          
          int br = std::round(shift / Br81);
          
          if (br > 0 &&
              br <= std::min(nBr,max_isp) &&
              std::abs(shift - br*Br81) <= ppm_scale)
          {
            
            parent.push_back(i);
            out_mz.push_back(isp_mz[idx]);
            out_time.push_back(isp_time[idx]);
            
            std::string label =
              "M_[+" + std::to_string(br) + "]";
            
            out_adduct.push_back(label);
            out_formula.push_back(
              f + "_[+" + std::to_string(br) + "]"
            );
            
            continue;
          }
          
          
          
          int s = std::round(shift / S34);
          
          if (s > 0 &&
              s <= std::min(nS,max_isp) &&
              std::abs(shift - s*S34) <= ppm_scale)
          {
            
            parent.push_back(i);
            out_mz.push_back(isp_mz[idx]);
            out_time.push_back(isp_time[idx]);
            
            std::string label =
              "M_[+" + std::to_string(s) + "]";
            
            out_adduct.push_back(label);
            out_formula.push_back(
              f + "_[+" + std::to_string(s) + "]"
            );
          }
          
        }
      }
    }
  }
  
  
  
  void join(const IsoWorker& rhs)
  {
    parent.insert(parent.end(),
                  rhs.parent.begin(),
                  rhs.parent.end());
    
    out_mz.insert(out_mz.end(),
                  rhs.out_mz.begin(),
                  rhs.out_mz.end());
    
    out_time.insert(out_time.end(),
                    rhs.out_time.begin(),
                    rhs.out_time.end());
    
    out_formula.insert(out_formula.end(),
                       rhs.out_formula.begin(),
                       rhs.out_formula.end());
    
    out_adduct.insert(out_adduct.end(),
                      rhs.out_adduct.begin(),
                      rhs.out_adduct.end());
  }
};




// [[Rcpp::export]]

DataFrame expand_isotopes_ultra_fast_cpp(
    
    NumericVector mz,
    NumericVector time,
    NumericVector mono_mz,
    IntegerVector module_id,
    
    IntegerVector group_start,
    IntegerVector group_end,
    
    NumericVector isp_mz,
    NumericVector isp_time,
    IntegerVector isp_module,
    
    CharacterVector formula,
    
    double max_diff_rt,
    int max_isp,
    double ppm_tol
)
{
  
  std::unordered_map<int,std::vector<int>> isp_index;
  
  for (int i=0;i<isp_module.length();i++)
    isp_index[ isp_module[i] ].push_back(i);
  
  
  
  IsoWorker worker(
      mz,
      time,
      mono_mz,
      module_id,
      group_start,
      group_end,
      isp_mz,
      isp_time,
      isp_index,
      formula,
      max_diff_rt,
      max_isp,
      ppm_tol
  );
  
  
  parallelReduce(
    0,
    group_start.length(),
    worker
  );
  
  
  return DataFrame::create(
    Named("parent_index") = worker.parent,
    Named("mz") = worker.out_mz,
    Named("time") = worker.out_time,
    Named("Formula") = worker.out_formula,
    Named("Adduct") = worker.out_adduct
  );
}