// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>
#include <algorithm>

using namespace Rcpp;
using namespace RcppParallel;

// -----------------------------
// Formula checks (fast)
// -----------------------------

inline bool golden_rule_valid(const std::string& formula) {
  
  int C=0,H=0,N=0,O=0;
  
  for (size_t i=0;i<formula.size();) {
    
    if (!isupper(formula[i])) { i++; continue; }
    
    char element = formula[i++];
    int count = 0;
    
    while (i<formula.size() && isdigit(formula[i])) {
      count = count*10 + (formula[i]-'0');
      i++;
    }
    if (count==0) count=1;
    
    switch(element) {
    case 'C': C+=count; break;
    case 'H': H+=count; break;
    case 'N': N+=count; break;
    case 'O': O+=count; break;
    default: break;
    }
  }
  
  if (C<=0) return false;
  if (H > 4*C + 2) return false;
  if (N > C) return false;
  if (O > 3*C) return false;
  
  return true;
}

inline bool oxygen_present(const std::string& formula) {
  return formula.find('O') != std::string::npos;
}

// -----------------------------
// WORKER
// -----------------------------

struct AnnotationWorker : public Worker {
  
  const RVector<double> data_mz;
  const RVector<double> data_time;
  
  const std::vector<double>& db_mz;
  const std::vector<std::string>& db_id;
  const std::vector<std::string>& db_name;
  const std::vector<std::string>& db_formula;
  const std::vector<double>& db_mono_mass;
  const std::vector<std::string>& db_adduct;
  const std::vector<double>& db_adduct_mass;
  const std::vector<bool>& golden_valid;
  const std::vector<bool>& oxygen_valid;
  
  const std::vector<std::vector<int>>& adduct_index;
  
  double ppm_tol;
  
  // thread-local output
  std::vector<double> out_theo_mz;
  std::vector<std::string> out_id;
  std::vector<std::string> out_name;
  std::vector<std::string> out_formula;
  std::vector<double> out_mono_mass;
  std::vector<std::string> out_adduct;
  std::vector<double> out_adduct_mass;
  std::vector<double> out_mz;
  std::vector<double> out_time;
  
  AnnotationWorker(
    NumericVector data_mz_,
    NumericVector data_time_,
    const std::vector<double>& db_mz_,
    const std::vector<std::string>& db_id_,
    const std::vector<std::string>& db_name_,
    const std::vector<std::string>& db_formula_,
    const std::vector<double>& db_mono_mass_,
    const std::vector<std::string>& db_adduct_,
    const std::vector<double>& db_adduct_mass_,
    const std::vector<bool>& golden_valid_,
    const std::vector<bool>& oxygen_valid_,
    const std::vector<std::vector<int>>& adduct_index_,
    double ppm_tol_
  ) :
    data_mz(data_mz_),
    data_time(data_time_),
    db_mz(db_mz_),
    db_id(db_id_),
    db_name(db_name_),
    db_formula(db_formula_),
    db_mono_mass(db_mono_mass_),
    db_adduct(db_adduct_),
    db_adduct_mass(db_adduct_mass_),
    golden_valid(golden_valid_),
    oxygen_valid(oxygen_valid_),
    adduct_index(adduct_index_),
    ppm_tol(ppm_tol_)
  {}
  
  // split constructor
  AnnotationWorker(const AnnotationWorker& other, Split) :
    data_mz(other.data_mz),
    data_time(other.data_time),
    db_mz(other.db_mz),
    db_id(other.db_id),
    db_name(other.db_name),
    db_formula(other.db_formula),
    db_mono_mass(other.db_mono_mass),
    db_adduct(other.db_adduct),
    db_adduct_mass(other.db_adduct_mass),
    golden_valid(other.golden_valid),
    oxygen_valid(other.oxygen_valid),
    adduct_index(other.adduct_index),
    ppm_tol(other.ppm_tol)
  {}
  
  void operator()(std::size_t begin, std::size_t end) {
    
    for (size_t i=begin;i<end;i++) {
      
      double mz_val = data_mz[i];
      double delta  = mz_val * ppm_tol / 1e6;
      double lower  = mz_val - delta;
      double upper  = mz_val + delta;
      
      for (const auto& vec_idx : adduct_index) {
        
        auto it = std::lower_bound(
          vec_idx.begin(), vec_idx.end(), lower,
          [&](int idx, double value){
            return db_mz[idx] < value;
          });
        
        for (; it != vec_idx.end(); ++it) {
          
          int j = *it;
          if (db_mz[j] > upper) break;
          
          if (!golden_valid[j]) continue;
          
          if ((db_adduct[j]=="M+H-H2O" ||
              db_adduct[j]=="M+H-2H2O" ||
              db_adduct[j]=="M-H2O-H") &&
              !oxygen_valid[j]) continue;
              
              out_theo_mz.push_back(db_mz[j]);
              out_id.push_back(db_id[j]);
              out_name.push_back(db_name[j]);
              out_formula.push_back(db_formula[j]);
              out_mono_mass.push_back(db_mono_mass[j]);
              out_adduct.push_back(db_adduct[j]);
              out_adduct_mass.push_back(db_adduct_mass[j]);
              out_mz.push_back(mz_val);
              out_time.push_back(data_time[i]);
        }
      }
    }
  }
  
  void join(const AnnotationWorker& rhs) {
    
    out_theo_mz.insert(out_theo_mz.end(), rhs.out_theo_mz.begin(), rhs.out_theo_mz.end());
    out_id.insert(out_id.end(), rhs.out_id.begin(), rhs.out_id.end());
    out_name.insert(out_name.end(), rhs.out_name.begin(), rhs.out_name.end());
    out_formula.insert(out_formula.end(), rhs.out_formula.begin(), rhs.out_formula.end());
    out_mono_mass.insert(out_mono_mass.end(), rhs.out_mono_mass.begin(), rhs.out_mono_mass.end());
    out_adduct.insert(out_adduct.end(), rhs.out_adduct.begin(), rhs.out_adduct.end());
    out_adduct_mass.insert(out_adduct_mass.end(), rhs.out_adduct_mass.begin(), rhs.out_adduct_mass.end());
    out_mz.insert(out_mz.end(), rhs.out_mz.begin(), rhs.out_mz.end());
    out_time.insert(out_time.end(), rhs.out_time.begin(), rhs.out_time.end());
  }
};

// -----------------------------
// EXPORT
// -----------------------------

// [[Rcpp::export]]
DataFrame run_ultra_annotation_engine(
    NumericVector data_mz,
    NumericVector data_time,
    NumericVector db_mz_R,
    CharacterVector db_id_R,
    CharacterVector db_name_R,
    CharacterVector db_formula_R,
    NumericVector db_mono_mass_R,
    CharacterVector db_adduct_R,
    NumericVector db_adduct_mass_R,
    CharacterVector query_adducts,
    double ppm_tol
) {
  
  int n = db_mz_R.size();
  
  std::vector<double> db_mz(n);
  std::vector<std::string> db_id(n);
  std::vector<std::string> db_name(n);
  std::vector<std::string> db_formula(n);
  std::vector<double> db_mono_mass(n);
  std::vector<std::string> db_adduct(n);
  std::vector<double> db_adduct_mass(n);
  std::vector<bool> golden_valid(n);
  std::vector<bool> oxygen_valid(n);
  
  for (int i=0;i<n;i++) {
    
    db_mz[i] = db_mz_R[i];
    db_id[i] = Rcpp::as<std::string>(db_id_R[i]);
    db_name[i] = Rcpp::as<std::string>(db_name_R[i]);
    db_formula[i] = Rcpp::as<std::string>(db_formula_R[i]);
    db_mono_mass[i] = db_mono_mass_R[i];
    db_adduct[i] = Rcpp::as<std::string>(db_adduct_R[i]);
    db_adduct_mass[i] = db_adduct_mass_R[i];
    
    golden_valid[i] = golden_rule_valid(db_formula[i]);
    oxygen_valid[i] = oxygen_present(db_formula[i]);
  }
  
  std::vector<std::vector<int>> adduct_index;
  
  for (int a=0;a<query_adducts.size();a++) {
    
    std::string q = Rcpp::as<std::string>(query_adducts[a]);
    std::vector<int> idx;
    
    for (int i=0;i<n;i++)
      if (db_adduct[i] == q)
        idx.push_back(i);
      
      std::sort(idx.begin(), idx.end(),
                [&](int a, int b){
                  return db_mz[a] < db_mz[b];
                });
      
      adduct_index.push_back(idx);
  }
  
  AnnotationWorker worker(
      data_mz, data_time,
      db_mz, db_id, db_name, db_formula,
      db_mono_mass, db_adduct, db_adduct_mass,
      golden_valid, oxygen_valid,
      adduct_index,
      ppm_tol
  );
  
  parallelReduce(0, data_mz.size(), worker);
  
  return DataFrame::create(
    Named("theoretical.mz") = worker.out_theo_mz,
    Named("chemical_ID") = worker.out_id,
    Named("Name") = worker.out_name,
    Named("Formula") = worker.out_formula,
    Named("MonoisotopicMass") = worker.out_mono_mass,
    Named("Adduct") = worker.out_adduct,
    Named("AdductMass") = worker.out_adduct_mass,
    Named("mz") = worker.out_mz,
    Named("time") = worker.out_time
  );
}