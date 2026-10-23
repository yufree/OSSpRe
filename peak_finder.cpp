#include <Rcpp.h>
// [[Rcpp::plugins(openmp)]]
#if defined(_OPENMP) && __has_include(<omp.h>)
#include <omp.h>
#else
inline int omp_get_max_threads() { return 1; }
inline int omp_get_thread_num() { return 0; }
#endif
#include <algorithm>
#include <numeric>
#include <vector>

using namespace Rcpp;

// The Grid2D structure remains a good approach for spatial partitioning.
// No changes needed here.
struct Grid2D {
  std::vector<std::vector<std::vector<int>>> cells;
  const NumericVector& mz_values;
  const NumericVector& ccs_values;
  double mz_min, mz_max, ccs_min, ccs_max;
  int mz_bins, ccs_bins;
  
  Grid2D(int mz_b, double mz_mn, double mz_mx,
         int ccs_b, double ccs_mn, double ccs_mx,
         const NumericVector& mz_vals,
         const NumericVector& ccs_vals)
    : mz_bins(mz_b), mz_min(mz_mn), mz_max(mz_mx),
      ccs_bins(ccs_b), ccs_min(ccs_mn), ccs_max(ccs_mx),
      mz_values(mz_vals), ccs_values(ccs_vals) {
    cells.resize(mz_bins, std::vector<std::vector<int>>(ccs_bins));
  }
  
  void add_point(int idx, double mz, double ccs) {
    int mz_idx = std::min(
      static_cast<int>((mz - mz_min) / (mz_max - mz_min) * mz_bins),
      mz_bins - 1
    );
    int ccs_idx = std::min(
      static_cast<int>((ccs - ccs_min) / (ccs_max - ccs_min) * ccs_bins),
      ccs_bins - 1
    );
    mz_idx = std::max(0, mz_idx);
    ccs_idx = std::max(0, ccs_idx);
    cells[mz_idx][ccs_idx].push_back(idx);
  }
  
  std::vector<int> get_nearby_points(
      double mz, double ccs,
      double mz_window_1x, double ccs_window_1x,
      bool use_2x_window = false
  ) const {
    std::vector<int> result;
    double mz_window = use_2x_window ? 2 * mz_window_1x : mz_window_1x;
    double ccs_window = use_2x_window ? 2 * ccs_window_1x : ccs_window_1x;
    
    int mz_start = std::max(0,
                            static_cast<int>((mz - mz_window - mz_min) / (mz_max - mz_min) * mz_bins) - 1
    );
    int mz_end = std::min(mz_bins - 1,
                          static_cast<int>((mz + mz_window - mz_min) / (mz_max - mz_min) * mz_bins) + 1
    );
    int ccs_start = std::max(0,
                             static_cast<int>((ccs - ccs_window - ccs_min) / (ccs_max - ccs_min) * ccs_bins) - 1
    );
    int ccs_end = std::min(ccs_bins - 1,
                           static_cast<int>((ccs + ccs_window - ccs_min) / (ccs_max - ccs_min) * ccs_bins) + 1
    );
    
    for (int i = mz_start; i <= mz_end; ++i) {
      for (int j = ccs_start; j <= ccs_end; ++j) {
        for (int idx : cells[i][j]) {
          double mz_diff = std::abs(mz_values[idx] - mz);
          double ccs_diff = std::abs(ccs_values[idx] - ccs);
          if (mz_diff <= mz_window && ccs_diff <= ccs_window) {
            result.push_back(idx);
          }
        }
      }
    }
    return result;
  }
  
  // Zero-allocation short-circuit local maximum test.
  // In mass spectrometry datasets, >98% of points are on peak slopes or noise.
  // Directly iterating the grid cells and checking intensity >= current_intensity
  // terminates in 1-3 iterations for non-peaks without allocating any std::vector.
  bool is_local_maximum(int current_idx, double mz, double ccs, double current_intensity,
                        double mz_window, double ccs_window,
                        const NumericVector& intensity) const {
    int mz_start = std::max(0,
                            static_cast<int>((mz - mz_window - mz_min) / (mz_max - mz_min) * mz_bins) - 1
    );
    int mz_end = std::min(mz_bins - 1,
                          static_cast<int>((mz + mz_window - mz_min) / (mz_max - mz_min) * mz_bins) + 1
    );
    int ccs_start = std::max(0,
                             static_cast<int>((ccs - ccs_window - ccs_min) / (ccs_max - ccs_min) * ccs_bins) - 1
    );
    int ccs_end = std::min(ccs_bins - 1,
                           static_cast<int>((ccs + ccs_window - ccs_min) / (ccs_max - ccs_min) * ccs_bins) + 1
    );
    
    for (int i = mz_start; i <= mz_end; ++i) {
      for (int j = ccs_start; j <= ccs_end; ++j) {
        for (int idx : cells[i][j]) {
          if (idx == current_idx) continue;
          if (intensity[idx] >= current_intensity) {
            double mz_diff = std::abs(mz_values[idx] - mz);
            if (mz_diff <= mz_window) {
              double ccs_diff = std::abs(ccs_values[idx] - ccs);
              if (ccs_diff <= ccs_window) {
                return false; // Found a strictly higher/equal neighbor: not a local maximum!
              }
            }
          }
        }
      }
    }
    return true;
  }
  
  // Directly mark neighbors as checked in 1x window to suppress shoulder bins
  void mark_checked(std::vector<char>& checked, double mz, double ccs,
                     double mz_window, double ccs_window) const {
    int mz_start = std::max(0,
                            static_cast<int>((mz - mz_window - mz_min) / (mz_max - mz_min) * mz_bins) - 1
    );
    int mz_end = std::min(mz_bins - 1,
                          static_cast<int>((mz + mz_window - mz_min) / (mz_max - mz_min) * mz_bins) + 1
    );
    int ccs_start = std::max(0,
                             static_cast<int>((ccs - ccs_window - ccs_min) / (ccs_max - ccs_min) * ccs_bins) - 1
    );
    int ccs_end = std::min(ccs_bins - 1,
                           static_cast<int>((ccs + ccs_window - ccs_min) / (ccs_max - ccs_min) * ccs_bins) + 1
    );
    
    for (int i = mz_start; i <= mz_end; ++i) {
      for (int j = ccs_start; j <= ccs_end; ++j) {
        for (int idx : cells[i][j]) {
          double mz_diff = std::abs(mz_values[idx] - mz);
          if (mz_diff <= mz_window) {
            double ccs_diff = std::abs(ccs_values[idx] - ccs);
            if (ccs_diff <= ccs_window) {
              checked[idx] = 1;
            }
          }
        }
      }
    }
  }
};

struct PeakResult {
  std::vector<double> mz;
  std::vector<double> ccs;
  std::vector<double> intensity;
};

// Reverted to the more robust median-based noise calculation with O(K) nth_element.
double compute_noise_2d(const std::vector<int>& nearby_points,
                        const NumericVector& intensity,
                        double mz,
                        double ccs,
                        double mz_window_2x,
                        double mz_window_1x,
                        double ccs_window_2x,
                        double ccs_window_1x,
                        const NumericVector& mz_values,
                        const NumericVector& ccs_values) {
  std::vector<double> noise_values;
  noise_values.reserve(nearby_points.size());
  for (int j : nearby_points) {
    double mz_diff = std::abs(mz_values[j] - mz);
    double ccs_diff = std::abs(ccs_values[j] - ccs);
    
    if ((mz_diff > mz_window_1x && mz_diff <= mz_window_2x) ||
        (ccs_diff > ccs_window_1x && ccs_diff <= ccs_window_2x)) {
      noise_values.push_back(intensity[j]);
    }
  }
  
  if (noise_values.empty()) {
    return 1.0; // Return a default low noise value if no points are in the noise region
  }
  
  size_t size = noise_values.size();
  size_t mid = size / 2;
  std::nth_element(noise_values.begin(), noise_values.begin() + mid, noise_values.end());
  
  if (size % 2 == 1) {
    return noise_values[mid];
  } else {
    double val1 = noise_values[mid];
    double val2 = *std::max_element(noise_values.begin(), noise_values.begin() + mid);
    return (val1 + val2) / 2.0;
  }
}

// No changes needed in the peak checking function.
bool check_if_peak_2d(const std::vector<int>& nearby_points,
                      const NumericVector& intensity,
                      double current_intensity,
                      const NumericVector& mz_values,
                      const NumericVector& ccs_values,
                      double current_mz,
                      double current_ccs,
                      double mz_window,
                      double ccs_window) {
  for (int j : nearby_points) {
    // Do not compare the point with itself
    if (mz_values[j] == current_mz && ccs_values[j] == current_ccs) continue;
    
    double mz_diff = std::abs(mz_values[j] - current_mz);
    double ccs_diff = std::abs(ccs_values[j] - current_ccs);
    
    if (mz_diff <= mz_window && ccs_diff <= ccs_window && intensity[j] > current_intensity) {
      return false;
    }
  }
  return true;
}

// [[Rcpp::export]]
DataFrame find_2d_peaks_spatial_parallel_openmp(NumericVector mz,
                                                NumericVector ccs,
                                                NumericVector intensity,
                                                double mz_ppm = 20.0,
                                                double ccs_window_factor = 0.05,
                                                double snr_threshold = 3.0,
                                                int mz_bins = 100,
                                                int ccs_bins = 100) {
  int n = mz.length();
  if (n == 0) {
    return DataFrame::create(Named("mz") = std::vector<double>(),
                             Named("ccs") = std::vector<double>(),
                             Named("intensity") = std::vector<double>());
  }
  
  double mz_min = *std::min_element(mz.begin(), mz.end());
  double mz_max = *std::max_element(mz.begin(), mz.end());
  double ccs_min = *std::min_element(ccs.begin(), ccs.end());
  double ccs_max = *std::max_element(ccs.begin(), ccs.end());
  
  // Guard against degenerate ranges (e.g. all m/z values identical after
  // rounding) that would make the grid bin computation divide by zero
  if (!(mz_max > mz_min)) mz_max = mz_min + 1.0;
  if (!(ccs_max > ccs_min)) ccs_max = ccs_min + 1.0;
  
  Grid2D grid(mz_bins, mz_min, mz_max, ccs_bins, ccs_min, ccs_max, mz, ccs);
  for (int i = 0; i < n; ++i) {
    grid.add_point(i, mz[i], ccs[i]);
  }
  
  std::vector<int> indices(n);
  std::iota(indices.begin(), indices.end(), 0);
  std::sort(indices.begin(), indices.end(),
            [&intensity](int i1, int i2) {
              return intensity[i1] > intensity[i2];
            });
  
  // Each thread will have its own result vector to avoid race conditions.
#ifdef _OPENMP
  std::vector<PeakResult> thread_results(omp_get_max_threads());
#else
  std::vector<PeakResult> thread_results(1);
#endif
  
  std::vector<char> checked(n, 0);
  
  // The main parallel loop with critical fixes for performance and correctness.
#pragma omp parallel for schedule(dynamic, 100)
  for (int i = 0; i < n; ++i) {
    int idx = indices[i];
    
    // atomic read avoids the formal data race with the writes that happen
    // inside the critical section below
    char is_checked = 0;
#ifdef _OPENMP
#pragma omp atomic read
#endif
    is_checked = checked[idx];
    if (is_checked) {
      continue;
    }
    
    double current_mz = mz[idx];
    double current_ccs = ccs[idx];
    double current_intensity = intensity[idx];
    
    // Guard against NaN/Inf or non-positive intensity
    if (std::isnan(current_mz) || std::isnan(current_ccs) || std::isnan(current_intensity) || current_intensity <= 0.0) {
      checked[idx] = 1;
      continue;
    }
    
    double mz_window_1x = current_mz * mz_ppm * 1e-6;
    double ccs_window_1x = current_ccs * ccs_window_factor;
    
    // FAST PATH: Test 1x local maximum directly with ZERO heap allocations.
    // In mass spectrometry data, >98% of points are not local maxima.
    // Testing grid cells directly and exiting on the first higher neighbor avoids
    // allocating millions of std::vector<int> instances.
    if (!grid.is_local_maximum(idx, current_mz, current_ccs, current_intensity,
                               mz_window_1x, ccs_window_1x, intensity)) {
      checked[idx] = 1;
      continue;
    }
    
    // SLOW PATH: Only candidate local maxima compute the 2x noise median.
    double mz_window_2x = 2 * mz_window_1x;
    double ccs_window_2x = 2 * ccs_window_1x;
    
    std::vector<int> nearby_points_2x = grid.get_nearby_points(current_mz, current_ccs, mz_window_1x, ccs_window_1x, true);
    double noise_median = compute_noise_2d(nearby_points_2x, intensity, current_mz, current_ccs, mz_window_2x, mz_window_1x, ccs_window_2x, ccs_window_1x, mz, ccs);
    
    if (current_intensity < snr_threshold * noise_median) {
      checked[idx] = 1;
      continue;
    }
    
    // Confirmed peak: synchronize across threads and register
#pragma omp critical
    {
      if (!checked[idx]) {
        int thread_num = omp_get_thread_num();
        thread_results[thread_num].mz.push_back(current_mz);
        thread_results[thread_num].ccs.push_back(current_ccs);
        thread_results[thread_num].intensity.push_back(current_intensity);
        
        // Suppress shoulder bins in 1x window to prevent duplicate peaks
        grid.mark_checked(checked, current_mz, current_ccs, mz_window_1x, ccs_window_1x);
      }
    }
  }
  
  // Merge results from all threads. Pre-allocating memory is slightly more efficient.
  size_t total_peaks = 0;
  for (const auto& res : thread_results) {
    total_peaks += res.mz.size();
  }
  PeakResult final_result;
  final_result.mz.reserve(total_peaks);
  final_result.ccs.reserve(total_peaks);
  final_result.intensity.reserve(total_peaks);
  
  for (const auto& thread_result : thread_results) {
    final_result.mz.insert(final_result.mz.end(), thread_result.mz.begin(), thread_result.mz.end());
    final_result.ccs.insert(final_result.ccs.end(), thread_result.ccs.begin(), thread_result.ccs.end());
    final_result.intensity.insert(final_result.intensity.end(), thread_result.intensity.begin(), thread_result.intensity.end());
  }
  
  // Sort final results by intensity
  std::vector<size_t> sort_indices(final_result.mz.size());
  std::iota(sort_indices.begin(), sort_indices.end(), 0);
  std::sort(sort_indices.begin(), sort_indices.end(),
            [&final_result](size_t i1, size_t i2) {
              return final_result.intensity[i1] > final_result.intensity[i2];
            });
  
  std::vector<double> sorted_mz, sorted_ccs, sorted_intensity;
  sorted_mz.reserve(final_result.mz.size());
  sorted_ccs.reserve(final_result.mz.size());
  sorted_intensity.reserve(final_result.mz.size());
  
  for (size_t idx : sort_indices) {
    sorted_mz.push_back(final_result.mz[idx]);
    sorted_ccs.push_back(final_result.ccs[idx]);
    sorted_intensity.push_back(final_result.intensity[idx]);
  }
  
  return DataFrame::create(Named("mz") = sorted_mz,
                           Named("ccs") = sorted_ccs,
                           Named("intensity") = sorted_intensity);
}