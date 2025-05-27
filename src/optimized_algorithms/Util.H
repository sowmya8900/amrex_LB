#include <AMReX_INT.H>
#include <AMReX_REAL.H>
#include <AMReX_BoxArray.H>
#include <AMReX_BoxList.H>
#include <type_traits>
#include <vector>
#include <memory>

template <class T, std::enable_if_t<std::is_arithmetic<T>::value, bool> = true>
std::vector<amrex::Long> scale_wgts(std::vector<T> wgts) {
    double wmax = *std::max_element(wgts.begin(), wgts.end());
    double scale = (wmax == 0) ? double(1e9) : double(1.e9) / wmax;
    std::vector<amrex::Long> swgts(wgts.size());
    for (long unsigned int i = 0; i < wgts.size(); ++i) {
        swgts[i] = amrex::Long(wgts[i] * scale) + 1L;
    }
    return swgts;
}

// inline amrex::Real computeEfficiency(const std::vector<amrex::Long>& wgts,
//                                    const std::vector<int>& dmap,
//                                    int nranks) {
//     BL_PROFILE("computeEfficiency");
    
//     if (wgts.empty() || nranks <= 0) return 0.0;

//     std::vector<amrex::Long> rank_weights(nranks, 0);
//     amrex::Long total_weight = 0;

//     for (size_t i = 0; i < wgts.size(); ++i) {
//         int rank = dmap[i];
//         if (rank >= 0 && rank < nranks) {
//             rank_weights[rank] += wgts[i];
//         }
//         total_weight += wgts[i];
//     }

//     if (total_weight == 0) return 0.0;

//     amrex::Long max_weight = *std::max_element(rank_weights.begin(), rank_weights.end());
//     amrex::Real avg_weight = static_cast<amrex::Real>(total_weight) / nranks;
    
//     return (max_weight > 0) ? avg_weight / max_weight : 0.0;
// }

// inline double calculateCommCost(const amrex::BoxArray& ba, const std::vector<int>& dmap) {
//     BL_PROFILE("calculateCommCost");
//     double total_cost = 0.0;
//     const int nboxes = ba.size();

//     for (int i = 0; i < nboxes; ++i) {
//         amrex::Box bx = ba[i]; // Create a non-const copy
//         bx.grow(1); // Grow the box first
        
//         // Get intersecting boxes and their indices
//         auto neighbors = ba.intersections(bx);

//         for (const auto& neighbor_pair : neighbors) {
//             int j = neighbor_pair.first;
//             const amrex::Box& neighbor = neighbor_pair.second;
            
//             if (i == j) continue; // Skip self
            
//             if (dmap[i] != dmap[j]) {
//                 const amrex::Box intersect = bx & neighbor;
//                 const amrex::IntVect size = intersect.size();

//                 if (size[0] == 1) { // x-face
//                     total_cost += size[1] * size[2];
//                 } else if (size[1] == 1) { // y-face
//                     total_cost += size[0] * size[2];
//                 } else if (size[2] == 1) { // z-face
//                     total_cost += size[0] * size[1];
//                 }
//             }
//         }
//     }
//     return total_cost;
// }