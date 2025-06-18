#include <AMReX.H>
#include <AMReX_Random.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ParallelDescriptor.H>
#include <random>
#include <string>
#include <cstring>
#include <time.h>  
#include <omp.h>
#include <fstream>
#include <numeric>
#include <cassert>
#include <map>
#include <tuple>
#include <set>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>

#include "Util.H"
#include "Knapsack.H"
#include "SFC.H"
#include "SFC_knapsack.H"
#include "painterPartition.H"
#include "HilbertSFC.H"

using namespace amrex;
void main_main();

int main(int argc, char* argv[]) {
    amrex::Initialize(argc, argv);
    main_main();
    amrex::Finalize();
}
bool are_face_neighbors(const amrex::Box& a, const amrex::Box& b);
// Structure to accumulate metrics across runs
struct MetricsAccumulator {
    double total_efficiency = 0.0;
    int total_edges = 0;
    int total_local_edges = 0;
    int total_cut_edges = 0;
    std::map<int, long long> total_halo_volumes;
    int run_count = 0;
    
    void add_efficiency(double eff) {
        total_efficiency += eff;
    }
    
    void add_graph_metrics(int edges, int local, int cut) {
        total_edges += edges;
        total_local_edges += local;
        total_cut_edges += cut;
    }
    
    void add_halo_volumes(const std::map<int, long long>& volumes) {
        for (const auto& [rank, volume] : volumes) {
            total_halo_volumes[rank] += volume;
        }
    }
    
    void increment_run() { run_count++; }
    
    double get_avg_efficiency() const { 
        return run_count > 0 ? total_efficiency / run_count : 0.0; 
    }
    
    double get_avg_cut_fraction() const {
        return total_edges > 0 ? (double)total_cut_edges / total_edges : 0.0;
    }
    
    void write_averaged_graph_cuts(const std::string& filename) const {
        std::ofstream ofs(filename);
        ofs << "TotalEdges " << total_edges / run_count << "\n";
        ofs << "LocalEdges " << total_local_edges / run_count << "\n";
        ofs << "CutEdges " << total_cut_edges / run_count << "\n";
        ofs << "LocalFraction " << (double)total_local_edges / total_edges << "\n";
        ofs << "CutFraction " << (double)total_cut_edges / total_edges << "\n";
        ofs.close();
    }
    
    void write_averaged_halo_exchange(const std::string& filename) const {
        std::ofstream outfile(filename);
        for (const auto& [rank, total_volume] : total_halo_volumes) {
            long long avg_volume = run_count > 0 ? total_volume / run_count : 0;
            outfile << "Neighbor rank: " << rank << " Halo exchange volume: " << avg_volume << std::endl;
        }
        outfile.close();
    }
};

// Modified halo exchange function that returns volumes instead of writing directly
std::map<int, long long> GetHaloExchangeVolumes(amrex::MultiFab& mf) {
    mf.FillBoundary();
    
    const auto& dm = mf.DistributionMap();
    const auto& ba = mf.boxArray();
    std::map<int, long long> halo_exchange_volume;

    for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
        const amrex::Box& valid_box = mfi.validbox();
        const amrex::Box& grown_box = mfi.growntilebox();
        int my_rank = amrex::ParallelDescriptor::MyProc();

        for (int j = 0; j < ba.size(); ++j) {
            if (dm[j] == my_rank) continue;
            const amrex::Box& neighbor_box = ba[j];
            amrex::Box overlap = grown_box & neighbor_box;
            if (!overlap.isEmpty() && !valid_box.contains(overlap)) {
                halo_exchange_volume[dm[j]] += overlap.numPts();
            }
        }
    }
    return halo_exchange_volume;
}

// Modified graph cut function that returns metrics instead of writing directly
std::tuple<int, int, int> GetGraphCutMetrics(const amrex::BoxArray& ba, const amrex::DistributionMapping& dm) {
    int nboxes = ba.size();
    int cut_edges = 0, local_edges = 0, total_edges = 0;
    std::set<std::pair<int, int>> counted_edges;

    for (int i = 0; i < nboxes; ++i) {
        const amrex::Box& box_i = ba[i];
        int rank_i = dm[i];
        for (int j = 0; j < nboxes; ++j) {
            if (i == j) continue;
            const amrex::Box& box_j = ba[j];
            int rank_j = dm[j];
            if (are_face_neighbors(box_i, box_j)) {
                auto edge = std::minmax(i, j);
                if (counted_edges.count(edge)) continue;
                counted_edges.insert(edge);
                total_edges++;
                if (rank_i == rank_j) local_edges++;
                else cut_edges++;
            }
        }
    }
    return std::make_tuple(total_edges, local_edges, cut_edges);
}

// Keep original functions for backward compatibility
void LogHaloExchange(amrex::MultiFab& mf, const std::string& filename) {
    auto volumes = GetHaloExchangeVolumes(mf);
    std::ofstream outfile(filename);
    for (const auto& [rank, volume] : volumes) {
        outfile << "Neighbor rank: " << rank << " Halo exchange volume: " << volume << std::endl;
    }
    outfile.close();
}

bool are_face_neighbors(const amrex::Box& a, const amrex::Box& b) {
    int n_touch = 0;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        if (a.bigEnd(d) + 1 == b.smallEnd(d) || b.bigEnd(d) + 1 == a.smallEnd(d)) {
            n_touch++;
        } else if (a.bigEnd(d) < b.smallEnd(d) - 1 || b.bigEnd(d) < a.smallEnd(d) - 1) {
            return false;
        }
    }
    return n_touch == 1;
}

void ComputeGraphCutMetrics(const amrex::BoxArray& ba, const amrex::DistributionMapping& dm, const std::string& filename) {
    auto [total_edges, local_edges, cut_edges] = GetGraphCutMetrics(ba, dm);
    std::ofstream ofs(filename);
    ofs << "TotalEdges " << total_edges << "\n";
    ofs << "LocalEdges " << local_edges << "\n";
    ofs << "CutEdges " << cut_edges << "\n";
    ofs << "LocalFraction " << (double)local_edges / total_edges << "\n";
    ofs << "CutFraction " << (double)cut_edges / total_edges << "\n";
    ofs.close();
}

void output_distribution_data(const amrex::BoxArray& ba, 
                            const std::vector<int>& dmap,
                            const std::vector<Long>& costs,
                            const std::string& filename,
                            const std::string& algorithm_name) {
    if (amrex::ParallelDescriptor::MyProc() != 0) return;

    AMREX_ALWAYS_ASSERT(ba.size() == dmap.size());
    AMREX_ALWAYS_ASSERT(ba.size() == costs.size());

    std::ofstream outfile(filename);
    if (!outfile.is_open()) {
        amrex::Print() << "Failed to open output file: " << filename << std::endl;
        return;
    }

    outfile << "# Distribution mapping data for " << algorithm_name << "\n";
    outfile << "# Format: box_id rank x y z cost\n";
    
    amrex::Print() << "\nWriting " << algorithm_name << " distribution to " << filename << "\n";
    amrex::Print() << "Number of boxes: " << ba.size() << "\n";
    amrex::Print() << "Distribution size: " << dmap.size() << "\n";
    amrex::Print() << "Costs size: " << costs.size() << "\n";
    
    std::map<int, std::vector<std::tuple<int, int, int, int, Long>>> rank_boxes;
    
    for (int i = 0; i < ba.size(); ++i) {
        const amrex::Box& box = ba[i];
        for (int z = box.smallEnd(2); z <= box.bigEnd(2); ++z) {
            for (int y = box.smallEnd(1); y <= box.bigEnd(1); ++y) {
                for (int x = box.smallEnd(0); x <= box.bigEnd(0); ++x) {
                    rank_boxes[dmap[i]].push_back(std::make_tuple(i, x, y, z, costs[i]));
                }
            }
        }
    }
    
    amrex::Print() << "Rank distribution:\n";
    for (const auto& [rank, boxes] : rank_boxes) {
        amrex::Print() << "Rank " << rank << ": " << boxes.size() << " cells\n";
    }
    
    for (const auto& [rank, boxes] : rank_boxes) {
        for (const auto& [box_id, x, y, z, cost] : boxes) {
            outfile << box_id << " " << rank << " " 
                   << x << " " << y << " " << z << " "
                   << cost << "\n";
        }
    }
    
    outfile.close();
    amrex::Print() << "Finished writing " << filename << "\n\n";
}

void main_main() {
    BL_PROFILE("main");

    int ncomp;
    double scaling = 0.0;
    std::string name = "fb";
    int nbins, nnodes, ranks_per_node;
    Real mean, stdev;
    int nruns = 1;
    IntVect d_size, mgs, nghost, piv;
    {
        ParmParse pp;
        pp.get("domain", d_size);
        pp.get("max_grid_size", mgs);
        pp.get("ncomp", ncomp);
        pp.get("nghost", nghost);
        pp.get("periodicity", piv);
        pp.get("nnodes", nnodes); 
        pp.get("mean", mean);
        pp.get("stdev",stdev);
        pp.get("ranks_per_node", ranks_per_node);  
        pp.query("nruns",nruns);
        pp.query("name", name);
        pp.query("scaling", scaling);
    }
    amrex::Print() << "Mean: " << mean << std::endl;
    amrex::Print() << "Stdev: " << stdev << std::endl;
    amrex::Print() << "No of Run: " << nruns << std::endl;

    amrex::ResetRandomSeed(rand());

    int nmax = std::numeric_limits<int>::max();
    Real k_eff = 0.0;
    Real s_eff = 0.0;

    Box domain(IntVect{0}, (d_size -= 1));
    BoxArray ba(domain);
    ba.maxSize(mgs);
    amrex::Print() << "domain       = " << domain << '\n';
    amrex::Print() << "max_grid_size= " << mgs    << '\n';
    amrex::Print() << "boxes.size() = " << ba.size() << '\n';

    int nitems = ba.size();
    int nranks = nnodes * ranks_per_node;

    amrex::Print() << "Number of nodes: " << nnodes << std::endl;
    amrex::Print() << "Number of boxes: " << nitems << std::endl;
    amrex::Print() << "Ranks per node: " << ranks_per_node << std::endl;
    amrex::Print() << "Number of ranks: " << nranks << std::endl;

    std::vector<amrex::Real> wgts(nitems);
    std::vector<Long> bytes;

    // Initialize metric accumulators
    MetricsAccumulator knapsack_metrics, sfc_metrics, hilbert_metrics;
    
    // Store the last run's distribution data for output
    std::vector<int> final_k_dmap, final_s_dmap, final_hilbert_dmap;
    std::vector<Long> final_scaled_wgts;

    for (int r = 0; r < nruns; r++) {
        amrex::Print() << "\n=== Starting Run " << r + 1 << " ===\n";
        amrex::ResetRandomSeed(rand());

        for (int i = 0; i < nitems; ++i) {
            wgts[i] = amrex::RandomNormal(mean, stdev);
            amrex::Print() << wgts[i] << " , ";
        }
        std::vector<Long> scaled_wgts = scale_wgts(wgts);
        amrex::Print() << " Scaled Weights: ";
        for (int i = 0; i < nitems; ++i) {
            amrex::Print() << scaled_wgts[i] << " , ";
        }

        amrex::Real sfc_eff = 0.0, knapsack_eff = 0.0, hilbertsfc_eff = 0.0;
        int node_size = 0;
        double time_start = 0;
        int ng = nghost[0];

        // KNAPSACK
        time_start = amrex::second();
        std::vector<int> k_dmap = KnapSackDoIt(scaled_wgts, nranks, k_eff, true, nmax, true, false, bytes);
        amrex::Print() << " Final Knapsack time: " << amrex::second() - time_start << std::endl << std::endl;
        
        // Accumulate knapsack metrics
        knapsack_metrics.add_efficiency(k_eff);
        amrex::MultiFab mf_knapsack(ba, amrex::DistributionMapping(amrex::Vector<int>(k_dmap.begin(), k_dmap.end())), 1, ng);
        auto knapsack_halo = GetHaloExchangeVolumes(mf_knapsack);
        knapsack_metrics.add_halo_volumes(knapsack_halo);
        auto [k_total, k_local, k_cut] = GetGraphCutMetrics(ba, amrex::DistributionMapping(amrex::Vector<int>(k_dmap.begin(), k_dmap.end())));
        knapsack_metrics.add_graph_metrics(k_total, k_local, k_cut);

        // SFC
        time_start = amrex::second();
        std::vector<int> s_dmap = SFCProcessorMapDoIt(ba, scaled_wgts, nranks, &s_eff, node_size, true, false, bytes);
        amrex::Print() << " Final SFC time: " << amrex::second() - time_start << std::endl << std::endl;
        
        // Accumulate SFC metrics
        sfc_metrics.add_efficiency(s_eff);
        amrex::MultiFab mf_sfc(ba, amrex::DistributionMapping(amrex::Vector<int>(s_dmap.begin(), s_dmap.end())), 1, ng);
        auto sfc_halo = GetHaloExchangeVolumes(mf_sfc);
        sfc_metrics.add_halo_volumes(sfc_halo);
        auto [s_total, s_local, s_cut] = GetGraphCutMetrics(ba, amrex::DistributionMapping(amrex::Vector<int>(s_dmap.begin(), s_dmap.end())));
        sfc_metrics.add_graph_metrics(s_total, s_local, s_cut);

        // Other algorithms (keeping original timing)
        time_start = amrex::second();
        std::vector<int> vec = painterPartition(ba, scaled_wgts, nranks);
        amrex::Print() << " Final SFC+Painter time: " << amrex::second() - time_start << std::endl << std::endl;

        time_start = amrex::second();
        std::vector<int> sfc_knapsack_dmap = SFCProcessorMapDoItCombined(ba, scaled_wgts, nnodes, ranks_per_node, &sfc_eff, &knapsack_eff, true, false, bytes);
        amrex::Print() << " Final SFC+Knapsack_Combined time: " << amrex::second() - time_start << std::endl << std::endl;
        
        time_start = amrex::second();
        std::vector<int> painter_knapsack_dmap = SFCProcessorMapDoItCombinedPainter(ba, scaled_wgts, nnodes, ranks_per_node, &sfc_eff, &knapsack_eff, true, false, bytes);
        amrex::Print() << " Final painter+Knapsack_Combined time: " << amrex::second() - time_start << std::endl;

        // HILBERT
        time_start = amrex::second();
        std::vector<int> hilbertsfc_dmap = HilbertProcessorMapDoIt(ba, scaled_wgts, nranks, &hilbertsfc_eff, ranks_per_node, true, false, bytes);
        amrex::Print() << " Final Hilbert SFC time: " << amrex::second() - time_start << std::endl;
        
        // Accumulate Hilbert metrics
        hilbert_metrics.add_efficiency(hilbertsfc_eff);
        amrex::MultiFab mf_hilbert(ba, amrex::DistributionMapping(amrex::Vector<int>(hilbertsfc_dmap.begin(), hilbertsfc_dmap.end())), 1, ng);
        auto hilbert_halo = GetHaloExchangeVolumes(mf_hilbert);
        hilbert_metrics.add_halo_volumes(hilbert_halo);
        auto [h_total, h_local, h_cut] = GetGraphCutMetrics(ba, amrex::DistributionMapping(amrex::Vector<int>(hilbertsfc_dmap.begin(), hilbertsfc_dmap.end())));
        hilbert_metrics.add_graph_metrics(h_total, h_local, h_cut);

        // Increment run counters
        knapsack_metrics.increment_run();
        sfc_metrics.increment_run();
        hilbert_metrics.increment_run();

        // Store final run data for distribution output
        if (r == nruns - 1) {
            final_k_dmap = k_dmap;
            final_s_dmap = s_dmap;
            final_hilbert_dmap = hilbertsfc_dmap;
            final_scaled_wgts = scaled_wgts;
        }

        amrex::Print() << "\n=== End of Run " << r + 1 << " ===\n";
        amrex::Print() << "======================================\n\n";
    }

    // Write averaged results
    amrex::Print() << "\n=== WRITING AVERAGED RESULTS ===\n";
    
    // Write distribution data from final run
    output_distribution_data(ba, final_k_dmap, final_scaled_wgts, "LBC_knapsack.txt", "Knapsack");
    output_distribution_data(ba, final_s_dmap, final_scaled_wgts, "LBC_sfc.txt", "SFC");
    output_distribution_data(ba, final_hilbert_dmap, final_scaled_wgts, "LBC_hilbert.txt", "Hilbert SFC");
    
    // Write averaged metrics
    knapsack_metrics.write_averaged_graph_cuts("LBC_knapsack_graph_cut.txt");
    knapsack_metrics.write_averaged_halo_exchange("LBC_knapsack_halo_exchange.txt");
    
    sfc_metrics.write_averaged_graph_cuts("LBC_sfc_graph_cut.txt");
    sfc_metrics.write_averaged_halo_exchange("LBC_sfc_halo_exchange.txt");
    
    hilbert_metrics.write_averaged_graph_cuts("LBC_hilbert_graph_cut.txt");
    hilbert_metrics.write_averaged_halo_exchange("LBC_hilbert_halo_exchange.txt");

    // Print summary
    amrex::Print() << "\n=== AVERAGED PERFORMANCE SUMMARY ===\n";
    amrex::Print() << "Knapsack - Avg Efficiency: " << knapsack_metrics.get_avg_efficiency() 
                   << ", Avg Cut Fraction: " << knapsack_metrics.get_avg_cut_fraction() << std::endl;
    amrex::Print() << "SFC - Avg Efficiency: " << sfc_metrics.get_avg_efficiency() 
                   << ", Avg Cut Fraction: " << sfc_metrics.get_avg_cut_fraction() << std::endl;
    amrex::Print() << "Hilbert - Avg Efficiency: " << hilbert_metrics.get_avg_efficiency() 
                   << ", Avg Cut Fraction: " << hilbert_metrics.get_avg_cut_fraction() << std::endl;
}
// #include <AMReX.H>
// #include <AMReX_Random.H>
// #include <AMReX_MultiFab.H>
// #include <AMReX_ParmParse.H>
// #include <AMReX_ParallelDescriptor.H>
// #include <random>
// #include <string>
// #include <cstring>
// #include <time.h>  
// #include <omp.h>
// #include <fstream>
// #include <numeric>
// #include <cassert>
// #include <map>
// #include <tuple>
// #include <set>
// #include <AMReX_BoxArray.H>
// #include <AMReX_DistributionMapping.H>


// #include "Util.H"
// #include "Knapsack.H"
// #include "SFC.H"
// #include "SFC_knapsack.H"
// #include "painterPartition.H"
// #include "HilbertSFC.H"


// // #if defined(AMREX_USE_MPI) || defined(AMREX_USE_GPU)
// // #error This is a serial test only.
// // #endif

// using namespace amrex;
// void main_main();

// int main(int argc, char* argv[]) {
//     amrex::Initialize(argc, argv);

//     main_main();

//     amrex::Finalize();
// }

// // halo volume exchange metrics

// void LogHaloExchange(amrex::MultiFab& mf, const std::string& filename) {
//     // Fill ghost cells (this is where the exchange happens)
//     mf.FillBoundary();

//     const auto& dm = mf.DistributionMap();
//     const auto& ba = mf.boxArray();

//     // Map: neighbor rank -> number of ghost cells exchanged
//     std::map<int, long long> halo_exchange_volume;

//     // For each box owned by this rank, check which ghost cells overlap with boxes owned by other ranks
//     for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
//         const amrex::Box& valid_box = mfi.validbox(); // owned region
//         const amrex::Box& grown_box = mfi.growntilebox(); // includes ghost region
//         int my_rank = amrex::ParallelDescriptor::MyProc();

//         // Loop over all boxes to find overlaps in the ghost region
//         for (int j = 0; j < ba.size(); ++j) {
//             if (dm[j] == my_rank) continue; // skip self
//             const amrex::Box& neighbor_box = ba[j];
//             amrex::Box overlap = grown_box & neighbor_box;
//             if (!overlap.isEmpty() && !valid_box.contains(overlap)) {
//                 // This overlap is in the ghost region and comes from another rank
//                 halo_exchange_volume[dm[j]] += overlap.numPts();
//             }
//         }
//     }

//     // Write to file
//     std::ofstream outfile(filename);
//     for (const auto& [rank, volume] : halo_exchange_volume) {
//         outfile << "Neighbor rank: " << rank << " Halo exchange volume: " << volume << std::endl;
//     }
//     outfile.close();
// }

// // edge cut metrics - from AMReX_Graph.cpp
// bool are_face_neighbors(const amrex::Box& a, const amrex::Box& b) {
//     int n_touch = 0;
//     for (int d = 0; d < AMREX_SPACEDIM; ++d) {
//         if (a.bigEnd(d) + 1 == b.smallEnd(d) || b.bigEnd(d) + 1 == a.smallEnd(d)) {
//             // They are adjacent along this axis
//             n_touch++;
//         } else if (a.bigEnd(d) < b.smallEnd(d) - 1 || b.bigEnd(d) < a.smallEnd(d) - 1) {
//             // They are not adjacent nor overlapping
//             return false;
//         }
//     }
//     // They must be adjacent along exactly one axis and overlap along the others
//     return n_touch == 1;
// }

// void ComputeGraphCutMetrics(const amrex::BoxArray& ba, const amrex::DistributionMapping& dm, const std::string& filename) {
//     int nboxes = ba.size();
//     int my_rank = amrex::ParallelDescriptor::MyProc();
//     int cut_edges = 0, local_edges = 0, total_edges = 0;

//     std::set<std::pair<int, int>> counted_edges; // avoid double-counting

//     for (int i = 0; i < nboxes; ++i) {
//         const amrex::Box& box_i = ba[i];
//         int rank_i = dm[i];
//         // Loop over all possible neighbors (face neighbors)
//         for (int j = 0; j < nboxes; ++j) {
//             if (i == j) continue;
//             const amrex::Box& box_j = ba[j];
//             int rank_j = dm[j];
//             if (are_face_neighbors(box_i, box_j)) {
//                 auto edge = std::minmax(i, j);
//                 if (counted_edges.count(edge)) continue;
//                 counted_edges.insert(edge);
//                 total_edges++;
//                 if (rank_i == rank_j) local_edges++;
//                 else cut_edges++;
//             }
//         }
//     }
//     // Output
//     std::ofstream ofs(filename);
//     ofs << "TotalEdges " << total_edges << "\n";
//     ofs << "LocalEdges " << local_edges << "\n";
//     ofs << "CutEdges " << cut_edges << "\n";
//     ofs << "LocalFraction " << (double)local_edges / total_edges << "\n";
//     ofs << "CutFraction " << (double)cut_edges / total_edges << "\n";
//     ofs.close();
// }

// // distribution mapping function
// void output_distribution_data(const amrex::BoxArray& ba, 
//                             const std::vector<int>& dmap,
//                             const std::vector<Long>& costs,
//                             const std::string& filename,
//                             const std::string& algorithm_name) {
//     if (amrex::ParallelDescriptor::MyProc() != 0) return;

//     AMREX_ALWAYS_ASSERT(ba.size() == dmap.size());
//     AMREX_ALWAYS_ASSERT(ba.size() == costs.size());

//     std::ofstream outfile(filename);
//     if (!outfile.is_open()) {
//         amrex::Print() << "Failed to open output file: " << filename << std::endl;
//         return;
//     }

//     // Write header
//     outfile << "# Distribution mapping data for " << algorithm_name << "\n";
//     outfile << "# Format: box_id rank x y z cost\n";
    
//     // Debug print
//     amrex::Print() << "\nWriting " << algorithm_name << " distribution to " << filename << "\n";
//     amrex::Print() << "Number of boxes: " << ba.size() << "\n";
//     amrex::Print() << "Distribution size: " << dmap.size() << "\n";
//     amrex::Print() << "Costs size: " << costs.size() << "\n";
    
//     // Create a map to store box coordinates by rank
//     std::map<int, std::vector<std::tuple<int, int, int, int, Long>>> rank_boxes;
    
//     // First collect all boxes for each rank
//     for (int i = 0; i < ba.size(); ++i) {
//         const amrex::Box& box = ba[i];
//         for (int z = box.smallEnd(2); z <= box.bigEnd(2); ++z) {
//             for (int y = box.smallEnd(1); y <= box.bigEnd(1); ++y) {
//                 for (int x = box.smallEnd(0); x <= box.bigEnd(0); ++x) {
//                     rank_boxes[dmap[i]].push_back(std::make_tuple(i, x, y, z, costs[i]));
//                 }
//             }
//         }
//     }
    
//     // Print rank distribution
//     amrex::Print() << "Rank distribution:\n";
//     for (const auto& [rank, boxes] : rank_boxes) {
//         amrex::Print() << "Rank " << rank << ": " << boxes.size() << " cells\n";
//     }
    
//     // Write data sorted by rank
//     for (const auto& [rank, boxes] : rank_boxes) {
//         for (const auto& [box_id, x, y, z, cost] : boxes) {
//             outfile << box_id << " " << rank << " " 
//                    << x << " " << y << " " << z << " "
//                    << cost << "\n";
//         }
//     }
    
//     outfile.close();
//     amrex::Print() << "Finished writing " << filename << "\n\n";
// }

// void main_main() {
//     BL_PROFILE("main");

//     int ncomp;
//     double scaling = 0.0;
//     std::string name = "fb";
//     int nbins, nnodes, ranks_per_node;
//     Real mean, stdev;
//     int nruns = 1;
//     IntVect d_size, mgs, nghost, piv;
//     {
//         ParmParse pp;
//         pp.get("domain", d_size);
//         pp.get("max_grid_size", mgs);
//         pp.get("ncomp", ncomp);
//         pp.get("nghost", nghost);
//         pp.get("periodicity", piv);
//         // pp.get("nbins", nbins);
//         pp.get("nnodes", nnodes); 
//         pp.get("mean", mean);
//         pp.get("stdev",stdev);
//         pp.get("ranks_per_node", ranks_per_node);  
//         pp.query("nruns",nruns);
//         pp.query("name", name);
//         pp.query("scaling", scaling);
//     }
//     amrex::Print() << "Mean: " << mean << std::endl;
//     amrex::Print() << "Stdev: " << stdev << std::endl;
//     amrex::Print() << "No of Run: " << nruns << std::endl;


//     // srand(time(NULL));
//     amrex::ResetRandomSeed(rand());

//     // amrex::ResetRandomSeed(27182182459045);

//     int nmax = std::numeric_limits<int>::max();
//     Real k_eff = 0.0;
//     Real s_eff = 0.0;

//     Box domain(IntVect{0}, (d_size -= 1));
//     BoxArray ba(domain);
//     ba.maxSize(mgs);
//     amrex::Print() << "domain       = " << domain << '\n';
//     amrex::Print() << "max_grid_size= " << mgs    << '\n';
//     amrex::Print() << "boxes.size() = " << ba.size() << '\n';

//     int nitems = ba.size();
//     int nranks = nnodes * ranks_per_node;

//     amrex::Print() << "Number of nodes: " << nnodes << std::endl;
//     amrex::Print() << "Number of boxes: " << nitems << std::endl;
//     amrex::Print() << "Ranks per node: " << ranks_per_node << std::endl;
//     amrex::Print() << "Number of ranks: " << nranks << std::endl;

//     std::vector<amrex::Real> wgts(nitems);
//     std::vector<Long> bytes;

//     // Real mean = 100000;
//     // Real stdev = 4523; // average case

//     // // Real stdev = 25231; //for the worst case 

//     // // Real stdev = 250; //for the best case

// for (int r = 0; r<nruns; r++) {

//     amrex::Print() << "\n=== Starting Run " << r + 1 << " ===\n";

//     amrex::ResetRandomSeed(rand());


//     for (int i = 0; i < nitems; ++i) {
//         wgts[i] = amrex::RandomNormal(mean, stdev);
//         amrex::Print()<<wgts[i]<<" , ";
//     }
//     std::vector<Long> scaled_wgts = scale_wgts(wgts);
//     amrex::Print()<<" Scaled Weights: ";
//       for (int i=0; i<nitems; ++i) {
        
//         amrex::Print()<<scaled_wgts[i]<<" , ";
       
//     }

//     amrex::Real sfc_eff = 0.0, knapsack_eff = 0.0, hilbertsfc_eff = 0.0;
//     int node_size = 0;
//     double time_start = 0;
//     int ng = nghost[0];

//     time_start = amrex::second();
//     std::vector<int> k_dmap = KnapSackDoIt(scaled_wgts, nranks, k_eff, true, nmax, true, false, bytes);
//     amrex::Print()<<" Final Knapsack time: " << amrex::second() - time_start << std::endl<<std::endl;
//     output_distribution_data(ba, k_dmap, scaled_wgts, "LBC_knapsack.txt", "Knapsack");
//     ComputeGraphCutMetrics(ba, amrex::DistributionMapping(amrex::Vector<int>(k_dmap.begin(), k_dmap.end())), "LBC_knapsack_graph_cut.txt");
//     amrex::MultiFab mf_knapsack(ba, amrex::DistributionMapping(amrex::Vector<int>(k_dmap.begin(), k_dmap.end())), 1, ng);
//     LogHaloExchange(mf_knapsack, "LBC_knapsack_halo_exchange.txt");

//     time_start = amrex::second();
//     std::vector<int> s_dmap = SFCProcessorMapDoIt(ba, scaled_wgts, nranks, &s_eff, node_size, true, false, bytes);
//     amrex::Print()<<" Final SFC time: " << amrex::second() - time_start << std::endl<<std::endl;
//     output_distribution_data(ba, s_dmap, scaled_wgts, "LBC_sfc.txt", "SFC");
//     ComputeGraphCutMetrics(ba, amrex::DistributionMapping(amrex::Vector<int>(s_dmap.begin(), s_dmap.end())), "LBC_sfc_graph_cut.txt");
//     amrex::MultiFab mf_sfc(ba, amrex::DistributionMapping(amrex::Vector<int>(s_dmap.begin(), s_dmap.end())), 1, ng);
//     LogHaloExchange(mf_sfc, "LBC_sfc_halo_exchange.txt");

//     time_start = amrex::second();
//     std::vector<int> vec=painterPartition(ba,scaled_wgts,nranks);
//     amrex::Print()<<" Final SFC+Painter time: " << amrex::second() - time_start << std::endl<<std::endl;

//     time_start = amrex::second();
//     std::vector<int> sfc_knapsack_dmap = SFCProcessorMapDoItCombined(ba, scaled_wgts, nnodes, ranks_per_node, &sfc_eff, &knapsack_eff, true, false, bytes);
//     amrex::Print()<<" Final SFC+Knapsack_Combined time: " << amrex::second() - time_start << std::endl<<std::endl;
    
//     time_start = amrex::second();
//     std::vector<int> painter_knapsack_dmap = SFCProcessorMapDoItCombinedPainter(ba, scaled_wgts, nnodes, ranks_per_node, &sfc_eff, &knapsack_eff, true, false, bytes);
//     amrex::Print()<<" Final painter+Knapsack_Combined time: " << amrex::second() - time_start << std::endl;

//     time_start = amrex::second();
//     std::vector<int> hilbertsfc_dmap = HilbertProcessorMapDoIt(ba, scaled_wgts, nranks, &hilbertsfc_eff, ranks_per_node, true, false, bytes);
//     amrex::Print()<<" Final Hilbert SFC time: " << amrex::second() - time_start << std::endl;
//     output_distribution_data(ba, hilbertsfc_dmap, scaled_wgts, "LBC_hilbert.txt", "Hilbert SFC");
//     ComputeGraphCutMetrics(ba, amrex::DistributionMapping(amrex::Vector<int>(hilbertsfc_dmap.begin(), hilbertsfc_dmap.end())), "LBC_hilbert_graph_cut.txt");
//     amrex::MultiFab mf_hilbert(ba, amrex::DistributionMapping(amrex::Vector<int>(hilbertsfc_dmap.begin(), hilbertsfc_dmap.end())), 1, ng);
//     LogHaloExchange(mf_hilbert, "LBC_hilbert_halo_exchange.txt");
    
//     amrex::Print() << "\n=== End of Run " << r + 1 << " ===\n";
//     amrex::Print() << "======================================\n\n";

// }

//     // Print SFC and Knapsack efficiencies
//     // amrex::Print() << "SFC Efficiency: " << sfc_eff << std::endl;
//     // amrex::Print() << "Knapsack Efficiency: " << knapsack_eff << std::endl;
// }

// Heterogeneous Load Balancing
// #include <AMReX.H>
// #include <AMReX_Random.H>
// #include <AMReX_MultiFab.H>
// #include <AMReX_ParmParse.H>
// #include <AMReX_ParallelDescriptor.H>
// #include <random>
// #include <string>
// #include <cstring>
// #include <time.h>  
// #include <omp.h>
// #include <fstream>
// #include <numeric>
// #include <cassert>

// #include "Util.H"
// #include "Knapsack.H"
// #include "SFC.H"
// #include "SFC_knapsack.H"
// #include "painterPartition.H"
// #include "HilbertSFC.H"
// #include "heterogeneous/HeterogeneousLB.H"

// using namespace amrex;

// // Forward declarations
// void main_main();
// void test_heterogeneous_lb();
// extern void example_main();

// int main(int argc, char* argv[]) {
//     amrex::Initialize(argc, argv);

//     ParmParse pp;
//     std::string test_type = "main";
//     pp.query("test_type", test_type);

//     if (test_type == "example") {
//         example_main();
//     } else if (test_type == "test") {
//         test_heterogeneous_lb();
//     } else {
//         main_main();
//     }

//     amrex::Finalize();
//     return 0;
// }

// void main_main() {
//     BL_PROFILE("main");

//     ParmParse pp;
    
//     // Domain size
//     Vector<int> n_cell(AMREX_SPACEDIM);
//     if (!pp.queryarr("domain", n_cell)) {
//         n_cell[0] = n_cell[1] = n_cell[2] = 64;
//     }
    
//     // Grid size
//     Vector<int> max_grid_size(AMREX_SPACEDIM);
//     if (!pp.queryarr("max_grid_size", max_grid_size)) {
//         max_grid_size[0] = max_grid_size[1] = max_grid_size[2] = 32;
//     }
    
//     // Create domain box
//     Box domain(IntVect(0), IntVect(n_cell[0]-1, n_cell[1]-1, n_cell[2]-1));
//     BoxArray ba(domain);
//     ba.maxSize(max_grid_size[0]);
    
//     // FIXED: Proper node configuration reading with command line override
//     int nnodes = 3;  // Default
//     pp.query("nnodes", nnodes);
    
//     Vector<std::string> node_types;
//     Vector<Real> perf_factors;
//     Vector<Real> memory_caps;
    
//     // Check for command line configuration
//     bool has_custom_config = pp.queryarr("node_types", node_types) &&
//                             pp.queryarr("perf_factors", perf_factors) &&
//                             pp.queryarr("memory_caps", memory_caps);
    
//     if (has_custom_config) {
//         // Use provided configuration
//         nnodes = std::min(nnodes, static_cast<int>(std::min({node_types.size(), 
//                                                             perf_factors.size(), 
//                                                             memory_caps.size()})));
//     } else {
//         // Create configuration based on nnodes
//         node_types.clear();
//         perf_factors.clear();
//         memory_caps.clear();
        
//         if (nnodes == 4) {
//             node_types = {"cpu0", "cpu1", "cpu2", "cpu3"};
//             perf_factors = {2.0, 1.0, 0.5, 0.25};
//             memory_caps = {512.0, 256.0, 128.0, 64.0};
//         } else if (nnodes == 6) {
//             node_types = {"cpu0", "cpu1", "cpu2", "cpu3", "cpu4", "cpu5"};
//             perf_factors = {1.0, 1.0, 0.9, 0.9, 0.8, 0.8};
//             memory_caps = {256.0, 256.0, 128.0, 128.0, 128.0, 128.0};
//         } else if (nnodes == 8) {
//             node_types = {"cpu0", "cpu1", "cpu2", "cpu3", "cpu4", "cpu5", "cpu6", "cpu7"};
//             perf_factors = {2.5, 2.5, 1.5, 1.5, 1.0, 1.0, 0.7, 0.7};
//             memory_caps = {512.0, 512.0, 256.0, 256.0, 128.0, 128.0, 64.0, 64.0};
//         } else {
//             // Default 3-node heterogeneous configuration
//             nnodes = 3;
//             node_types = {"cpu0", "cpu1", "cpu2"};
//             perf_factors = {1.0, 0.9, 0.8};
//             memory_caps = {256.0, 128.0, 128.0};
//         }
//     }
    
//     // Create node info
//     Vector<ComputeNodeInfo> nodes;
//     for (int i = 0; i < nnodes; ++i) {
//         nodes.push_back({i, node_types[i], perf_factors[i], memory_caps[i]});
//     }
    
//     // Create and initialize load balancer
//     HeterogeneousLB lb;
//     lb.InitializeNodes(nodes);
    
//     // Create initial distribution mapping with correct number of processors
//     DistributionMapping dm(ba, nnodes);
    
//     // FIXED: Reasonable weight scaling to avoid huge makespan values
//     MultiFab weights(ba, dm, 1, 0);
//     for (MFIter mfi(weights); mfi.isValid(); ++mfi) {
//         const Box& box = mfi.validbox();
//         IntVect center = box.smallEnd() + box.size() / 2;
//         Real distance = std::sqrt(
//             std::pow(center[0] - n_cell[0]/2.0, 2) +
//             std::pow(center[1] - n_cell[1]/2.0, 2) +
//             std::pow(center[2] - n_cell[2]/2.0, 2)
//         );
//         // FIXED: Normalized weight scaling
//         Real normalized_distance = distance / (n_cell[0]/2.0);
//         Real weight = 1.0 + 0.5 * std::exp(-normalized_distance);  // Range [1.0, 1.5]
//         weights[mfi].setVal(weight, box, 0, 1);
//     }
    
//     // Test strategies
//     amrex::Print() << "\n" << std::string(80, '=') << "\n";
//     amrex::Print() << "HETEROGENEOUS LOAD BALANCER TEST CONFIGURATION\n";
//     amrex::Print() << std::string(80, '=') << "\n";
//     amrex::Print() << "Domain size: " << n_cell[0] << "x" << n_cell[1] << "x" << n_cell[2] << "\n";
//     amrex::Print() << "Number of boxes: " << ba.size() << "\n";
//     amrex::Print() << "Number of nodes: " << nnodes << "\n";
//     amrex::Print() << "Node configuration:\n";
//     for (const auto& node : nodes) {
//         amrex::Print() << "  Node " << node.node_id << ": " << node.node_type 
//                        << " (perf=" << std::fixed << std::setprecision(2) << node.performance_factor 
//                        << ", mem=" << node.memory_capacity << "GB)\n";
//     }
    
//     std::string strategy = "all";
//     pp.query("strategy", strategy);
    
//     Vector<std::string> strategies_to_test;
//     if (strategy == "all") {
//         strategies_to_test = {"rij", "knapsack", "sfc", "grouped"};
//     } else {
//         strategies_to_test.push_back(strategy);
//     }
    
//     Vector<LoadBalanceMetrics> all_metrics;
//     Vector<std::string> completed_strategies;
    
//     for (const auto& strat : strategies_to_test) {
//         amrex::Print() << "\n" << std::string(60, '=') << "\n";
//         amrex::Print() << "Testing " << strat << " strategy\n";
//         amrex::Print() << std::string(60, '=') << "\n";
        
//         try {
//             auto new_dm = lb.BalanceLoad(ba, weights, strat);
//             auto metrics = lb.EvaluateBalance(new_dm, weights);
//             all_metrics.push_back(metrics);
//             completed_strategies.push_back(strat);
            
//             amrex::Print() << "\nResults for " << strat << " strategy:\n";
//             amrex::Print() << "  Makespan: " << std::fixed << std::setprecision(1) << metrics.makespan << "\n";
//             amrex::Print() << "  Efficiency: " << std::fixed << std::setprecision(3) << metrics.efficiency << "\n";
//             amrex::Print() << "  Load Imbalance: " << std::fixed << std::setprecision(2) 
//                            << (metrics.load_imbalance * 100.0) << "%\n";
//             amrex::Print() << "  Total Work: " << static_cast<long>(metrics.total_work) << "\n";
            
//             lb.PrintStats();
            
//             // FIXED: More realistic performance simulation
//             Vector<double> simulated_timings;
//             std::random_device rd;
//             std::mt19937 gen(rd());
//             std::uniform_real_distribution<> dis(0.98, 1.02); // Very small variation
            
//             for (const auto& node : nodes) {
//                 double base_time = 1.0 / node.performance_factor;
//                 simulated_timings.push_back(base_time * dis(gen));
//             }
            
//             lb.UpdatePerformanceMetrics(simulated_timings);
            
//         } catch (const std::exception& e) {
//             amrex::Print() << "Error testing " << strat << " strategy: " << e.what() << "\n";
//         }
//     }
    
//     // Summary comparison
//     if (completed_strategies.size() > 1) {
//         amrex::Print() << "\n" << std::string(80, '=') << "\n";
//         amrex::Print() << "FINAL STRATEGY COMPARISON SUMMARY\n";
//         amrex::Print() << std::string(80, '=') << "\n";
//         amrex::Print() << std::setw(12) << "Strategy"
//                        << std::setw(15) << "Makespan"
//                        << std::setw(12) << "Efficiency"
//                        << std::setw(15) << "Imbalance %"
//                        << std::setw(15) << "Total Work" << "\n";
//         amrex::Print() << std::string(69, '-') << "\n";
        
//         for (size_t i = 0; i < completed_strategies.size(); ++i) {
//             amrex::Print() << std::setw(12) << completed_strategies[i]
//                            << std::setw(15) << std::fixed << std::setprecision(1) << all_metrics[i].makespan
//                            << std::setw(12) << std::fixed << std::setprecision(3) << all_metrics[i].efficiency
//                            << std::setw(15) << std::fixed << std::setprecision(2) 
//                            << (all_metrics[i].load_imbalance * 100.0)
//                            << std::setw(15) << static_cast<long>(all_metrics[i].total_work) << "\n";
//         }
        
//         // Find best strategies
//         auto best_makespan = std::min_element(all_metrics.begin(), all_metrics.end(),
//             [](const LoadBalanceMetrics& a, const LoadBalanceMetrics& b) {
//                 return a.makespan < b.makespan;
//             });
        
//         auto best_efficiency = std::max_element(all_metrics.begin(), all_metrics.end(),
//             [](const LoadBalanceMetrics& a, const LoadBalanceMetrics& b) {
//                 return a.efficiency < b.efficiency;
//             });
        
//         auto best_balance = std::min_element(all_metrics.begin(), all_metrics.end(),
//             [](const LoadBalanceMetrics& a, const LoadBalanceMetrics& b) {
//                 return a.load_imbalance < b.load_imbalance;
//             });
        
//         if (best_makespan != all_metrics.end()) {
//             size_t idx = std::distance(all_metrics.begin(), best_makespan);
//             amrex::Print() << "\nBest makespan: " << completed_strategies[idx] 
//                            << " (" << std::fixed << std::setprecision(1) << best_makespan->makespan << ")\n";
//         }
        
//         if (best_efficiency != all_metrics.end()) {
//             size_t idx = std::distance(all_metrics.begin(), best_efficiency);
//             amrex::Print() << "Best efficiency: " << completed_strategies[idx] 
//                            << " (" << std::fixed << std::setprecision(3) << best_efficiency->efficiency << ")\n";
//         }
        
//         if (best_balance != all_metrics.end()) {
//             size_t idx = std::distance(all_metrics.begin(), best_balance);
//             amrex::Print() << "Best load balance: " << completed_strategies[idx] 
//                            << " (" << std::fixed << std::setprecision(2) << (best_balance->load_imbalance * 100.0) << "%)\n";
//         }
//     }
// }
