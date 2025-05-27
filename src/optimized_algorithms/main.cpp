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


#include "Util.H"
#include "Knapsack.H"
#include "SFC.H"
#include "SFC_knapsack.H"
#include "painterPartition.H"
#include "HilbertSFC.H"


// #if defined(AMREX_USE_MPI) || defined(AMREX_USE_GPU)
// #error This is a serial test only.
// #endif

using namespace amrex;
void main_main();

int main(int argc, char* argv[]) {
    amrex::Initialize(argc, argv);

    main_main();

    amrex::Finalize();
}

void output_distribution_data(const amrex::BoxArray& ba, 
                            const std::vector<int>& dmap,
                            const std::vector<Long>& costs,
                            const std::string& filename,
                            const std::string& algorithm_name) {
    std::ofstream outfile(filename);
    if (!outfile.is_open()) {
        amrex::Print() << "Failed to open output file: " << filename << std::endl;
        return;
    }

    // Write header
    outfile << "# Distribution mapping data for " << algorithm_name << "\n";
    outfile << "# Format: box_id rank x y z cost\n";
    
    // Debug print
    amrex::Print() << "\nWriting " << algorithm_name << " distribution to " << filename << "\n";
    amrex::Print() << "Number of boxes: " << ba.size() << "\n";
    amrex::Print() << "Distribution size: " << dmap.size() << "\n";
    amrex::Print() << "Costs size: " << costs.size() << "\n";
    
    // Create a map to store box coordinates by rank
    std::map<int, std::vector<std::tuple<int, int, int, int, Long>>> rank_boxes;
    
    // First collect all boxes for each rank
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
    
    // Print rank distribution
    amrex::Print() << "Rank distribution:\n";
    for (const auto& [rank, boxes] : rank_boxes) {
        amrex::Print() << "Rank " << rank << ": " << boxes.size() << " cells\n";
    }
    
    // Write data sorted by rank
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
        // pp.get("nbins", nbins);
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


    // srand(time(NULL));
    amrex::ResetRandomSeed(rand());

    // amrex::ResetRandomSeed(27182182459045);

    int nmax = std::numeric_limits<int>::max();
    Real k_eff = 0.0;
    Real s_eff = 0.0;

    Box domain(IntVect{0}, (d_size -= 1));
    BoxArray ba(domain);
    ba.maxSize(mgs);

    int nitems = ba.size();
    int nranks = nnodes * ranks_per_node;

    amrex::Print() << "Number of nodes: " << nnodes << std::endl;
    amrex::Print() << "Number of boxes: " << nitems << std::endl;
    amrex::Print() << "Ranks per node: " << ranks_per_node << std::endl;
    amrex::Print() << "Number of ranks: " << nranks << std::endl;

    std::vector<amrex::Real> wgts(nitems);
    std::vector<Long> bytes;

    // Real mean = 100000;
    // Real stdev = 4523; // average case

    // // Real stdev = 25231; //for the worst case 

    // // Real stdev = 250; //for the best case

for (int r = 0; r<nruns; r++) {

    amrex::Print() << "\n=== Starting Run " << r + 1 << " ===\n";

    amrex::ResetRandomSeed(rand());


    for (int i = 0; i < nitems; ++i) {
        wgts[i] = amrex::RandomNormal(mean, stdev);
        amrex::Print()<<wgts[i]<<" , ";
    }
    std::vector<Long> scaled_wgts = scale_wgts(wgts);
    amrex::Print()<<" Scaled Weights: ";
      for (int i=0; i<nitems; ++i) {
        
        amrex::Print()<<scaled_wgts[i]<<" , ";
       
    }

    amrex::Real sfc_eff = 0.0, knapsack_eff = 0.0, hilbertsfc_eff = 0.0;
    int node_size = 0;
    double time_start=0;

    time_start = amrex::second();
    std::vector<int> k_dmap = KnapSackDoIt(scaled_wgts, nranks, k_eff, true, nmax, true, false, bytes);
    amrex::Print()<<" Final Knapsack time: " << amrex::second() - time_start << std::endl<<std::endl;
    output_distribution_data(ba, k_dmap, scaled_wgts, "LBC_knapsack.txt", "Knapsack");

    time_start = amrex::second();
    std::vector<int> s_dmap = SFCProcessorMapDoIt(ba, scaled_wgts, nranks, &s_eff, node_size, true, false, bytes);
    amrex::Print()<<" Final SFC time: " << amrex::second() - time_start << std::endl<<std::endl;
    output_distribution_data(ba, s_dmap, scaled_wgts, "LBC_sfc.txt", "SFC");

    std::set<int> unique_ranks(s_dmap.begin(), s_dmap.end());
    amrex::Print() << "Unique ranks in SFC mapping: ";
    for (auto r : unique_ranks) amrex::Print() << r << " ";
    amrex::Print() << "\n";

    time_start = amrex::second();
    std::vector<int> vec=painterPartition(ba,scaled_wgts,nranks);
    amrex::Print()<<" Final SFC+Painter time: " << amrex::second() - time_start << std::endl<<std::endl;

    time_start = amrex::second();
    std::vector<int> sfc_knapsack_dmap = SFCProcessorMapDoItCombined(ba, scaled_wgts, nnodes, ranks_per_node, &sfc_eff, &knapsack_eff, true, false, bytes);
    amrex::Print()<<" Final SFC+Knapsack_Combined time: " << amrex::second() - time_start << std::endl<<std::endl;
    
    time_start = amrex::second();
    std::vector<int> painter_knapsack_dmap = SFCProcessorMapDoItCombinedPainter(ba, scaled_wgts, nnodes, ranks_per_node, &sfc_eff, &knapsack_eff, true, false, bytes);
    amrex::Print()<<" Final painter+Knapsack_Combined time: " << amrex::second() - time_start << std::endl;

    time_start = amrex::second();
    std::vector<int> hilbertsfc_dmap = HilbertProcessorMapDoIt(ba, scaled_wgts, nranks, &hilbertsfc_eff, ranks_per_node, true, false, bytes);
    amrex::Print()<<" Final Hilbert SFC time: " << amrex::second() - time_start << std::endl;
    output_distribution_data(ba, hilbertsfc_dmap, scaled_wgts, "LBC_hilbert.txt", "Hilbert SFC");

    amrex::Print() << "\n=== End of Run " << r + 1 << " ===\n";
    amrex::Print() << "======================================\n\n";

}

    // Print SFC and Knapsack efficiencies
    // amrex::Print() << "SFC Efficiency: " << sfc_eff << std::endl;
    // amrex::Print() << "Knapsack Efficiency: " << knapsack_eff << std::endl;
}

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
// // #include "TopologyAware.H"
// #include "HilbertSFC.H"
// #include "HilbertSFC_Knapsack.H"
// #include "heterogeneous/HeterogeneousLB.H"

// // #if defined(AMREX_USE_MPI) || defined(AMREX_USE_GPU)
// // #error This is a serial test only.
// // #endif

// using namespace amrex;

// // Forward declarations
// void main_main();
// extern void example_main();  // Declare as external

// int main(int argc, char* argv[]) {
//     amrex::Initialize(argc, argv);

//     // Check if we should run the example or the main test
//     ParmParse pp;
//     std::string test_type = "main";
//     pp.query("test_type", test_type);

//     if (test_type == "example") {
//         example_main();
//     } else {
//         main_main();
//     }

//     amrex::Finalize();
//     return 0;
// }

// void main_main() {
//     BL_PROFILE("main");

//     // Read parameters from inputs file
//     ParmParse pp;
    
//     // Domain size
//     Vector<int> n_cell(AMREX_SPACEDIM);
//     if (!pp.queryarr("domain", n_cell)) {
//         // Default domain size if not specified
//         n_cell[0] = n_cell[1] = n_cell[2] = 64;
//     }
    
//     // Grid size
//     Vector<int> max_grid_size(AMREX_SPACEDIM);
//     if (!pp.queryarr("max_grid_size", max_grid_size)) {
//         // Default grid size if not specified
//         max_grid_size[0] = max_grid_size[1] = max_grid_size[2] = 32;
//     }
    
//     // Create domain box
//     Box domain(IntVect(0), IntVect(n_cell[0]-1, n_cell[1]-1, n_cell[2]-1));
//     BoxArray ba(domain);
//     ba.maxSize(max_grid_size[0]);
    
//     // Read node configuration
//     int nnodes = 6;  // Default number of nodes
//     pp.query("nnodes", nnodes);
    
//     Vector<std::string> node_types;
//     Vector<Real> perf_factors;
//     Vector<Real> memory_caps;
    
//     if (!pp.queryarr("node_types", node_types) ||
//         !pp.queryarr("perf_factors", perf_factors) ||
//         !pp.queryarr("memory_caps", memory_caps)) {
//         // Default node configuration if not specified
//         node_types = {"CPU_AMD_EPYC_7763", "CPU_AMD_EPYC_7763",
//                      "CPU_Intel_8380", "CPU_Intel_8380",
//                      "CPU_AMD_EPYC_7542", "CPU_AMD_EPYC_7542"};
//         perf_factors = {1.0, 1.0, 0.9, 0.9, 0.8, 0.8};
//         memory_caps = {256.0, 256.0, 128.0, 128.0, 128.0, 128.0};
//     }
    
//     // Create node info
//     Vector<ComputeNodeInfo> nodes;
//     for (int i = 0; i < nnodes && i < node_types.size(); ++i) {
//         nodes.push_back({i, node_types[i], perf_factors[i], memory_caps[i]});
//     }
    
//     // Create and initialize load balancer
//     HeterogeneousLB lb;
//     lb.InitializeNodes(nodes);
    
//     // Create initial distribution mapping
//     DistributionMapping dm(ba, nodes.size());
    
//     // Create weights with spatial variation
//     MultiFab weights(ba, dm, 1, 0);
//     for (MFIter mfi(weights); mfi.isValid(); ++mfi) {
//         const Box& box = mfi.validbox();
//         IntVect center = box.smallEnd() + box.size() / 2;
//         Real distance = std::sqrt(
//             std::pow(center[0] - n_cell[0]/2.0, 2) +
//             std::pow(center[1] - n_cell[1]/2.0, 2) +
//             std::pow(center[2] - n_cell[2]/2.0, 2)
//         );
//         weights[mfi].setVal(box.numPts() * (1.0 + std::exp(-distance/(n_cell[0]/4.0))), box, 0, 1);
//     }
    
//     // Test different load balancing strategies
//     amrex::Print() << "\nTesting different load balancing strategies:\n";
    
//     std::string strategy = "knapsack";  // Default strategy
//     pp.query("strategy", strategy);
    
//     auto new_dm = lb.BalanceLoad(ba, weights, strategy);
//     amrex::Print() << "\nAfter " << strategy << " balancing:\n";
//     lb.PrintStats();
// }
