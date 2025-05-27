#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include "HeterogeneousLB.H"

using namespace amrex;

// Forward declarations
void example_main();
void test_heterogeneous_lb();

int main(int argc, char* argv[])
{
    amrex::Initialize(argc, argv);

    // Parse test type from command line or inputs file
    ParmParse pp;
    std::string test_type = "example";
    pp.query("test_type", test_type);
    
    if (test_type == "test") {
        amrex::Print() << "\nRunning heterogeneous load balancer tests...\n";
        test_heterogeneous_lb();
    } else {
        example_main();
    }

    amrex::Finalize();
    return 0;
} 