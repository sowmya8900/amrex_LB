#include <AMReX.H>
#include <AMReX_Print.H>

// Forward declaration of test function
void test_heterogeneous_lb();

int main(int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    
    amrex::Print() << "\nStarting Heterogeneous Load Balancer Tests...\n";
    
    test_heterogeneous_lb();
    
    amrex::Print() << "\nTests completed.\n";
    
    amrex::Finalize();
    return 0;
} 