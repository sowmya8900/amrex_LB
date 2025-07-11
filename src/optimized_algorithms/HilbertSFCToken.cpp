#include "HilbertSFCToken.H"
#include "HilbertSFC.H"
#include <AMReX_IntVect.H>
#include <algorithm>

bool HilbertSFCToken::Compare::operator()(const HilbertSFCToken& lhs,
                                        const HilbertSFCToken& rhs) const {
    return lhs.m_hilbert_key < rhs.m_hilbert_key;
}

HilbertSFCToken makeHilbertSFCToken(int box_index, amrex::IntVect const& iv) {
    HilbertSFCToken token;
    token.m_box = box_index;

    // Use order 10 for reasonable performance (1024^3 grid)
    // constexpr int order = 10;
    constexpr int order = 8;  // 256^3 grid for better performance in large simulations
    constexpr int imin = -(1 << (order-1));  // -512
    constexpr int imax = (1 << (order-1)) - 1;  // 511
    
    // Clamp coordinates to valid range and shift to positive
    uint32_t x = std::max(imin, std::min(imax, iv[0])) - imin;
    uint32_t y = std::max(imin, std::min(imax, iv[1])) - imin;
    uint32_t z = std::max(imin, std::min(imax, iv[2])) - imin;

    token.m_hilbert_key = coords_to_hilbert(x, y, z, order);
    return token;
}