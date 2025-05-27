#include "HilbertSFCToken.H"
#include "HilbertSFC.H"
#include <AMReX_IntVect.H>

bool HilbertSFCToken::Compare::operator()(const HilbertSFCToken& lhs,
                                        const HilbertSFCToken& rhs) const {
    return lhs.m_hilbert_key < rhs.m_hilbert_key;
}

HilbertSFCToken makeHilbertSFCToken(int box_index, amrex::IntVect const& iv) {
    HilbertSFCToken token;
    token.m_box = box_index;

    constexpr int imin = -(1 << 20);
    AMREX_ASSERT_WITH_MESSAGE(iv[0] >= imin && iv[0] < -imin &&
                            iv[1] >= imin && iv[1] < -imin &&
                            iv[2] >= imin && iv[2] < -imin,
                          "HilbertSFCToken: index out of range");

    uint32_t x = iv[0] - imin;
    uint32_t y = iv[1] - imin;
    uint32_t z = iv[2] - imin;

    token.m_hilbert_key = coords_to_hilbert(x, y, z, 21);
    return token;
}