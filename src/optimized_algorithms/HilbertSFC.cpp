#include <AMReX_Extension.H>
#include <AMReX_IntVect.H>
#include <AMReX_Print.H>
#include <AMReX_ParallelContext.H>
#include <AMReX.H>
#include <array>
#include <algorithm>
#include <fstream>

#include "HilbertSFC.H"
#include "Knapsack.H"
#include "LeastUsed.H"

namespace {
    // Helper function to get sign of a number
    inline int sgn(int x) {
        return (x < 0) ? -1 : ((x > 0) ? 1 : 0);
    }

    // Core recursive function for 3D Hilbert curve generation
    uint64_t generate3d(int x, int y, int z,
                       int ax, int ay, int az,
                       int bx, int by, int bz,
                       int cx, int cy, int cz,
                       int bits) {
        int w = std::abs(ax + ay + az);
        int h = std::abs(bx + by + bz);
        int d = std::abs(cx + cy + cz);

        // Unit direction vectors
        int dax = sgn(ax), day = sgn(ay), daz = sgn(az);  // major direction ("right")
        int dbx = sgn(bx), dby = sgn(by), dbz = sgn(bz);  // ortho direction ("forward")
        int dcx = sgn(cx), dcy = sgn(cy), dcz = sgn(cz);  // ortho direction ("up")

        // Base cases for trivial dimensions
        if (h == 1 && d == 1) {
            return x | (y << bits) | (z << (2 * bits));
        }
        if (w == 1 && d == 1) {
            return x | (y << bits) | (z << (2 * bits));
        }
        if (w == 1 && h == 1) {
            return x | (y << bits) | (z << (2 * bits));
        }

        // Split dimensions
        int ax2 = ax/2, ay2 = ay/2, az2 = az/2;
        int bx2 = bx/2, by2 = by/2, bz2 = bz/2;
        int cx2 = cx/2, cy2 = cy/2, cz2 = cz/2;

        int w2 = std::abs(ax2 + ay2 + az2);
        int h2 = std::abs(bx2 + by2 + bz2);
        int d2 = std::abs(cx2 + cy2 + cz2);

        // Prefer even steps
        if ((w2 % 2) && (w > 2)) {
            ax2 += dax; ay2 += day; az2 += daz;
        }
        if ((h2 % 2) && (h > 2)) {
            bx2 += dbx; by2 += dby; bz2 += dbz;
        }
        if ((d2 % 2) && (d > 2)) {
            cx2 += dcx; cy2 += dcy; cz2 += dcz;
        }

        uint64_t index;
        
        // Wide case: split in w only
        if ((2*w > 3*h) && (2*w > 3*d)) {
            if (x < ax2 + ay2 + az2) {
                index = generate3d(x, y, z,
                                 ax2, ay2, az2,
                                 bx, by, bz,
                                 cx, cy, cz,
                                 bits);
            } else {
                index = generate3d(x-(ax2+ay2+az2), y-(ay2+ay2+az2), z-(az2+az2+az2),
                                 ax-ax2, ay-ay2, az-az2,
                                 bx, by, bz,
                                 cx, cy, cz,
                                 bits);
            }
        }
        // Do not split in d
        else if (3*h > 4*d) {
            if (y < by2 + bz2) {
                index = generate3d(x, y, z,
                                 bx2, by2, bz2,
                                 cx, cy, cz,
                                 ax2, ay2, az2,
                                 bits);
            } else if (x < ax + bx2) {
                index = generate3d(x-bx2, y-by2, z-bz2,
                                 ax, ay, az,
                                 bx-bx2, by-by2, bz-bz2,
                                 cx, cy, cz,
                                 bits);
            } else {
                index = generate3d(x-(ax-dax)-(bx2-dbx),
                                 y-(ay-day)-(by2-dby),
                                 z-(az-daz)-(bz2-dbz),
                                 -bx2, -by2, -bz2,
                                 cx, cy, cz,
                                 -(ax-ax2), -(ay-ay2), -(az-az2),
                                 bits);
            }
        }
        // Do not split in h
        else if (3*d > 4*h) {
            if (z < cz2) {
                index = generate3d(x, y, z,
                                 cx2, cy2, cz2,
                                 ax2, ay2, az2,
                                 bx, by, bz,
                                 bits);
            } else if (x < ax + cx2) {
                index = generate3d(x-cx2, y-cy2, z-cz2,
                                 ax, ay, az,
                                 bx, by, bz,
                                 cx-cx2, cy-cy2, cz-cz2,
                                 bits);
            } else {
                index = generate3d(x-(ax-dax)-(cx2-dcx),
                                 y-(ay-day)-(cy2-dcy),
                                 z-(az-daz)-(cz2-dcz),
                                 -cx2, -cy2, -cz2,
                                 -(ax-ax2), -(ay-ay2), -(az-az2),
                                 bx, by, bz,
                                 bits);
            }
        }
        // Regular case: split in all dimensions
        else {
            if (z < cz2 && y < by2) {
                index = generate3d(x, y, z,
                                 bx2, by2, bz2,
                                 cx2, cy2, cz2,
                                 ax2, ay2, az2,
                                 bits);
            } else if (z < cz + by2) {
                index = generate3d(x-bx2, y-by2, z-bz2,
                                 cx, cy, cz,
                                 ax2, ay2, az2,
                                 bx-bx2, by-by2, bz-bz2,
                                 bits);
            } else if (x < ax + bx2 + cx) {
                index = generate3d(x-(bx2-dbx)-(cx-dcx),
                                 y-(by2-dby)-(cy-dcy),
                                 z-(bz2-dbz)-(cz-dcz),
                                 ax, ay, az,
                                 -bx2, -by2, -bz2,
                                 -(cx-cx2), -(cy-cy2), -(cz-cz2),
                                 bits);
            } else if (y < by + cy) {
                index = generate3d(x-(ax-dax)-bx2-(cx-dcx),
                                 y-(ay-day)-by2-(cy-dcy),
                                 z-(az-daz)-bz2-(cz-dcz),
                                 -cx, -cy, -cz,
                                 -(ax-ax2), -(ay-ay2), -(az-az2),
                                 bx-bx2, by-by2, bz-bz2,
                                 bits);
            } else {
                index = generate3d(x-(ax-dax)-(bx2-dbx),
                                 y-(ay-day)-(by2-dby),
                                 z-(az-daz)-(bz2-dbz),
                                 -bx2, -by2, -bz2,
                                 cx2, cy2, cz2,
                                 -(ax-ax2), -(ay-ay2), -(az-az2),
                                 bits);
            }
        }

        return index;
    }
} // anonymous namespace

uint64_t coords_to_hilbert(uint32_t x, uint32_t y, uint32_t z, int order) {
    // Calculate the size of the grid (2^order)
    uint32_t size = 1 << order;
    
    // Ensure coordinates are within bounds
    x = std::min(x, size - 1);
    y = std::min(y, size - 1);
    z = std::min(z, size - 1);

    // Call the Princeton LIPS implementation
    return generate3d(x, y, z,
                     size, 0, 0,  // width direction
                     0, size, 0,  // height direction
                     0, 0, size,  // depth direction
                     order);      // bits per dimension
}

void HilbertDistribute(const std::vector<HilbertSFCToken>& tokens,
                      const std::vector<amrex::Long>& wgts,
                      int nprocs,
                      amrex::Real volpercpu,
                      std::vector<std::vector<int>>& v,
                      bool flag_verbose_mapper) {
    BL_PROFILE("HilbertDistribute()");

    BL_ASSERT(static_cast<int>(v.size()) == nprocs);

    int K = 0;
    amrex::Real totalvol = 0;

    for (int i = 0; i < nprocs; ++i) {
        int cnt = 0;
        amrex::Real vol = 0;

        for (int TSZ = static_cast<int>(tokens.size());
             K < TSZ && (i == (nprocs - 1) || (vol < volpercpu));
             ++K) {
            vol += wgts[tokens[K].m_box];
            ++cnt;
            v[i].push_back(tokens[K].m_box);
        }

        totalvol += vol;

        if ((totalvol / (i + 1)) > volpercpu &&  // Too much for this bin
            cnt > 1 &&                          // More than one box in this bin
            i < nprocs - 1) {                    // Not the last bin
            --K;
            v[i].pop_back();
            totalvol -= wgts[tokens[K].m_box];
        }
    }

#ifdef AMREX_DEBUG
    int cnt = 0;
    for (int i = 0; i < nprocs; ++i) {
        cnt += v[i].size();
    }
    BL_ASSERT(cnt == static_cast<int>(tokens.size()));
#endif
}

std::vector<int>
HilbertProcessorMapDoIt(const amrex::BoxArray& boxes,
                       const std::vector<amrex::Long>& wgts,
                       int nprocs,
                       amrex::Real* eff,
                       int node_size,
                       bool flag_verbose_mapper,
                       bool sort,
                       const std::vector<amrex::Long>& bytes) {
    if (flag_verbose_mapper) {
        amrex::Print() << "DM: HilbertProcessorMapDoIt called..." << std::endl;
    }

    BL_PROFILE("HilbertProcessorMapDoIt()");

    int nteams = nprocs;
    int nworkers = 1;

    if (node_size > 0) {
        nteams = nprocs / node_size;
        nworkers = node_size;
        if (nworkers * nteams != nprocs) {
            nteams = nprocs;
            nworkers = 1;
        }
    }

    if (flag_verbose_mapper) {
        amrex::Print() << "  (nprocs, nteams, nworkers) = ("
                       << nprocs << ", " << nteams << ", " << nworkers << ")\n";
    }

    const int N = boxes.size();
    std::vector<HilbertSFCToken> tokens;
    tokens.reserve(N);
    for (int i = 0; i < N; ++i) {
        const amrex::Box& bx = boxes[i];
        tokens.push_back(makeHilbertSFCToken(i, bx.smallEnd()));
    }

    // Sort tokens along Hilbert curve
    std::sort(tokens.begin(), tokens.end(), HilbertSFCToken::Compare());

    // Split them up as equitably as possible per team
    amrex::Real volperteam = 0;
    for (amrex::Long wt : wgts) {
        volperteam += wt;
    }
    volperteam /= nteams;

    std::vector<std::vector<int>> vec(nteams);
    HilbertDistribute(tokens, wgts, nteams, volperteam, vec);

    tokens.clear();

    std::vector<LIpair> LIpairV;
    LIpairV.reserve(nteams);

    for (int i = 0; i < nteams; ++i) {
        amrex::Long wgt = 0;
        const std::vector<int>& vi = vec[i];
        for (int j = 0, M = vi.size(); j < M; ++j) {
            wgt += wgts[vi[j]];
        }
        LIpairV.push_back(LIpair(wgt, i));
    }

    if (sort) Sort(LIpairV, true);

    if (flag_verbose_mapper) {
        for (const auto &p : LIpairV) {
            amrex::Print() << "  Bucket " << p.second << " contains " << p.first << std::endl;
        }
    }

    amrex::Vector<int> ord;
    amrex::Vector<amrex::Vector<int>> wrkerord;

    if (nteams == nprocs) {
        if (sort) {
            LeastUsedCPUs(nprocs, bytes, ord, flag_verbose_mapper);
        } else {
            ord.resize(nprocs);
            std::iota(ord.begin(), ord.end(), 0);
        }
    } else {
        if (sort) {
            // LeastUsedTeams(ord, wrkerord, nteams, nworkers);
        } else {
            ord.resize(nteams);
            std::iota(ord.begin(), ord.end(), 0);
            wrkerord.resize(nteams);
            for (auto& v : wrkerord) {
                v.resize(nworkers);
                std::iota(v.begin(), v.end(), 0);
            }
        }
    }

    std::vector<int> result(wgts.size());

    for (int i = 0; i < nteams; ++i) {
        const int tid = ord[i];
        const int ivec = LIpairV[i].second;
        const std::vector<int>& vi = vec[ivec];
        const int Nbx = vi.size();

        if (nteams == nprocs) {
            for (int j = 0; j < Nbx; ++j) {
                result[vi[j]] = amrex::ParallelContext::local_to_global_rank(tid);
            }
        } else {
            std::vector<amrex::Long> local_wgts;
            for (int j = 0; j < Nbx; ++j) {
                local_wgts.push_back(wgts[vi[j]]);
            }

            std::vector<std::vector<int>> kpres;
            amrex::Real kpeff;
            knapsack(local_wgts, nworkers, kpres, kpeff, true, N);

            std::vector<LIpair> ww;
            for (int w = 0; w < nworkers; ++w) {
                amrex::Long wgt = 0;
                for (std::vector<int>::const_iterator it = kpres[w].begin();
                     it != kpres[w].end(); ++it) {
                    wgt += local_wgts[*it];
                }
                ww.push_back(LIpair(wgt, w));
            }
            Sort(ww, true);

            const amrex::Vector<int>& sorted_workers = wrkerord[i];
            const int leadrank = tid * nworkers;

            for (int w = 0; w < nworkers; ++w) {
                const int cpu = leadrank + sorted_workers[w];
                int ikp = ww[w].second;
                const std::vector<int>& js = kpres[ikp];
                for (std::vector<int>::const_iterator it = js.begin(); it != js.end(); ++it) {
                    result[vi[*it]] = cpu;
                }
            }
        }
    }

    if (eff || flag_verbose_mapper) {
        amrex::Real sum_wgt = 0, max_wgt = 0;
        for (int i = 0; i < nteams; ++i) {
            const amrex::Long W = LIpairV[i].first;
            if (W > max_wgt) max_wgt = W;
            sum_wgt += W;
        }
        amrex::Real efficiency = (sum_wgt / (nteams * max_wgt));
        if (eff) *eff = efficiency;

        if (flag_verbose_mapper) {
            amrex::Print() << "Hilbert SFC efficiency: " << efficiency << '\n';
        }
    }

    return result;
}