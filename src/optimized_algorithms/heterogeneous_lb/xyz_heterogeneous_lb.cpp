#include "heterogeneous_lb.h"
#include <limits>
#include <vector>
#include <iostream>
#include "../Knapsack.H"

HeterogeneousLB::HeterogeneousLB(const std::vector<Node>& nodes) : nodes_(nodes) {}

int HeterogeneousLB::assign_task(const Task& task) {
    int best_node = -1;
    double min_time = std::numeric_limits<double>::max();
    for (size_t i = 0; i < nodes_.size(); ++i) {
        double proj_time = nodes_[i].current_time + task.base_time / nodes_[i].performance_factor;
        if (proj_time < min_time) {
            min_time = proj_time;
            best_node = static_cast<int>(i);
        }
    }
    if (best_node >= 0) {
        nodes_[best_node].current_time += task.base_time / nodes_[best_node].performance_factor;
        return nodes_[best_node].id;
    }
    return -1;
}

std::vector<double> HeterogeneousLB::projected_times(const Task& task) const {
    std::vector<double> times;
    for (const auto& node : nodes_) {
        times.push_back(node.current_time + task.base_time / node.performance_factor);
    }
    return times;
}

const std::vector<Node>& HeterogeneousLB::get_nodes() const {
    return nodes_;
}

void HeterogeneousLB::reset() {
    for (auto& node : nodes_) node.current_time = 0.0;
}

// Knapsack wrapper for this project
std::vector<int> knapsack_assignment(const std::vector<Task>& tasks, int n_nodes) {
    std::vector<long> wgts;
    for (const auto& t : tasks) wgts.push_back(static_cast<long>(t.base_time));
    double eff = 0.0;
    return KnapSackDoIt(wgts, n_nodes, eff, true, std::numeric_limits<int>::max(), false, false, std::vector<amrex::Long>());
} 