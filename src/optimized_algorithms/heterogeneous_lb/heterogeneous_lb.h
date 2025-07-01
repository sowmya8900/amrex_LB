#ifndef HETEROGENEOUS_LB_H
#define HETEROGENEOUS_LB_H

#include <vector>
#include <string>
#include <AMReX_INT.H>

struct Node {
    int id;
    std::string name;
    double performance_factor;
    double current_load;
    
    Node(int _id, const std::string& _name, double _perf, double _load = 0.0)
        : id(_id), name(_name), performance_factor(_perf), current_load(_load) {}
};

struct Task {
    int id;
    double base_time;
    
    Task(int _id, double _time) : id(_id), base_time(_time) {}
};

struct Migration {
    int source_node;
    int target_node;
    int task_id;
    double amount;
    
    Migration(int _source, int _target, int _task_id, double _amount)
        : source_node(_source), target_node(_target), task_id(_task_id), amount(_amount) {}
};

// Simplified system metrics structure
struct SystemMetrics {
    double makespan;
    double load_balance_index;
    double system_utilization;
    double efficiency;
    double heterogeneity_factor;
    std::vector<double> node_loads;
};

// Core functions
std::vector<std::vector<double>> compute_rij(const std::vector<Node>& nodes);

// Three main assignment approaches
std::vector<int> homogeneous_knapsack_assignment(const std::vector<Task>& tasks, const std::vector<Node>& nodes);
std::vector<int> knapsack_without_rij_assignment(const std::vector<Task>& tasks, const std::vector<Node>& nodes);
std::vector<int> knapsack_with_rij_assignment(const std::vector<Task>& tasks, const std::vector<Node>& nodes, const std::vector<std::vector<double>>& rij);

// System metrics calculation
SystemMetrics calculate_system_metrics(const std::vector<Node>& nodes, double total_work, const std::string& approach_name);

// Main functions
void run_performance_comparison();
void demonstrate_key_differences();

#endif