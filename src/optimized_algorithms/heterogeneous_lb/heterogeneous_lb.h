// #ifndef HETEROGENEOUS_LB_H
// #define HETEROGENEOUS_LB_H

// #include <vector>
// #include <string>

// struct Node {
//     int id;
//     std::string name;
//     double performance_factor; // Higher is better (e.g., 1.0 = baseline, 2.0 = twice as fast)
//     double current_time;       // Accumulated time of assigned tasks
// };

// struct Task {
//     int id;
//     double base_time; // Time required on baseline node
// };

// class HeterogeneousLB {
// public:
//     HeterogeneousLB(const std::vector<Node>& nodes);
//     int assign_task(const Task& task); // Returns node id
//     std::vector<double> projected_times(const Task& task) const;
//     const std::vector<Node>& get_nodes() const;
//     void reset();
// private:
//     std::vector<Node> nodes_;
// };

// #endif // HETEROGENEOUS_LB_H 

#ifndef HETEROGENEOUS_LB_H
#define HETEROGENEOUS_LB_H

#include <vector>
#include <string>

struct Node {
    int id;
    std::string name;
    double performance_factor;
    double current_load;  // Add this field
    
    // Constructor
    Node(int _id, const std::string& _name, double _perf, double _load = 0.0) 
        : id(_id), name(_name), performance_factor(_perf), current_load(_load) {}
};

struct Task {
    int id;
    double base_time;
    
    // Constructor
    Task(int _id, double _time) : id(_id), base_time(_time) {}
};

// Function declarations
std::vector<std::vector<double>> compute_rij(const std::vector<Node>& nodes);
std::vector<int> knapsack_with_rij_assignment(const std::vector<Task>& tasks, 
                                             const std::vector<Node>& nodes,
                                             const std::vector<std::vector<double>>& rij);
std::vector<int> greedy_with_rij_assignment(const std::vector<Task>& tasks, 
                                           const std::vector<Node>& nodes,
                                           const std::vector<std::vector<double>>& rij);

#endif