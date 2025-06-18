#include "heterogeneous_lb.h"
#include <iostream>
#include <vector>
#include <random>
#include <fstream>
#include <iomanip>
#include "../Knapsack.H"
#include <AMReX_INT.H>
#include <limits>
#include <algorithm>

// Compute rij matrix: rij[i][j] = perf[j] / perf[i]
std::vector<std::vector<double>> compute_rij(const std::vector<Node>& nodes) {
    size_t n = nodes.size();
    std::vector<std::vector<double>> rij(n, std::vector<double>(n, 1.0));
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < n; ++j)
            rij[i][j] = nodes[j].performance_factor / nodes[i].performance_factor;
    return rij;
}

// Enhanced knapsack using rij matrix for weight calculation
std::vector<int> knapsack_with_rij_assignment(const std::vector<Task>& tasks, 
                                             const std::vector<Node>& nodes,
                                             const std::vector<std::vector<double>>& rij) {
    std::vector<amrex::Long> rij_adjusted_weights;
    
    for (const auto& task : tasks) {
        double min_rij = std::numeric_limits<double>::max();
        for (size_t i = 0; i < nodes.size(); ++i) {
            if (rij[0][i] < min_rij) {
                min_rij = rij[0][i];
            }
        }
        double worst_case_time = task.base_time / min_rij;
        rij_adjusted_weights.push_back(static_cast<amrex::Long>(worst_case_time * 100));
    }
    
    double eff = 0.0;
    return KnapSackDoIt(rij_adjusted_weights, nodes.size(), eff, true, 
                       std::numeric_limits<int>::max(), false, false, 
                       std::vector<amrex::Long>());
}

// Enhanced greedy assignment explicitly using rij matrix
std::vector<int> greedy_with_rij_assignment(const std::vector<Task>& tasks, 
                                           const std::vector<Node>& nodes,
                                           const std::vector<std::vector<double>>& rij) {
    std::vector<int> assignments;
    std::vector<double> node_times(nodes.size(), 0.0);
    
    for (size_t i = 0; i < nodes.size(); ++i) {
        node_times[i] = nodes[i].current_load;
    }
    
    for (const auto& task : tasks) {
        int best_node = -1;
        double min_completion_time = std::numeric_limits<double>::max();
        
        for (size_t i = 0; i < nodes.size(); ++i) {
            double exec_time = task.base_time / rij[0][i];
            double completion_time = node_times[i] + exec_time;
            
            if (completion_time < min_completion_time) {
                min_completion_time = completion_time;
                best_node = i;
            }
        }
        
        assignments.push_back(best_node);
        node_times[best_node] += task.base_time / rij[0][best_node];
    }
    
    return assignments;
}

// Original dynamic knapsack (for comparison)
std::vector<int> dynamic_knapsack_assignment(const std::vector<Task>& tasks, const std::vector<Node>& nodes) {
    std::vector<amrex::Long> adj_wgts;
    for (const auto& t : tasks) {
        double min_time = std::numeric_limits<double>::max();
        for (const auto& n : nodes) {
            double t_adj = t.base_time / n.performance_factor;
            if (t_adj < min_time) min_time = t_adj;
        }
        adj_wgts.push_back(static_cast<amrex::Long>(min_time * 100));
    }
    double eff = 0.0;
    return KnapSackDoIt(adj_wgts, nodes.size(), eff, true, std::numeric_limits<int>::max(), false, false, std::vector<amrex::Long>());
}

// Function 1: Single new task analysis
void analyze_single_new_task() {
    std::cout << "\n=== SINGLE NEW TASK ANALYSIS ===\n";
    
    std::vector<Node> nodes = {
        Node(0, "cpu0", 1.0, 0.0),
        Node(1, "cpu1", 0.8, 0.0),
        Node(2, "cpu2", 1.2, 0.0),
        Node(3, "cpu3", 1.3, 0.0),
        Node(4, "cpu4", 0.9, 0.0),
        Node(5, "cpu5", 1.2, 50.0)  // Pre-loaded
    };
    
    auto rij = compute_rij(nodes);
    
    // Generate tasks
    std::vector<Task> tasks;
    std::mt19937 gen(std::random_device{}());
    std::normal_distribution<> d(10.0, 2.0);
    for (int i = 0; i < 30; ++i) {
        double t = std::max(1.0, d(gen));
        tasks.push_back(Task(i, t));
    }

    // Run strategies
    auto greedy_rij_assign = greedy_with_rij_assignment(tasks, nodes, rij);
    auto knap_rij_assign = knapsack_with_rij_assignment(tasks, nodes, rij);
    
    // Calculate loads
    std::vector<double> node_times_greedy_rij(nodes.size(), 0.0);
    std::vector<double> node_times_knap_rij(nodes.size(), 0.0);
    
    for (size_t i = 0; i < nodes.size(); ++i) {
        node_times_greedy_rij[i] = nodes[i].current_load;
        node_times_knap_rij[i] = nodes[i].current_load;
    }
    
    for (size_t i = 0; i < tasks.size(); ++i) {
        int greedy_node = greedy_rij_assign[i];
        int knap_node = knap_rij_assign[i];
        node_times_greedy_rij[greedy_node] += tasks[i].base_time / rij[0][greedy_node];
        node_times_knap_rij[knap_node] += tasks[i].base_time / rij[0][knap_node];
    }

    // Generate new task for analysis
    std::mt19937 gen_plot(std::random_device{}());
    std::normal_distribution<> d_plot(10.0, 2.0);
    double t = std::max(1.0, d_plot(gen_plot));
    Task new_task(static_cast<int>(tasks.size()), t);
    
    // Output single task analysis
    std::ofstream fplot("results/heterogeneous_lb_rij_plot.csv");
    fplot << "node_id,perf_factor,greedy_rij_current,greedy_rij_added,greedy_rij_proj,knap_rij_current,knap_rij_added,knap_rij_proj\n";
    for (size_t i = 0; i < nodes.size(); ++i) {
        double greedy_curr = node_times_greedy_rij[i];
        double greedy_add = new_task.base_time / rij[0][i];
        double greedy_proj = greedy_curr + greedy_add;
        
        double knap_curr = node_times_knap_rij[i];
        double knap_add = new_task.base_time / rij[0][i];
        double knap_proj = knap_curr + knap_add;
        
        fplot << i << "," << nodes[i].performance_factor << ","
              << std::fixed << std::setprecision(3)
              << greedy_curr << "," << greedy_add << "," << greedy_proj << ","
              << knap_curr << "," << knap_add << "," << knap_proj << "\n";
    }
    fplot.close();
    
    std::cout << "Single task analysis saved to 'results/heterogeneous_lb_rij_plot.csv'\n";
    std::cout << "New task base_time: " << new_task.base_time << "\n";
}

// Function 2: Multiple new tasks analysis
void analyze_multiple_new_tasks() {
    std::cout << "\n=== MULTIPLE NEW TASKS ANALYSIS ===\n";
    
    std::vector<Node> nodes = {
        Node(0, "cpu0", 1.0, 0.0),
        Node(1, "cpu1", 0.8, 0.0),
        Node(2, "cpu2", 1.2, 0.0),
        Node(3, "cpu3", 1.1, 0.0),
        Node(4, "cpu4", 0.9, 0.0),
        Node(5, "cpu5", 1.3, 0.0)
    };
    
    auto rij = compute_rij(nodes);
    
    // Generate base tasks
    std::vector<Task> tasks;
    std::mt19937 gen(42); // Fixed seed for reproducibility
    std::normal_distribution<> d(10.0, 2.0);
    for (int i = 0; i < 30; ++i) {
        double t = std::max(1.0, d(gen));
        tasks.push_back(Task(i, t));
    }

    // Run strategies
    auto greedy_rij_assign = greedy_with_rij_assignment(tasks, nodes, rij);
    auto knap_rij_assign = knapsack_with_rij_assignment(tasks, nodes, rij);
    
    // Calculate loads
    std::vector<double> node_times_greedy_rij(nodes.size(), 0.0);
    std::vector<double> node_times_knap_rij(nodes.size(), 0.0);
    
    for (size_t i = 0; i < tasks.size(); ++i) {
        int greedy_node = greedy_rij_assign[i];
        int knap_node = knap_rij_assign[i];
        node_times_greedy_rij[greedy_node] += tasks[i].base_time / rij[0][greedy_node];
        node_times_knap_rij[knap_node] += tasks[i].base_time / rij[0][knap_node];
    }

    // Generate 5 different new tasks
    std::mt19937 gen_new_tasks(std::random_device{}());
    std::normal_distribution<> d_new_tasks(10.0, 2.0);
    
    std::vector<Task> new_tasks;
    for (int i = 0; i < 5; ++i) {
        double t = std::max(1.0, d_new_tasks(gen_new_tasks));
        new_tasks.push_back(Task(static_cast<int>(tasks.size()) + i, t));
    }

    // Output multiple tasks analysis
    std::ofstream fplot("results/heterogeneous_lb_multiple_new_tasks.csv");
    fplot << "new_task_id,task_base_time,node_id,perf_factor,greedy_rij_current,greedy_rij_added,greedy_rij_proj,knap_rij_current,knap_rij_added,knap_rij_proj\n";
    
    std::cout << "TaskID,BaseTime,Greedy_Best_Node,Greedy_Best_Time,Knap_Best_Node,Knap_Best_Time,Winner\n";
    
    for (size_t task_idx = 0; task_idx < new_tasks.size(); ++task_idx) {
        const Task& new_task = new_tasks[task_idx];
        
        int greedy_best_node = -1, knap_best_node = -1;
        double greedy_best_time = std::numeric_limits<double>::max();
        double knap_best_time = std::numeric_limits<double>::max();
        
        for (size_t i = 0; i < nodes.size(); ++i) {
            double greedy_curr = node_times_greedy_rij[i];
            double greedy_add = new_task.base_time / rij[0][i];
            double greedy_proj = greedy_curr + greedy_add;
            
            double knap_curr = node_times_knap_rij[i];
            double knap_add = new_task.base_time / rij[0][i];
            double knap_proj = knap_curr + knap_add;
            
            if (greedy_proj < greedy_best_time) {
                greedy_best_time = greedy_proj;
                greedy_best_node = i;
            }
            
            if (knap_proj < knap_best_time) {
                knap_best_time = knap_proj;
                knap_best_node = i;
            }
            
            fplot << new_task.id << "," << std::fixed << std::setprecision(3) << new_task.base_time << ","
                  << i << "," << nodes[i].performance_factor << ","
                  << greedy_curr << "," << greedy_add << "," << greedy_proj << ","
                  << knap_curr << "," << knap_add << "," << knap_proj << "\n";
        }
        
        std::string winner = (greedy_best_time < knap_best_time) ? "Greedy" : 
                            (knap_best_time < greedy_best_time) ? "Knapsack" : "Tie";
        
        std::cout << new_task.id << "," << std::fixed << std::setprecision(2) << new_task.base_time << ","
                  << greedy_best_node << "," << greedy_best_time << ","
                  << knap_best_node << "," << knap_best_time << ","
                  << winner << "\n";
    }
    
    fplot.close();
    std::cout << "Multiple tasks analysis saved to 'results/heterogeneous_lb_multiple_new_tasks.csv'\n";
}

// Function 3: Comprehensive scenario analysis
void test_comprehensive_scenarios() {
    std::cout << "\n=== COMPREHENSIVE SCENARIO ANALYSIS ===\n";
    
    std::vector<std::vector<double>> perf_profiles = {
        {1.0, 0.8, 1.2, 1.3, 0.9, 1.2},  // Moderate heterogeneity
        {1.0, 0.5, 1.5, 2.0, 0.7, 1.8},  // High heterogeneity
        {1.0, 1.0, 1.0, 1.0, 1.0, 1.0},  // Homogeneous
        {0.8, 1.0, 1.2, 1.4, 1.6, 1.8},  // Gradual increase
        {2.0, 0.5, 1.8, 0.6, 1.5, 0.7}   // Mixed performance
    };
    
    std::vector<std::vector<double>> preload_scenarios = {
        {0, 0, 0, 0, 0, 0},      // No pre-load
        {0, 0, 0, 0, 0, 50},     // Node 5 overloaded
        {30, 0, 0, 0, 0, 0},     // Node 0 overloaded
        {0, 0, 40, 0, 0, 0},     // Node 2 overloaded
        {10, 10, 10, 0, 0, 0},   // First 3 nodes pre-loaded
        {5, 15, 25, 35, 45, 55}  // Increasing load pattern
    };
    
    std::vector<int> seeds = {42, 123, 456, 789, 999};
    std::vector<double> test_task_sizes = {5.0, 10.0, 15.0, 20.0, 30.0};
    
    std::ofstream results_file("results/comprehensive_scenario_analysis.csv");
    results_file << "Scenario_ID,Perf_Profile,Preload_Profile,Seed,Task_Size,"
                 << "Greedy_Best_Node,Greedy_Best_Time,Knap_Best_Node,Knap_Best_Time,"
                 << "Winner,Time_Diff,Greedy_Makespan,Knap_Makespan\n";
    
    int scenario_id = 0;
    
    for (size_t perf_idx = 0; perf_idx < perf_profiles.size(); ++perf_idx) {
        for (size_t preload_idx = 0; preload_idx < preload_scenarios.size(); ++preload_idx) {
            for (int seed : seeds) {
                std::vector<Node> nodes;
                for (int i = 0; i < 6; ++i) {
                    nodes.push_back(Node(i, "cpu" + std::to_string(i), 
                                        perf_profiles[perf_idx][i], 
                                        preload_scenarios[preload_idx][i]));
                }
                
                auto rij = compute_rij(nodes);
                
                std::vector<Task> tasks;
                std::mt19937 gen(seed);
                std::normal_distribution<> d(10.0, 2.0);
                for (int i = 0; i < 30; ++i) {
                    double t = std::max(1.0, d(gen));
                    tasks.push_back(Task(i, t));
                }
                
                auto greedy_assign = greedy_with_rij_assignment(tasks, nodes, rij);
                auto knap_assign = knapsack_with_rij_assignment(tasks, nodes, rij);
                
                std::vector<double> node_times_greedy(nodes.size(), 0.0);
                std::vector<double> node_times_knap(nodes.size(), 0.0);
                
                for (size_t i = 0; i < nodes.size(); ++i) {
                    node_times_greedy[i] = nodes[i].current_load;
                    node_times_knap[i] = nodes[i].current_load;
                }
                
                for (size_t i = 0; i < tasks.size(); ++i) {
                    int greedy_node = greedy_assign[i];
                    int knap_node = knap_assign[i];
                    node_times_greedy[greedy_node] += tasks[i].base_time / rij[0][greedy_node];
                    node_times_knap[knap_node] += tasks[i].base_time / rij[0][knap_node];
                }
                
                double greedy_makespan = *std::max_element(node_times_greedy.begin(), node_times_greedy.end());
                double knap_makespan = *std::max_element(node_times_knap.begin(), node_times_knap.end());
                
                for (double task_size : test_task_sizes) {
                    int greedy_best = -1, knap_best = -1;
                    double greedy_min = std::numeric_limits<double>::max();
                    double knap_min = std::numeric_limits<double>::max();
                    
                    for (size_t i = 0; i < nodes.size(); ++i) {
                        double greedy_proj = node_times_greedy[i] + task_size / rij[0][i];
                        double knap_proj = node_times_knap[i] + task_size / rij[0][i];
                        
                        if (greedy_proj < greedy_min) {
                            greedy_min = greedy_proj;
                            greedy_best = i;
                        }
                        if (knap_proj < knap_min) {
                            knap_min = knap_proj;
                            knap_best = i;
                        }
                    }
                    
                    std::string winner = (greedy_min < knap_min) ? "Greedy" : "Knapsack";
                    double time_diff = std::abs(greedy_min - knap_min);
                    
                    results_file << scenario_id << "," << perf_idx << "," << preload_idx << "," 
                                << seed << "," << task_size << ","
                                << greedy_best << "," << std::fixed << std::setprecision(3) << greedy_min << ","
                                << knap_best << "," << knap_min << ","
                                << winner << "," << time_diff << ","
                                << greedy_makespan << "," << knap_makespan << "\n";
                    
                    scenario_id++;
                }
            }
        }
    }
    
    results_file.close();
    std::cout << "Comprehensive analysis completed. Results saved to 'results/comprehensive_scenario_analysis.csv'\n";
    std::cout << "Total scenarios tested: " << scenario_id << "\n";
}

// Function 4: Detailed scenario analysis
void analyze_detailed_scenario() {
    std::cout << "\n=== DETAILED SCENARIO ANALYSIS ===\n";
    
    std::vector<Node> nodes = {
        Node(0, "cpu0", 1.0, 0.0),
        Node(1, "cpu1", 0.5, 0.0),
        Node(2, "cpu2", 1.5, 0.0),
        Node(3, "cpu3", 2.0, 0.0),
        Node(4, "cpu4", 0.7, 0.0),
        Node(5, "cpu5", 1.8, 50.0)  // High performance but pre-loaded
    };
    
    auto rij = compute_rij(nodes);
    
    std::vector<Task> tasks;
    std::mt19937 gen(42);
    std::normal_distribution<> d(10.0, 2.0);
    for (int i = 0; i < 30; ++i) {
        double t = std::max(1.0, d(gen));
        tasks.push_back(Task(i, t));
    }
    
    auto greedy_assign = greedy_with_rij_assignment(tasks, nodes, rij);
    auto knap_assign = knapsack_with_rij_assignment(tasks, nodes, rij);
    
    std::vector<double> node_times_greedy(nodes.size(), 0.0);
    std::vector<double> node_times_knap(nodes.size(), 0.0);
    
    for (size_t i = 0; i < nodes.size(); ++i) {
        node_times_greedy[i] = nodes[i].current_load;
        node_times_knap[i] = nodes[i].current_load;
    }
    
    for (size_t i = 0; i < tasks.size(); ++i) {
        int greedy_node = greedy_assign[i];
        int knap_node = knap_assign[i];
        node_times_greedy[greedy_node] += tasks[i].base_time / rij[0][greedy_node];
        node_times_knap[knap_node] += tasks[i].base_time / rij[0][knap_node];
    }
    
    double new_task_size = 15.0;
    
    std::ofstream detail_file("results/detailed_scenario_analysis.csv");
    detail_file << "node_id,perf_factor,preload,greedy_current,knap_current,"
                << "greedy_proj,knap_proj,greedy_rank,knap_rank\n";
    
    std::vector<std::pair<double, int>> greedy_projections, knap_projections;
    for (size_t i = 0; i < nodes.size(); ++i) {
        double greedy_proj = node_times_greedy[i] + new_task_size / rij[0][i];
        double knap_proj = node_times_knap[i] + new_task_size / rij[0][i];
        greedy_projections.push_back({greedy_proj, i});
        knap_projections.push_back({knap_proj, i});
    }
    
    std::sort(greedy_projections.begin(), greedy_projections.end());
    std::sort(knap_projections.begin(), knap_projections.end());
    
    for (size_t i = 0; i < nodes.size(); ++i) {
        int greedy_rank = 0, knap_rank = 0;
        for (size_t j = 0; j < greedy_projections.size(); ++j) {
            if (greedy_projections[j].second == (int)i) greedy_rank = j + 1;
            if (knap_projections[j].second == (int)i) knap_rank = j + 1;
        }
        
        double greedy_proj = node_times_greedy[i] + new_task_size / rij[0][i];
        double knap_proj = node_times_knap[i] + new_task_size / rij[0][i];
        
        detail_file << i << "," << nodes[i].performance_factor << "," 
                    << nodes[i].current_load << ","
                    << std::fixed << std::setprecision(3) << node_times_greedy[i] << ","
                    << node_times_knap[i] << "," << greedy_proj << "," << knap_proj << ","
                    << greedy_rank << "," << knap_rank << "\n";
    }
    
    detail_file.close();
    std::cout << "Detailed analysis saved to 'results/detailed_scenario_analysis.csv'\n";
}

int main() {
    std::cout << "=== COMPREHENSIVE HETEROGENEOUS LOAD BALANCING ANALYSIS ===\n";
    
    // Run all analyses
    analyze_single_new_task();
    analyze_multiple_new_tasks();
    test_comprehensive_scenarios();
    analyze_detailed_scenario();
    
    std::cout << "\n=== ALL ANALYSES COMPLETED ===\n";
    std::cout << "Generated files:\n";
    std::cout << "1. results/heterogeneous_lb_rij_plot.csv - Single new task analysis\n";
    std::cout << "2. results/heterogeneous_lb_multiple_new_tasks.csv - Multiple new tasks analysis\n";
    std::cout << "3. results/comprehensive_scenario_analysis.csv - Comprehensive scenarios\n";
    std::cout << "4. results/detailed_scenario_analysis.csv - Detailed scenario analysis\n";
    std::cout << "\nRun 'python3 plot_heterogeneous_lb.py' to generate all plots.\n";
    
    return 0;
}