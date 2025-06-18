import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os

def plot_single_task_analysis():
    """Plot single new task analysis"""
    print("Creating single task analysis plots...")
    
    try:
        df = pd.read_csv('results/heterogeneous_lb_rij_plot.csv')
        
        # Create comprehensive visualization
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))
        
        bar_width = 0.35
        indices = np.arange(len(df['node_id']))
        
        # Plot 1: Stacked bar comparison with arrows and annotations
        ax1.bar(indices - bar_width/2, df['greedy_rij_current'], bar_width, 
                label='Greedy Current', color='lightblue', alpha=0.8)
        ax1.bar(indices - bar_width/2, df['greedy_rij_added'], bar_width, 
                bottom=df['greedy_rij_current'], label='Greedy Added', color='darkblue', alpha=0.8)
        
        ax1.bar(indices + bar_width/2, df['knap_rij_current'], bar_width, 
                label='Knapsack Current', color='lightgreen', alpha=0.8)
        ax1.bar(indices + bar_width/2, df['knap_rij_added'], bar_width, 
                bottom=df['knap_rij_current'], label='Knapsack Added', color='darkgreen', alpha=0.8)
        
        greedy_best = df['greedy_rij_proj'].idxmin()
        knap_best = df['knap_rij_proj'].idxmin()
        
        ax1.bar(indices[greedy_best] - bar_width/2, df['greedy_rij_proj'][greedy_best], 
                bar_width, fill=False, edgecolor='red', linewidth=3)
        ax1.bar(indices[knap_best] + bar_width/2, df['knap_rij_proj'][knap_best], 
                bar_width, fill=False, edgecolor='red', linewidth=3)
        
        # Add arrows and annotations for best choices
        ax1.annotate(f'Best Greedy\n(Node {greedy_best})\nTime: {df["greedy_rij_proj"][greedy_best]:.2f}',
                     xy=(indices[greedy_best] - bar_width/2, df['greedy_rij_proj'][greedy_best]),
                     xytext=(indices[greedy_best] - bar_width/2, df['greedy_rij_proj'][greedy_best] + 8),
                     ha='center', color='red', fontweight='bold', fontsize=9,
                     bbox=dict(boxstyle="round,pad=0.3", facecolor="yellow", alpha=0.7),
                     arrowprops=dict(facecolor='red', shrink=0.05, width=1, headwidth=6))
        
        ax1.annotate(f'Best Knapsack\n(Node {knap_best})\nTime: {df["knap_rij_proj"][knap_best]:.2f}',
                     xy=(indices[knap_best] + bar_width/2, df['knap_rij_proj'][knap_best]),
                     xytext=(indices[knap_best] + bar_width/2, df['knap_rij_proj'][knap_best] + 8),
                     ha='center', color='red', fontweight='bold', fontsize=9,
                     bbox=dict(boxstyle="round,pad=0.3", facecolor="yellow", alpha=0.7),
                     arrowprops=dict(facecolor='red', shrink=0.05, width=1, headwidth=6))
        
        ax1.set_xlabel('Node ID')
        ax1.set_ylabel('Execution Time')
        ax1.set_title("Projected Time if New Task Added \n(Greedy vs Knapsack with rij Matrix)")
        ax1.set_xticks(indices)
        ax1.set_xticklabels([f'cpu{i}' for i in df['node_id']])
        ax1.legend(loc='upper left', fontsize=8)
        ax1.grid(True, alpha=0.3)
        
        # Plot 2: Performance vs Load scatter with all labels
        colors = plt.cm.viridis(np.linspace(0, 1, len(df)))
        ax2.scatter(df['perf_factor'], df['greedy_rij_current'], 
                   c=colors, s=100, alpha=0.7, label='Greedy rij')
        ax2.scatter(df['perf_factor'], df['knap_rij_current'], 
                   c=colors, s=100, alpha=0.7, marker='s', label='Knapsack rij')
        
        # Add labels for ALL points
        for i, (x, y1, y2) in enumerate(zip(df['perf_factor'], df['greedy_rij_current'], df['knap_rij_current'])):
            # Label for greedy points
            ax2.annotate(f'cpu{i}', (x, y1), xytext=(5, 5), 
                        textcoords='offset points', fontsize=8, color='blue',
                        bbox=dict(boxstyle="round,pad=0.2", facecolor="lightblue", alpha=0.7))
            # Label for knapsack points
            ax2.annotate(f'cpu{i}', (x, y2), xytext=(5, -15), 
                        textcoords='offset points', fontsize=8, color='green',
                        bbox=dict(boxstyle="round,pad=0.2", facecolor="lightgreen", alpha=0.7))
        
        ax2.set_xlabel('Performance Factor (Higher = Faster)')
        ax2.set_ylabel('Current Load Time (rij-adjusted)')
        ax2.set_title('Performance Factor vs Current Load\n(Shows Heterogeneity Impact)')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        
        # Plot 3: Load distribution
        x_pos = np.arange(len(df))
        width = 0.35
        
        ax3.bar(x_pos - width/2, df['greedy_rij_current'], width, 
               label='Greedy rij', color='lightblue', alpha=0.8)
        ax3.bar(x_pos + width/2, df['knap_rij_current'], width, 
               label='Knapsack rij', color='lightgreen', alpha=0.8)
        
        ax3.set_xlabel('Node ID')
        ax3.set_ylabel('Current Load Time')
        ax3.set_title('Current Load Distribution Comparison\n(Before Adding New Task)')
        ax3.set_xticks(x_pos)
        ax3.set_xticklabels([f'cpu{i}' for i in df['node_id']])
        ax3.legend()
        ax3.grid(True, alpha=0.3)
        
        # Plot 4: Decision analysis with enhanced annotations
        projected_times_greedy = df['greedy_rij_proj'].values
        projected_times_knap = df['knap_rij_proj'].values
        node_labels = [f'cpu{i}\n(perf={p:.1f})' for i, p in zip(df['node_id'], df['perf_factor'])]
        
        x_pos = np.arange(len(df))
        ax4.plot(x_pos, projected_times_greedy, 'o-', color='blue', linewidth=2, 
                markersize=8, label='Greedy rij Projected', alpha=0.8)
        ax4.plot(x_pos, projected_times_knap, 's-', color='green', linewidth=2, 
                markersize=8, label='Knapsack rij Projected', alpha=0.8)
        
        # Highlight minimum points with enhanced styling
        ax4.plot(greedy_best, projected_times_greedy[greedy_best], 'o', 
                color='red', markersize=12, markerfacecolor='none', markeredgewidth=3)
        ax4.plot(knap_best, projected_times_knap[knap_best], 's', 
                color='red', markersize=12, markerfacecolor='none', markeredgewidth=3)
        
        ax4.set_xlabel('Node ID')
        ax4.set_ylabel('Projected Completion Time')
        ax4.set_title('Load Balancing Decision Analysis\n(Which Node to Choose for New Task?)')
        ax4.set_xticks(x_pos)
        ax4.set_xticklabels(node_labels, fontsize=8)
        ax4.legend()
        ax4.grid(True, alpha=0.3)
        
        # Add performance factor as secondary information
        ax4_twin = ax4.twinx()
        ax4_twin.bar(x_pos, df['perf_factor'], alpha=0.2, color='gray', width=0.5)
        ax4_twin.set_ylabel('Performance Factor', color='gray')
        ax4_twin.tick_params(axis='y', labelcolor='gray')
        
        plt.tight_layout()
        plt.savefig('outputs/single_task_analysis.png', dpi=300, bbox_inches='tight')
        plt.show()
        
        # Create additional focused side-by-side comparison
        fig, (ax_left, ax_right) = plt.subplots(1, 2, figsize=(16, 6))
        
        # Left plot: Greedy strategy with enhanced annotations
        bar_width = 0.6
        x_pos = np.arange(len(df))
        
        ax_left.bar(x_pos, df['greedy_rij_current'], bar_width, 
                   label='Current Load', color='lightblue', alpha=0.8)
        ax_left.bar(x_pos, df['greedy_rij_added'], bar_width, 
                   bottom=df['greedy_rij_current'], label='New Task Time', color='darkblue', alpha=0.8)
        
        # Highlight the best choice
        best_idx = df['greedy_rij_proj'].idxmin()
        ax_left.bar(best_idx, df['greedy_rij_proj'][best_idx], bar_width, 
                   fill=False, edgecolor='red', linewidth=4)
        
        ax_left.set_title('Greedy Strategy with rij Matrix\n(Chooses Node with Minimum Projected Time)', 
                         fontsize=12, fontweight='bold')
        ax_left.set_xlabel('Node ID')
        ax_left.set_ylabel('Execution Time (rij-adjusted)')
        ax_left.set_xticks(x_pos)
        ax_left.set_xticklabels([f'cpu{i}\n(perf={p:.1f})' for i, p in zip(df['node_id'], df['perf_factor'])], 
                               fontsize=10)
        ax_left.legend()
        ax_left.grid(True, alpha=0.3)
        
        # Add annotation for best choice
        ax_left.annotate(f'BEST CHOICE\nNode {best_idx}\nTime: {df["greedy_rij_proj"][best_idx]:.2f}',
                        xy=(best_idx, df['greedy_rij_proj'][best_idx]),
                        xytext=(best_idx + 1, df['greedy_rij_proj'][best_idx] + 10),
                        ha='center', color='red', fontweight='bold', fontsize=10,
                        bbox=dict(boxstyle="round,pad=0.3", facecolor="yellow", alpha=0.7),
                        arrowprops=dict(facecolor='red', shrink=0.05, width=2, headwidth=8))
        
        # Right plot: Knapsack strategy with enhanced annotations
        ax_right.bar(x_pos, df['knap_rij_current'], bar_width, 
                    label='Current Load', color='lightgreen', alpha=0.8)
        ax_right.bar(x_pos, df['knap_rij_added'], bar_width, 
                    bottom=df['knap_rij_current'], label='New Task Time', color='darkgreen', alpha=0.8)
        
        # Highlight the best choice
        best_idx_knap = df['knap_rij_proj'].idxmin()
        ax_right.bar(best_idx_knap, df['knap_rij_proj'][best_idx_knap], bar_width, 
                    fill=False, edgecolor='red', linewidth=4)
        
        ax_right.set_title('Knapsack Strategy with rij Matrix\n(Global Load Balancing Optimization)', 
                          fontsize=12, fontweight='bold')
        ax_right.set_xlabel('Node ID')
        ax_right.set_ylabel('Execution Time (rij-adjusted)')
        ax_right.set_xticks(x_pos)
        ax_right.set_xticklabels([f'cpu{i}\n(perf={p:.1f})' for i, p in zip(df['node_id'], df['perf_factor'])], 
                                fontsize=10)
        ax_right.legend()
        ax_right.grid(True, alpha=0.3)
        
        # Add annotation for best choice
        ax_right.annotate(f'BEST CHOICE\nNode {best_idx_knap}\nTime: {df["knap_rij_proj"][best_idx_knap]:.2f}',
                         xy=(best_idx_knap, df['knap_rij_proj'][best_idx_knap]),
                         xytext=(best_idx_knap + 1, df['knap_rij_proj'][best_idx_knap] + 10),
                         ha='center', color='red', fontweight='bold', fontsize=10,
                         bbox=dict(boxstyle="round,pad=0.3", facecolor="yellow", alpha=0.7),
                         arrowprops=dict(facecolor='red', shrink=0.05, width=2, headwidth=8))
        
        plt.tight_layout()
        plt.savefig('outputs/single_task_strategy_comparison.png', dpi=300, bbox_inches='tight')
        plt.show()
        
    except FileNotFoundError:
        print("Single task analysis file not found. Skipping...")

def plot_multiple_tasks_analysis():
    """Plot multiple new tasks analysis"""
    print("Creating multiple tasks analysis plots...")
    
    try:
        df = pd.read_csv('results/heterogeneous_lb_multiple_new_tasks.csv')
        
        unique_tasks = df['new_task_id'].unique()
        n_tasks = len(unique_tasks)
        
        # Create subplots for each new task
        fig, axes = plt.subplots(2, 3, figsize=(20, 12))
        axes = axes.flatten()
        
        for idx, task_id in enumerate(unique_tasks):
            if idx >= 6:
                break
                
            task_data = df[df['new_task_id'] == task_id]
            ax = axes[idx]
            
            bar_width = 0.35
            indices = np.arange(len(task_data))
            
            # Greedy strategy bars
            ax.bar(indices - bar_width/2, task_data['greedy_rij_current'], bar_width, 
                   label='Greedy Current' if idx == 0 else "", color='lightblue', alpha=0.8)
            ax.bar(indices - bar_width/2, task_data['greedy_rij_added'], bar_width, 
                   bottom=task_data['greedy_rij_current'], 
                   label='Greedy Added' if idx == 0 else "", color='darkblue', alpha=0.8)
            
            # Knapsack strategy bars
            ax.bar(indices + bar_width/2, task_data['knap_rij_current'], bar_width, 
                   label='Knapsack Current' if idx == 0 else "", color='lightgreen', alpha=0.8)
            ax.bar(indices + bar_width/2, task_data['knap_rij_added'], bar_width, 
                   bottom=task_data['knap_rij_current'], 
                   label='Knapsack Added' if idx == 0 else "", color='darkgreen', alpha=0.8)
            
            # Highlight best choices
            greedy_best_idx = task_data['greedy_rij_proj'].idxmin()
            knap_best_idx = task_data['knap_rij_proj'].idxmin()
            
            greedy_local_idx = list(task_data.index).index(greedy_best_idx)
            knap_local_idx = list(task_data.index).index(knap_best_idx)
            
            ax.bar(greedy_local_idx - bar_width/2, task_data.loc[greedy_best_idx, 'greedy_rij_proj'], 
                   bar_width, fill=False, edgecolor='red', linewidth=3)
            ax.bar(knap_local_idx + bar_width/2, task_data.loc[knap_best_idx, 'knap_rij_proj'], 
                   bar_width, fill=False, edgecolor='red', linewidth=3)
            
            # Get task info
            task_time = task_data['task_base_time'].iloc[0]
            greedy_best_node = task_data.loc[greedy_best_idx, 'node_id']
            knap_best_node = task_data.loc[knap_best_idx, 'node_id']
            greedy_best_time = task_data.loc[greedy_best_idx, 'greedy_rij_proj']
            knap_best_time = task_data.loc[knap_best_idx, 'knap_rij_proj']
            
            # Determine winner
            if greedy_best_time < knap_best_time:
                winner = "Greedy Wins"
                winner_color = 'blue'
            elif knap_best_time < greedy_best_time:
                winner = "Knapsack Wins"
                winner_color = 'green'
            else:
                winner = "Tie"
                winner_color = 'orange'
            
            ax.set_title(f'Task {int(task_id)} (time={task_time:.2f})\n{winner}\nGreedy→cpu{int(greedy_best_node)} vs Knap→cpu{int(knap_best_node)}', 
                        fontsize=10, color=winner_color, fontweight='bold')
            ax.set_xlabel('Node ID')
            ax.set_ylabel('Execution Time')
            ax.set_xticks(indices)
            ax.set_xticklabels([f'cpu{int(i)}' for i in task_data['node_id']])
            ax.grid(True, alpha=0.3)
            
            if idx == 0:
                ax.legend(loc='upper left', fontsize=8)
        
        # Hide unused subplots
        for idx in range(n_tasks, len(axes)):
            axes[idx].set_visible(False)
        
        plt.tight_layout()
        plt.savefig('outputs/multiple_tasks_analysis.png', dpi=300, bbox_inches='tight')
        plt.show()
        
        # Create summary plot
        summary_data = []
        for task_id in unique_tasks:
            task_data = df[df['new_task_id'] == task_id]
            task_time = task_data['task_base_time'].iloc[0]
            
            greedy_best_time = task_data['greedy_rij_proj'].min()
            knap_best_time = task_data['knap_rij_proj'].min()
            greedy_best_node = task_data.loc[task_data['greedy_rij_proj'].idxmin(), 'node_id']
            knap_best_node = task_data.loc[task_data['knap_rij_proj'].idxmin(), 'node_id']
            
            summary_data.append({
                'task_id': int(task_id),
                'task_time': task_time,
                'greedy_best_time': greedy_best_time,
                'knap_best_time': knap_best_time,
                'greedy_best_node': int(greedy_best_node),
                'knap_best_node': int(knap_best_node),
                'greedy_wins': greedy_best_time < knap_best_time,
                'time_difference': abs(greedy_best_time - knap_best_time)
            })
        
        summary_df = pd.DataFrame(summary_data)
        
        # Summary comparison plot
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
        
        x_pos = np.arange(len(summary_df))
        width = 0.35
        
        bars1 = ax1.bar(x_pos - width/2, summary_df['greedy_best_time'], width, 
                        label='Greedy Best Time', color='lightblue', alpha=0.8)
        bars2 = ax1.bar(x_pos + width/2, summary_df['knap_best_time'], width, 
                        label='Knapsack Best Time', color='lightgreen', alpha=0.8)
        
        # Highlight winners
        for i, (greedy_time, knap_time, greedy_wins) in enumerate(zip(summary_df['greedy_best_time'], 
                                                                     summary_df['knap_best_time'],
                                                                     summary_df['greedy_wins'])):
            if greedy_wins:
                bars1[i].set_edgecolor('gold')
                bars1[i].set_linewidth(3)
            else:
                bars2[i].set_edgecolor('gold')
                bars2[i].set_linewidth(3)
        
        ax1.set_xlabel('Task ID')
        ax1.set_ylabel('Best Completion Time')
        ax1.set_title('Best Completion Time Comparison\n(Gold border = Winner)')
        ax1.set_xticks(x_pos)
        ax1.set_xticklabels([f'T{int(tid)}' for tid in summary_df['task_id']])
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        
        # Node selection comparison
        ax2.scatter(summary_df['greedy_best_node'], summary_df['task_time'], 
                   s=100, alpha=0.7, label='Greedy Choice', color='blue')
        ax2.scatter(summary_df['knap_best_node'], summary_df['task_time'], 
                   s=100, alpha=0.7, label='Knapsack Choice', color='green', marker='s')
        
        ax2.set_xlabel('Selected Node ID')
        ax2.set_ylabel('Task Base Time')
        ax2.set_title('Node Selection Strategy Comparison')
        ax2.set_xticks(range(6))
        ax2.set_xticklabels([f'cpu{i}' for i in range(6)])
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig('outputs/multiple_tasks_summary.png', dpi=300, bbox_inches='tight')
        plt.show()
        
    except FileNotFoundError:
        print("Multiple tasks analysis file not found. Skipping...")

def plot_comprehensive_analysis():
    """Plot comprehensive scenario analysis"""
    print("Creating comprehensive analysis plots...")
    
    try:
        df = pd.read_csv('results/comprehensive_scenario_analysis.csv')
        
        # Create performance profile names
        perf_names = ['Moderate', 'High_Hetero', 'Homogeneous', 'Gradual', 'Mixed']
        preload_names = ['None', 'Node5_Heavy', 'Node0_Heavy', 'Node2_Heavy', 'First3_Heavy', 'Increasing']
        
        df['Perf_Name'] = df['Perf_Profile'].map(lambda x: perf_names[x])
        df['Preload_Name'] = df['Preload_Profile'].map(lambda x: preload_names[x])
        
        # Strategy comparison heatmap
        fig, axes = plt.subplots(2, 3, figsize=(18, 12))
        axes = axes.flatten()
        
        for i, perf_name in enumerate(['Moderate', 'High_Hetero', 'Homogeneous', 'Gradual', 'Mixed']):
            if i >= 6:
                break
                
            perf_data = df[df['Perf_Name'] == perf_name]
            
            pivot_data = perf_data.pivot_table(
                values='Winner', 
                index='Preload_Name', 
                columns='Task_Size',
                aggfunc=lambda x: (x == 'Knapsack').mean()
            )
            
            ax = axes[i]
            sns.heatmap(pivot_data, annot=True, fmt='.2f', cmap='RdYlBu_r', 
                       center=0.5, ax=ax, cbar_kws={'label': 'Knapsack Win Rate'})
            ax.set_title(f'{perf_name} Performance Profile\n(1.0 = Knapsack always wins)')
            ax.set_xlabel('Task Size')
            ax.set_ylabel('Preload Scenario')
        
        if len(axes) > 5:
            axes[5].set_visible(False)
        
        plt.tight_layout()
        plt.savefig('outputs/comprehensive_heatmap.png', dpi=300, bbox_inches='tight')
        plt.show()
        
        # Performance analysis plots
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))
        
        # Time difference distribution
        ax1.hist(df['Time_Diff'], bins=50, alpha=0.7, edgecolor='black')
        ax1.set_xlabel('Time Difference (|Greedy - Knapsack|)')
        ax1.set_ylabel('Frequency')
        ax1.set_title('Distribution of Performance Differences')
        ax1.axvline(df['Time_Diff'].mean(), color='red', linestyle='--', 
                    label=f'Mean: {df["Time_Diff"].mean():.2f}')
        ax1.legend()
        
        # Performance by task size
        task_size_analysis = df.groupby('Task_Size').agg({
            'Time_Diff': 'mean',
            'Winner': lambda x: (x == 'Knapsack').mean()
        }).reset_index()
        
        ax2_twin = ax2.twinx()
        ax2.bar(task_size_analysis['Task_Size'], task_size_analysis['Time_Diff'], 
                alpha=0.7, color='skyblue', label='Avg Time Diff')
        ax2_twin.plot(task_size_analysis['Task_Size'], task_size_analysis['Winner'], 
                      'ro-', linewidth=2, label='Knapsack Win Rate')
        
        ax2.set_xlabel('Task Size')
        ax2.set_ylabel('Average Time Difference', color='blue')
        ax2_twin.set_ylabel('Knapsack Win Rate', color='red')
        ax2.set_title('Performance vs Task Size')
        ax2.legend(loc='upper left')
        ax2_twin.legend(loc='upper right')
        
        # Performance by heterogeneity
        hetero_analysis = df.groupby('Perf_Name').agg({
            'Time_Diff': 'mean',
            'Winner': lambda x: (x == 'Knapsack').mean(),
            'Greedy_Best_Time': 'mean',
            'Knap_Best_Time': 'mean'
        }).reset_index()
        
        x_pos = np.arange(len(hetero_analysis))
        width = 0.35
        
        ax3.bar(x_pos - width/2, hetero_analysis['Greedy_Best_Time'], width, 
                label='Greedy Avg Time', alpha=0.8, color='lightcoral')
        ax3.bar(x_pos + width/2, hetero_analysis['Knap_Best_Time'], width, 
                label='Knapsack Avg Time', alpha=0.8, color='lightgreen')
        
        ax3.set_xlabel('Performance Profile')
        ax3.set_ylabel('Average Best Time')
        ax3.set_title('Strategy Performance by Heterogeneity Level')
        ax3.set_xticks(x_pos)
        ax3.set_xticklabels(hetero_analysis['Perf_Name'], rotation=45)
        ax3.legend()
        
        # Node selection patterns
        node_selection = df.groupby(['Perf_Name', 'Greedy_Best_Node']).size().unstack(fill_value=0)
        node_selection_pct = node_selection.div(node_selection.sum(axis=1), axis=0) * 100
        
        node_selection_pct.plot(kind='bar', stacked=True, ax=ax4, alpha=0.7)
        ax4.set_xlabel('Performance Profile')
        ax4.set_ylabel('Selection Percentage')
        ax4.set_title('Node Selection Patterns (Greedy Strategy)')
        ax4.legend(title='Node ID', bbox_to_anchor=(1.05, 1), loc='upper left')
        ax4.tick_params(axis='x', rotation=45)
        
        plt.tight_layout()
        plt.savefig('outputs/comprehensive_performance.png', dpi=300, bbox_inches='tight')
        plt.show()
        
    except FileNotFoundError:
        print("Comprehensive analysis file not found. Skipping...")

def plot_detailed_scenario():
    """Plot detailed scenario analysis"""
    print("Creating detailed scenario analysis plots...")
    
    try:
        df = pd.read_csv('results/detailed_scenario_analysis.csv')
        
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 10))
        
        x_pos = np.arange(len(df))
        width = 0.35
        
        # Current loads comparison
        ax1.bar(x_pos - width/2, df['greedy_current'], width, 
                label='Greedy Current', alpha=0.8, color='lightblue')
        ax1.bar(x_pos + width/2, df['knap_current'], width, 
                label='Knapsack Current', alpha=0.8, color='lightgreen')
        
        ax1.set_xlabel('Node ID')
        ax1.set_ylabel('Current Load')
        ax1.set_title('Current Load Distribution')
        ax1.set_xticks(x_pos)
        ax1.set_xticklabels([f'cpu{i}\n(perf={p:.1f})' for i, p in 
                            zip(df['node_id'], df['perf_factor'])])
        ax1.legend()
        
        # Projected loads comparison
        ax2.bar(x_pos - width/2, df['greedy_proj'], width, 
                label='Greedy Projected', alpha=0.8, color='darkblue')
        ax2.bar(x_pos + width/2, df['knap_proj'], width, 
                label='Knapsack Projected', alpha=0.8, color='darkgreen')
        
        greedy_best = df['greedy_rank'].idxmin()
        knap_best = df['knap_rank'].idxmin()
        
        ax2.bar(greedy_best - width/2, df.loc[greedy_best, 'greedy_proj'], width,
                fill=False, edgecolor='red', linewidth=3)
        ax2.bar(knap_best + width/2, df.loc[knap_best, 'knap_proj'], width,
                fill=False, edgecolor='red', linewidth=3)
        
        ax2.set_xlabel('Node ID')
        ax2.set_ylabel('Projected Load')
        ax2.set_title('Projected Load After Adding New Task')
        ax2.set_xticks(x_pos)
        ax2.set_xticklabels([f'cpu{i}' for i in df['node_id']])
        ax2.legend()
        
        # Performance vs load scatter
        ax3.scatter(df['perf_factor'], df['greedy_current'], 
                   s=100, alpha=0.7, label='Greedy', color='blue')
        ax3.scatter(df['perf_factor'], df['knap_current'], 
                   s=100, alpha=0.7, label='Knapsack', color='green', marker='s')
        
        for i, (x, y1, y2) in enumerate(zip(df['perf_factor'], 
                                           df['greedy_current'], 
                                           df['knap_current'])):
            ax3.annotate(f'cpu{i}', (x, y1), xytext=(5, 5), 
                        textcoords='offset points', fontsize=8, color='blue')
            ax3.annotate(f'cpu{i}', (x, y2), xytext=(5, -15), 
                        textcoords='offset points', fontsize=8, color='green')
        
        ax3.set_xlabel('Performance Factor')
        ax3.set_ylabel('Current Load')
        ax3.set_title('Performance vs Load Trade-off')
        ax3.legend()
        ax3.grid(True, alpha=0.3)
        
        # Ranking comparison
        strategies = ['Greedy', 'Knapsack']
        ranks = [df['greedy_rank'].values, df['knap_rank'].values]
        colors = ['lightblue', 'lightgreen']
        
        for i, (strategy, rank_data, color) in enumerate(zip(strategies, ranks, colors)):
            ax4.bar(x_pos + i*width, rank_data, width, label=strategy, 
                   color=color, alpha=0.8)
        
        ax4.set_xlabel('Node ID')
        ax4.set_ylabel('Rank (1=Best Choice)')
        ax4.set_title('Node Ranking by Strategy')
        ax4.set_xticks(x_pos + width/2)
        ax4.set_xticklabels([f'cpu{i}' for i in df['node_id']])
        ax4.legend()
        ax4.set_ylim(0, 7)
        
        plt.tight_layout()
        plt.savefig('outputs/detailed_scenario.png', dpi=300, bbox_inches='tight')
        plt.show()
        
    except FileNotFoundError:
        print("Detailed scenario file not found. Skipping...")

def print_comprehensive_summary():
    """Print comprehensive analysis summary"""
    try:
        df = pd.read_csv('results/comprehensive_scenario_analysis.csv')
        
        perf_names = ['Moderate', 'High_Hetero', 'Homogeneous', 'Gradual', 'Mixed']
        df['Perf_Name'] = df['Perf_Profile'].map(lambda x: perf_names[x])
        
        print("\n=== COMPREHENSIVE ANALYSIS SUMMARY ===")
        print(f"Total scenarios tested: {len(df)}")
        
        overall_winner = df['Winner'].value_counts()
        print(f"\nOverall strategy performance:")
        for strategy, count in overall_winner.items():
            percentage = (count / len(df)) * 100
            print(f"  {strategy}: {count} wins ({percentage:.1f}%)")
        
        print(f"\nAverage performance difference: {df['Time_Diff'].mean():.2f}")
        print(f"Maximum performance difference: {df['Time_Diff'].max():.2f}")
        
        print(f"\nPerformance by heterogeneity level:")
        hetero_summary = df.groupby('Perf_Name')['Winner'].apply(
            lambda x: f"Knapsack: {(x=='Knapsack').mean():.2f}, Greedy: {(x=='Greedy').mean():.2f}"
        )
        for perf_type, summary in hetero_summary.items():
            print(f"  {perf_type}: {summary}")
            
    except FileNotFoundError:
        print("Comprehensive analysis file not found. Cannot print summary.")

def main():
    """Main function to run all visualizations"""
    print("=== COMPREHENSIVE HETEROGENEOUS LOAD BALANCING VISUALIZATION ===\n")
    
    # Create all plots
    plot_single_task_analysis()
    plot_multiple_tasks_analysis()
    plot_comprehensive_analysis()
    plot_detailed_scenario()
    
    # Print summary
    print_comprehensive_summary()
    
    print("\n=== ALL VISUALIZATIONS COMPLETED ===")
    print("Generated plots:")
    print("1. outputs/single_task_analysis.png - Single new task analysis")
    print("2. outputs/multiple_tasks_analysis.png - Multiple new tasks comparison")
    print("3. outputs/multiple_tasks_summary.png - Multiple tasks summary")
    print("4. outputs/comprehensive_heatmap.png - Strategy comparison heatmap")
    print("5. outputs/comprehensive_performance.png - Performance analysis")
    print("6. outputs/detailed_scenario.png - Detailed scenario analysis")

if __name__ == "__main__":
    main()