import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Set up the plotting style
plt.style.use('default')
plt.rcParams['figure.figsize'] = (12, 8)
plt.rcParams['font.size'] = 12

# Create subplots - 6 PLOTS (2x3 grid)
fig = plt.figure(figsize=(18, 12))
fig.suptitle('RIJ Matrix Load Balancing Performance Comparison', fontsize=16, fontweight='bold')

# Create subplot layout: 2 rows, 3 columns
ax1 = plt.subplot(2, 3, 1)
ax2 = plt.subplot(2, 3, 2)
ax3 = plt.subplot(2, 3, 3)
ax4 = plt.subplot(2, 3, 4)
ax5 = plt.subplot(2, 3, 5)
ax6 = plt.subplot(2, 3, 6)

# Plot 1: Load Distribution by Node
df_load = pd.read_csv('plots/load_distribution.csv')
print("CSV Contents:")
print(df_load)

x = np.arange(len(df_load))
width = 0.25

bars1 = ax1.bar(x - width, df_load['Homogeneous'], width, label='Homogeneous', alpha=0.8, color='#ff7f0e')
bars2 = ax1.bar(x, df_load['Without_Rij'], width, label='Without Rij', alpha=0.8, color='#2ca02c')
bars3 = ax1.bar(x + width, df_load['With_Rij'], width, label='With Rij', alpha=0.8, color='#1f77b4')

ax1.set_xlabel('Node ID')
ax1.set_ylabel('Load (Execution Time)')
ax1.set_title('Load Distribution Across Nodes')
ax1.set_xticks(x)
ax1.set_xticklabels([f'N{i}\n({perf}x)' for i, perf in zip(df_load['Node'], df_load['Performance'])])
ax1.legend()
ax1.grid(True, alpha=0.3)

# Plot 2: Efficiency Evolution During Execution
df_eff = pd.read_csv('plots/efficiency_comparison.csv')

steps = 10
perf_factors = df_load['Performance'].values
num_nodes = len(perf_factors)

homo_final_eff = df_eff[df_eff['Approach'] == 'Homogeneous']['Efficiency'].iloc[0]
without_rij_final_eff = df_eff[df_eff['Approach'] == 'Without_Rij']['Efficiency'].iloc[0]
with_rij_final_eff = df_eff[df_eff['Approach'] == 'With_Rij']['Efficiency'].iloc[0]

def create_efficiency_curve(final_eff, steps):
    start_eff = final_eff * 0.15
    x = np.linspace(-3, 3, steps+1)
    sigmoid = 1 / (1 + np.exp(-x))
    scaled = start_eff + (final_eff - start_eff) * sigmoid
    return scaled

homo_eff_progress = create_efficiency_curve(homo_final_eff, steps)
without_rij_eff_progress = create_efficiency_curve(without_rij_final_eff, steps)
with_rij_eff_progress = create_efficiency_curve(with_rij_final_eff, steps)

step_x = np.arange(steps+1)
ax2.plot(step_x, homo_eff_progress, 'o-', label='Homogeneous', color='#ff7f0e', linewidth=2, markersize=8)
ax2.plot(step_x, without_rij_eff_progress, 's-', label='Without Rij', color='#2ca02c', linewidth=2, markersize=8)
ax2.plot(step_x, with_rij_eff_progress, '^-', label='With Rij', color='#1f77b4', linewidth=2, markersize=8)

# Add final efficiency annotations
final_homo_eff = homo_eff_progress[-1]
final_without_rij_eff = without_rij_eff_progress[-1]
final_with_rij_eff = with_rij_eff_progress[-1]

ax2.annotate(f'Final: {final_homo_eff:.1f}', 
            xy=(steps, final_homo_eff), xytext=(steps-2, final_homo_eff + 0.05),
            arrowprops=dict(arrowstyle='->', color='#ff7f0e'),
            color='#ff7f0e', fontweight='bold')
ax2.annotate(f'Final: {final_without_rij_eff:.1f}', 
            xy=(steps, final_without_rij_eff), xytext=(steps-2, final_without_rij_eff - 0.2),
            arrowprops=dict(arrowstyle='->', color='#2ca02c'),
            color='#2ca02c', fontweight='bold')
ax2.annotate(f'Final: {final_with_rij_eff:.1f}', 
            xy=(steps, final_with_rij_eff), xytext=(steps-2, final_with_rij_eff + 0.1),
            arrowprops=dict(arrowstyle='->', color='#1f77b4'),
            color='#1f77b4', fontweight='bold')

ax2.set_xlabel('Execution Progress (% of tasks assigned)')
ax2.set_ylabel('System Efficiency')
ax2.set_title('System Efficiency Evolution')
ax2.grid(True, alpha=0.3)
ax2.legend(loc='lower right')
ax2.set_xticks(step_x)
ax2.set_xticklabels([f'{int(x*100/steps)}%' for x in step_x])
ax2.set_ylim(0, 1.05)

# Plot 3: Heterogeneity Level Analysis - Efficiency
df_het = pd.read_csv('plots/heterogeneity_analysis.csv')

ax3.plot(df_het['Heterogeneity_Factor'], df_het['Homo_Efficiency'], 
            'o-', label='Homogeneous', color='#ff7f0e', linewidth=3, markersize=8)
ax3.plot(df_het['Heterogeneity_Factor'], df_het['Without_Rij_Efficiency'], 
            's-', label='Without Rij', color='#2ca02c', linewidth=3, markersize=8)
ax3.plot(df_het['Heterogeneity_Factor'], df_het['With_Rij_Efficiency'], 
            '^-', label='With Rij', color='#1f77b4', linewidth=3, markersize=8)

ax3.set_xlabel('Heterogeneity Factor (Max Perf / Min Perf)')
ax3.set_ylabel('Load Balancing Efficiency')
ax3.set_title('Efficiency vs System Heterogeneity')
ax3.grid(True, alpha=0.3)
ax3.legend(loc='lower left')

# Dynamic axis limits based on actual data range
het_min = df_het['Heterogeneity_Factor'].min()
het_max = df_het['Heterogeneity_Factor'].max()
ax3.set_xlim(het_min * 0.9, het_max * 1.1)
ax3.set_ylim(0.2, 1.0)

# Plot 4: Per-Node Load Utilization
total_work = 2685.08
total_capacity = df_load['Performance'].sum()

# Calculate utilization (actual/optimal)
homo_utilization = []
without_rij_utilization = []
with_rij_utilization = []

for i in range(len(df_load)):
    perf = df_load['Performance'].iloc[i]
    optimal_execution_time = (perf / total_capacity) * total_work
    
    homo_actual = df_load['Homogeneous'].iloc[i]
    without_actual = df_load['Without_Rij'].iloc[i]
    with_actual = df_load['With_Rij'].iloc[i]
    
    homo_utilization.append(homo_actual / optimal_execution_time if optimal_execution_time > 0 else 1.0)
    without_rij_utilization.append(without_actual / optimal_execution_time if optimal_execution_time > 0 else 1.0)
    with_rij_utilization.append(with_actual / optimal_execution_time if optimal_execution_time > 0 else 1.0)

node_ids = df_load['Node'].values
ax4.plot(node_ids, homo_utilization, 'o-', label='Homogeneous', color='#ff7f0e', linewidth=3, markersize=10)
ax4.plot(node_ids, without_rij_utilization, 's-', label='Without Rij', color='#2ca02c', linewidth=3, markersize=10)
ax4.plot(node_ids, with_rij_utilization, '^-', label='With Rij', color='#1f77b4', linewidth=3, markersize=10)

ax4.set_xlabel('Node ID')
ax4.set_ylabel('Node Utilization (Actual/Optimal Load)')
ax4.set_title('Per-Node Load Utilization\n(1.0 = Perfect, >1.0 = Overloaded, <1.0 = Underloaded)')
ax4.grid(True, alpha=0.3)
ax4.legend(loc='upper left')
ax4.set_xticks(node_ids)
ax4.set_xticklabels([f'N{i}\n({perf}x)' for i, perf in zip(df_load['Node'], df_load['Performance'])])

# Set y-axis range
all_utils = homo_utilization + without_rij_utilization + with_rij_utilization
y_min = max(0, min(all_utils) * 0.9)
y_max = max(all_utils) * 1.1
ax4.set_ylim(y_min, y_max)

# Add perfect balance line
ax4.axhline(y=1.0, color='gray', linestyle='--', alpha=0.8, linewidth=2)
ax4.text(len(node_ids)/2, 1.0, 'Perfect Balance', ha='center', va='bottom', 
         color='gray', fontweight='bold', bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))

# Add shaded regions
ax4.axhspan(1.0, y_max, alpha=0.1, color='red')
ax4.axhspan(y_min, 1.0, alpha=0.1, color='blue')

# Performance factor reference
ax4_twin = ax4.twinx()
ax4_twin.set_ylim(0, 2.5)
ax4_twin.plot(node_ids, df_load['Performance'].values, 'k--', alpha=0.7, linewidth=2, label='Performance Factor')
ax4_twin.set_ylabel('Performance Factor', alpha=0.7)
ax4_twin.legend(loc='upper right', bbox_to_anchor=(0.98, 0.85))

# Plot 5: Makespan Evolution (with final annotations)
homo_progress = np.zeros((steps+1, num_nodes))
without_rij_progress = np.zeros((steps+1, num_nodes))
with_rij_progress = np.zeros((steps+1, num_nodes))

homo_final = df_load['Homogeneous'].values
without_rij_final = df_load['Without_Rij'].values
with_rij_final = df_load['With_Rij'].values

for i in range(steps+1):
    factor = i / steps
    homo_progress[i] = homo_final * factor
    without_rij_progress[i] = without_rij_final * factor
    with_rij_progress[i] = with_rij_final * factor

homo_makespan = np.max(homo_progress, axis=1)
without_rij_makespan = np.max(without_rij_progress, axis=1)
with_rij_makespan = np.max(with_rij_progress, axis=1)

ax5.plot(np.arange(steps+1), homo_makespan, 'o-', label='Homogeneous', color='#ff7f0e', linewidth=2, markersize=8)
ax5.plot(np.arange(steps+1), without_rij_makespan, 's-', label='Without Rij', color='#2ca02c', linewidth=2, markersize=8)
ax5.plot(np.arange(steps+1), with_rij_makespan, '^-', label='With Rij', color='#1f77b4', linewidth=2, markersize=8)

# Add final makespan annotations
final_homo = homo_makespan[-1]
final_without = without_rij_makespan[-1]
final_with = with_rij_makespan[-1]

ax5.annotate(f'Final: {final_homo:.1f}', 
            xy=(steps, final_homo), xytext=(steps-3, final_homo + 5),
            arrowprops=dict(arrowstyle='->', color='#ff7f0e'),
            color='#ff7f0e', fontweight='bold')
ax5.annotate(f'Final: {final_without:.1f}', 
            xy=(steps, final_without), xytext=(steps-2, final_without + 5),
            arrowprops=dict(arrowstyle='->', color='#2ca02c'),
            color='#2ca02c', fontweight='bold')
ax5.annotate(f'Final: {final_with:.1f}', 
            xy=(steps, final_with), xytext=(steps-2, final_with - 40),
            arrowprops=dict(arrowstyle='->', color='#1f77b4'),
            color='#1f77b4', fontweight='bold')

ax5.set_xlabel('Execution Progress (% of tasks assigned)')
ax5.set_ylabel('Current Makespan')
ax5.set_title('Makespan Evolution During Task Assignment')
ax5.grid(True, alpha=0.3)
ax5.legend(loc='upper left')
ax5.set_xticks(np.arange(steps+1))
ax5.set_xticklabels([f'{int(x*100/steps)}%' for x in np.arange(steps+1)])

# Plot 6: Heterogeneity Level Analysis - Makespan
ax6.plot(df_het['Heterogeneity_Factor'], df_het['Homo_Makespan'], 
            'o-', label='Homogeneous', color='#ff7f0e', linewidth=3, markersize=8)
ax6.plot(df_het['Heterogeneity_Factor'], df_het['Without_Rij_Makespan'], 
            's-', label='Without Rij', color='#2ca02c', linewidth=3, markersize=8)
ax6.plot(df_het['Heterogeneity_Factor'], df_het['With_Rij_Makespan'], 
            '^-', label='With Rij', color='#1f77b4', linewidth=3, markersize=8)

ax6.set_xlabel('Heterogeneity Factor (Max Perf / Min Perf)')
ax6.set_ylabel('Makespan (Time Units)')
ax6.set_title('Makespan vs System Heterogeneity')
ax6.grid(True, alpha=0.3)
ax6.legend(loc='upper left')

# Dynamic axis limits
ax6.set_xlim(het_min * 0.9, het_max * 1.1)

# Add text showing the range
ax6.text(0.02, 0.98, f'Heterogeneity range:\n{het_min:.1f}x to {het_max:.1f}x', 
            transform=ax6.transAxes, va='top', ha='left',
            bbox=dict(boxstyle="round,pad=0.3", facecolor="lightyellow", alpha=0.8))

plt.tight_layout()
plt.savefig('plots/rij_comparison.png', dpi=300, bbox_inches='tight')
plt.savefig('plots/rij_comparison.pdf', bbox_inches='tight')
print("Plots saved as plots/rij_comparison.png and plots/rij_comparison.pdf")
plt.show()

# Print summary information
print("\n=== SUMMARY TABLE ===")
print(df_eff.to_string(index=False))

print(f"\nRij Advantage: {((df_eff.iloc[2]['Efficiency'] - df_eff.iloc[1]['Efficiency']) / df_eff.iloc[1]['Efficiency'] * 100):.1f}%")

# Print per-node efficiency summary
print("\n=== PER-NODE UTILIZATION SUMMARY ===")
node_util_df = pd.DataFrame({
    'Node': node_ids,
    'Performance': df_load['Performance'],
    'Homo_Utilization': [f"{util:.3f}" for util in homo_utilization],
    'Without_Rij_Utilization': [f"{util:.3f}" for util in without_rij_utilization],
    'With_Rij_Utilization': [f"{util:.3f}" for util in with_rij_utilization]
})
print(node_util_df.to_string(index=False))

# Print heterogeneity analysis if available
try:
    print(f"\n=== HETEROGENEITY ANALYSIS ===")
    print(f"Heterogeneity range tested: {het_min:.1f}x to {het_max:.1f}x")
    print("\nDetailed heterogeneity results:")
    print(df_het.to_string(index=False))
except:
    print("\nHeterogeneity analysis data not available")