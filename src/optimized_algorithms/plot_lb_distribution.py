import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib.patches import Rectangle

# Read and clean the data
def clean_data(filename):
    """Clean the CSV data by handling extra commas"""
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    cleaned_lines = []
    for line in lines:
        # Remove extra commas by splitting and rejoining properly
        parts = [part.strip() for part in line.strip().split(',') if part.strip()]
        
        # Should have exactly 5 parts: Config, Algorithm, Efficiency, CutFraction, HaloExchange
        if len(parts) == 5:
            cleaned_lines.append(','.join(parts))
        elif len(parts) == 4:  # Missing one field, likely CutFraction
            # Insert empty CutFraction
            cleaned_lines.append(','.join(parts[:2] + [''] + parts[2:]))
    
    # Write cleaned data to temporary file
    with open('temp_cleaned.csv', 'w') as f:
        f.write('Config,Algorithm,Efficiency,CutFraction,HaloExchange\n')
        for line in cleaned_lines[1:]:  # Skip original header
            f.write(line + '\n')
    
    return pd.read_csv('temp_cleaned.csv')

# Use the cleaned data from the actual file
df = clean_data('summary_metrics.csv')

# Convert to numeric (should work properly now)
df['CutFraction'] = pd.to_numeric(df['CutFraction'], errors='coerce')
df['Efficiency'] = pd.to_numeric(df['Efficiency'], errors='coerce')
df['HaloExchange'] = pd.to_numeric(df['HaloExchange'], errors='coerce')

# Set up the plotting style
plt.style.use('default')
sns.set_palette("husl")
colors = {'knapsack': '#e74c3c', 'sfc': '#3498db', 'hilbert': '#2ecc71'}

def create_comprehensive_analysis():
    """Create a comprehensive analysis of all configurations"""
    
    # Create figure with subplots
    fig = plt.figure(figsize=(20, 16))
    
    # 1. Load Balance Efficiency Comparison (Top Left)
    ax1 = plt.subplot(3, 3, 1)
    configs = df['Config'].unique()
    algorithms = df['Algorithm'].unique()
    
    x = np.arange(len(configs))
    width = 0.25
    
    for i, alg in enumerate(algorithms):
        data = df[df['Algorithm'] == alg]['Efficiency']
        bars = ax1.bar(x + i*width, data, width, label=alg.title(), 
                      color=colors[alg], alpha=0.8)
        
        # Add value labels on bars
        for bar, val in zip(bars, data):
            if not pd.isna(val):
                ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01,
                        f'{val:.3f}', ha='center', va='bottom', fontsize=8, rotation=45)
    
    ax1.set_xlabel('Configuration')
    ax1.set_ylabel('Load Balance Efficiency')
    ax1.set_title('Load Balance Efficiency by Configuration', fontweight='bold')
    ax1.set_xticks(x + width)
    ax1.set_xticklabels(configs, rotation=45, ha='right')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    ax1.set_ylim(0, 1.1)
    
    # 2. Cut Fraction Comparison (Top Center)
    ax2 = plt.subplot(3, 3, 2)
    
    for i, alg in enumerate(algorithms):
        data = df[df['Algorithm'] == alg]['CutFraction']
        bars = ax2.bar(x + i*width, data, width, label=alg.title(), 
                      color=colors[alg], alpha=0.8)
        
        # Add value labels on bars
        for bar, val in zip(bars, data):
            if not pd.isna(val):
                ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02,
                        f'{val:.3f}', ha='center', va='bottom', fontsize=8, rotation=45)
    
    ax2.set_xlabel('Configuration')
    ax2.set_ylabel('Cut Fraction (Lower is Better)')
    ax2.set_title('Graph Cut Fraction by Configuration', fontweight='bold')
    ax2.set_xticks(x + width)
    ax2.set_xticklabels(configs, rotation=45, ha='right')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # 3. Halo Exchange Comparison (Top Right)
    ax3 = plt.subplot(3, 3, 3)
    
    for i, alg in enumerate(algorithms):
        data = df[df['Algorithm'] == alg]['HaloExchange']
        bars = ax3.bar(x + i*width, data, width, label=alg.title(), 
                      color=colors[alg], alpha=0.8)
        
        # Add value labels on bars
        for bar, val in zip(bars, data):
            if not pd.isna(val):
                ax3.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(data)*0.02,
                        f'{int(val):,}', ha='center', va='bottom', fontsize=8, rotation=45)
    
    ax3.set_xlabel('Configuration')
    ax3.set_ylabel('Total Halo Exchange Volume')
    ax3.set_title('Communication Overhead by Configuration', fontweight='bold')
    ax3.set_xticks(x + width)
    ax3.set_xticklabels(configs, rotation=45, ha='right')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    # 4. Algorithm Performance Radar Chart (Middle Left)
    ax4 = plt.subplot(3, 3, 4, projection='polar')
    
    # Normalize metrics for radar chart (higher is better for all)
    metrics = ['Load Balance', 'Spatial Locality', 'Communication']
    
    for alg in algorithms:
        alg_data = df[df['Algorithm'] == alg]
        
        # Calculate average performance (normalize so higher is better)
        avg_efficiency = alg_data['Efficiency'].mean()
        avg_locality = 1 - alg_data['CutFraction'].mean()  # Invert cut fraction
        avg_comm = 1 - (alg_data['HaloExchange'].mean() / alg_data['HaloExchange'].max())  # Invert halo exchange
        
        values = [avg_efficiency, avg_locality, avg_comm]
        values += values[:1]  # Complete the circle
        
        angles = np.linspace(0, 2*np.pi, len(metrics), endpoint=False).tolist()
        angles += angles[:1]
        
        ax4.plot(angles, values, 'o-', linewidth=2, label=alg.title(), color=colors[alg])
        ax4.fill(angles, values, alpha=0.25, color=colors[alg])
    
    ax4.set_xticks(angles[:-1])
    ax4.set_xticklabels(metrics)
    ax4.set_ylim(0, 1)
    ax4.set_title('Algorithm Performance Profile', fontweight='bold', pad=20)
    ax4.legend(loc='upper right', bbox_to_anchor=(1.3, 1.0))
    
    # 5. Configuration Impact Analysis (Middle Center)
    ax5 = plt.subplot(3, 3, 5)
    
    # Create heatmap of algorithm performance across configurations
    pivot_efficiency = df.pivot(index='Config', columns='Algorithm', values='Efficiency')
    sns.heatmap(pivot_efficiency, annot=True, fmt='.3f', cmap='RdYlGn', 
                ax=ax5, cbar_kws={'label': 'Efficiency'})
    ax5.set_title('Load Balance Efficiency Heatmap', fontweight='bold')
    ax5.set_xlabel('Algorithm')
    ax5.set_ylabel('Configuration')
    
    # 6. Trade-off Analysis (Middle Right)
    ax6 = plt.subplot(3, 3, 6)
    
    for alg in algorithms:
        alg_data = df[df['Algorithm'] == alg]
        ax6.scatter(alg_data['CutFraction'], alg_data['Efficiency'], 
                   label=alg.title(), color=colors[alg], s=100, alpha=0.7)
        
        # Add configuration labels
        for i, row in alg_data.iterrows():
            ax6.annotate(row['Config'], (row['CutFraction'], row['Efficiency']),
                        xytext=(5, 5), textcoords='offset points', fontsize=8)
    
    ax6.set_xlabel('Cut Fraction (Lower is Better)')
    ax6.set_ylabel('Load Balance Efficiency (Higher is Better)')
    ax6.set_title('Load Balance vs Spatial Locality Trade-off', fontweight='bold')
    ax6.legend()
    ax6.grid(True, alpha=0.3)
    
    # 7. Configuration Ranking (Bottom Left)
    ax7 = plt.subplot(3, 3, 7)
    
    # Calculate composite score for each algorithm-config combination
    df_norm = df.copy()
    df_norm['Efficiency_norm'] = df_norm['Efficiency']
    df_norm['Locality_norm'] = 1 - df_norm['CutFraction']  # Invert so higher is better
    df_norm['Comm_norm'] = 1 - (df_norm['HaloExchange'] / df_norm['HaloExchange'].max())
    df_norm['Composite_Score'] = (df_norm['Efficiency_norm'] + 
                                 df_norm['Locality_norm'] + 
                                 df_norm['Comm_norm']) / 3
    
    # Plot composite scores
    for i, alg in enumerate(algorithms):
        data = df_norm[df_norm['Algorithm'] == alg]['Composite_Score']
        bars = ax7.bar(x + i*width, data, width, label=alg.title(), 
                      color=colors[alg], alpha=0.8)
        
        # Add value labels
        for bar, val in zip(bars, data):
            if not pd.isna(val):
                ax7.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01,
                        f'{val:.3f}', ha='center', va='bottom', fontsize=8, rotation=45)
    
    ax7.set_xlabel('Configuration')
    ax7.set_ylabel('Composite Performance Score')
    ax7.set_title('Overall Algorithm Ranking', fontweight='bold')
    ax7.set_xticks(x + width)
    ax7.set_xticklabels(configs, rotation=45, ha='right')
    ax7.legend()
    ax7.grid(True, alpha=0.3)
    
    # 8. Detailed Metrics Table (Bottom Center)
    ax8 = plt.subplot(3, 3, 8)
    ax8.axis('off')
    
    # Create summary statistics table
    summary_stats = []
    for alg in algorithms:
        alg_data = df[df['Algorithm'] == alg]
        stats = {
            'Algorithm': alg.title(),
            'Avg Efficiency': f"{alg_data['Efficiency'].mean():.3f}",
            'Avg Cut Fraction': f"{alg_data['CutFraction'].mean():.3f}",
            'Avg Halo Exchange': f"{alg_data['HaloExchange'].mean():,.0f}",
        }
        summary_stats.append(stats)
    
    summary_df = pd.DataFrame(summary_stats)
    
    # Create table
    table = ax8.table(cellText=summary_df.values,
                     colLabels=summary_df.columns,
                     cellLoc='center',
                     loc='center',
                     colColours=['lightgray']*len(summary_df.columns))
    
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1, 2)
    ax8.set_title('Algorithm Summary Statistics', fontweight='bold', pad=20)
    
    # 9. Key Insights (Bottom Right)
    ax9 = plt.subplot(3, 3, 9)
    ax9.axis('off')
    
    # Find best performing algorithm for each metric
    best_efficiency = df.loc[df['Efficiency'].idxmax()]
    best_locality = df.loc[df['CutFraction'].idxmin()]
    best_communication = df.loc[df['HaloExchange'].idxmin()]
    
    insights_text = f"""KEY INSIGHTS

BEST PERFORMERS:
• Load Balance: {best_efficiency['Algorithm'].title()} 
  ({best_efficiency['Config']}: {best_efficiency['Efficiency']:.3f})

• Spatial Locality: {best_locality['Algorithm'].title()}
  ({best_locality['Config']}: {best_locality['CutFraction']:.3f})

• Communication: {best_communication['Algorithm'].title()}
  ({best_communication['Config']}: {best_communication['HaloExchange']:,.0f})

ALGORITHM STRENGTHS:
• Knapsack: Perfect load balance
• SFC: Best spatial locality
• Hilbert: Balanced performance

CONFIGURATION IMPACT:
• SmallGrid: Increases communication
• Thin: Challenges load balancing
• FewerRanks: Reduces efficiency
"""
    
    ax9.text(0.05, 0.95, insights_text, transform=ax9.transAxes, 
            fontsize=11, verticalalignment='top', fontfamily='monospace',
            bbox=dict(boxstyle="round,pad=0.5", facecolor="lightblue", alpha=0.8))
    
    plt.tight_layout()
    plt.savefig('load_balance_sweep_analysis.png', dpi=300, bbox_inches='tight')
    print("Comprehensive analysis saved as 'load_balance_sweep_analysis.png'")
    plt.show()

def print_summary_table():
    """Print a formatted summary table"""
    print("\n" + "="*80)
    print("LOAD BALANCING ALGORITHM PERFORMANCE SUMMARY")
    print("="*80)
    
    for config in df['Config'].unique():
        print(f"\n{config.upper()} CONFIGURATION:")
        print("-" * 50)
        config_data = df[df['Config'] == config]
        
        for _, row in config_data.iterrows():
            print(f"{row['Algorithm'].title():>10}: "
                  f"Efficiency={row['Efficiency']:.3f}, "
                  f"CutFraction={row['CutFraction']:.3f}, "
                  f"HaloExchange={row['HaloExchange']:>8,.0f}")

if __name__ == "__main__":
    # Create all visualizations
    create_comprehensive_analysis()
    print_summary_table()
    
    # Save corrected data
    df.to_csv('summary_metrics.csv', index=False)
    print("\nCorrected data saved as 'summary_metrics.csv'")

# import os
# import numpy as np
# import matplotlib.pyplot as plt
# from mpl_toolkits.axes_grid1 import make_axes_locatable

# class LoadBalanceAnalyzer:
#     """Comprehensive load balance analysis and visualization"""
    
#     def __init__(self):
#         self.algorithms = ['Knapsack', 'SFC', 'Hilbert']
#         self.colors = ['#e74c3c', '#3498db', '#2ecc71']  # Red, Blue, Green
#         self.data = {}
        
#     def load_distribution_data(self, filename):
#         """Load distribution mapping data"""
#         try:
#             raw_data = np.loadtxt(filename, comments='#')
#             # box_id, rank, x, y, z, cost
#             ranks = raw_data[:, 1].astype(int)
#             costs = raw_data[:, 5]
            
#             # Calculate load balance efficiency
#             rank_costs = {}
#             for rank, cost in zip(ranks, costs):
#                 rank_costs[rank] = rank_costs.get(rank, 0) + cost
            
#             cost_values = np.array(list(rank_costs.values()))
#             efficiency = cost_values.mean() / cost_values.max() if cost_values.max() > 0 else 0
            
#             # For 3D visualization
#             x_coords = raw_data[:, 2].astype(int)
#             y_coords = raw_data[:, 3].astype(int)
#             z_coords = raw_data[:, 4].astype(int)
            
#             return {
#                 'ranks': ranks,
#                 'costs': costs,
#                 'efficiency': efficiency,
#                 'coords': (x_coords, y_coords, z_coords),
#                 'rank_costs': rank_costs
#             }
#         except Exception as e:
#             print(f"Error loading {filename}: {e}")
#             return None
    
#     def load_communication_data(self):
#         """Load graph cut and halo exchange data"""
#         comm_data = {}
        
#         for i, alg in enumerate(self.algorithms):
#             # Graph cuts
#             cut_file = f'LBC_{alg.lower()}_graph_cut.txt'
#             try:
#                 with open(cut_file) as f:
#                     for line in f:
#                         if line.strip().startswith('CutFraction'):
#                             comm_data[alg] = {'cut_fraction': float(line.strip().split()[1])}
#                             break
#             except:
#                 comm_data[alg] = {'cut_fraction': 0.0}
            
#             # Halo exchange
#             halo_file = f'LBC_{alg.lower()}_halo_exchange.txt'
#             try:
#                 volumes = []
#                 with open(halo_file) as f:
#                     for line in f:
#                         if 'Halo exchange volume' in line:
#                             volumes.append(int(line.strip().split()[-1]))
#                 comm_data[alg]['halo_total'] = sum(volumes) if volumes else 0
#                 comm_data[alg]['halo_volumes'] = volumes
#             except:
#                 comm_data[alg]['halo_total'] = 0
#                 comm_data[alg]['halo_volumes'] = []
        
#         return comm_data
    
#     def plot_3d_distribution(self, data, title, ax):
#         """Plot 3D distribution with better styling"""
#         if data is None:
#             ax.text(0.5, 0.5, 0.5, 'No Data', ha='center', va='center')
#             return
            
#         x, y, z = data['coords']
#         ranks = data['ranks']
        
#         # Use consistent colormap
#         scatter = ax.scatter(x, y, z, c=ranks, cmap='tab20', s=8, alpha=0.7)
#         ax.set_xlabel('X')
#         ax.set_ylabel('Y') 
#         ax.set_zlabel('Z')
#         ax.set_title(title, fontsize=12, fontweight='bold')
        
#         # Clean up the plot
#         ax.grid(True, alpha=0.3)
        
#         return scatter
    
#     def create_comprehensive_plot(self):
#         """Create a comprehensive analysis plot"""
#         # Load all data
#         distribution_data = {}
#         for alg in self.algorithms:
#             filename = f'LBC_{alg.lower()}.txt'
#             distribution_data[alg] = self.load_distribution_data(filename)
        
#         comm_data = self.load_communication_data()
        
#         # Create figure with subplots
#         fig = plt.figure(figsize=(16, 12))
        
#         # 3D Distribution plots (top row)
#         for i, alg in enumerate(self.algorithms):
#             ax = fig.add_subplot(3, 3, i+1, projection='3d')
#             scatter = self.plot_3d_distribution(
#                 distribution_data[alg], 
#                 f'{alg} Distribution', 
#                 ax
#             )
#             if i == 2 and scatter:  # Add colorbar to last plot
#                 plt.colorbar(scatter, ax=ax, shrink=0.5, label='Rank')
        
#         # Load Balance Efficiency (middle left)
#         ax4 = fig.add_subplot(3, 3, 4)
#         efficiencies = [distribution_data[alg]['efficiency'] if distribution_data[alg] else 0 for alg in self.algorithms]
#         bars = ax4.bar(self.algorithms, efficiencies, color=self.colors, alpha=0.7)
#         ax4.set_ylabel('Load Balance Efficiency')
#         ax4.set_title('Load Balance Efficiency', fontweight='bold')
#         ax4.set_ylim(0, 1)
#         ax4.grid(True, alpha=0.3)
        
#         # Add value labels on bars
#         for bar, eff in zip(bars, efficiencies):
#             ax4.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01, 
#                     f'{eff:.3f}', ha='center', va='bottom', fontweight='bold')
        
#         # Graph Cut Fraction (middle center)
#         ax5 = fig.add_subplot(3, 3, 5)
#         cut_fractions = [comm_data[alg]['cut_fraction'] if alg in comm_data and not np.isnan(comm_data[alg]['cut_fraction']) else 0 for alg in self.algorithms]
#         bars = ax5.bar(self.algorithms, cut_fractions, color=self.colors, alpha=0.7)
#         ax5.set_ylabel('Cut Fraction')
#         ax5.set_title('Graph Cut Fraction', fontweight='bold')
#         ax5.set_ylim(0, 1)
#         ax5.grid(True, alpha=0.3)
        
#         # Add value labels
#         for bar, frac in zip(bars, cut_fractions):
#             ax5.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01, 
#                     f'{frac:.3f}', ha='center', va='bottom', fontweight='bold')
        
#         # Total Halo Exchange (middle right)
#         ax6 = fig.add_subplot(3, 3, 6)
#         halo_totals = [comm_data[alg]['halo_total'] if alg in comm_data else 0 for alg in self.algorithms]
#         bars = ax6.bar(self.algorithms, halo_totals, color=self.colors, alpha=0.7)
#         ax6.set_ylabel('Total Halo Exchange Volume')
#         ax6.set_title('Total Halo Exchange', fontweight='bold')
#         ax6.grid(True, alpha=0.3)
        
#         # Add value labels
#         for bar, total in zip(bars, halo_totals):
#             ax6.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(halo_totals)*0.01, 
#                     f'{total:,}', ha='center', va='bottom', fontweight='bold', rotation=45)
        
#         # Load Distribution Comparison (bottom left)
#         ax7 = fig.add_subplot(3, 3, 7)
#         for i, alg in enumerate(self.algorithms):
#             if distribution_data[alg]:
#                 costs = list(distribution_data[alg]['rank_costs'].values())
#                 ax7.hist(costs, bins=20, alpha=0.6, label=alg, color=self.colors[i])
#         ax7.set_xlabel('Load per Rank')
#         ax7.set_ylabel('Frequency')
#         ax7.set_title('Load Distribution', fontweight='bold')
#         ax7.legend()
#         ax7.grid(True, alpha=0.3)
        
#         # Halo Exchange Distribution (bottom center)
#         ax8 = fig.add_subplot(3, 3, 8)
#         halo_data = [comm_data[alg]['halo_volumes'] for alg in self.algorithms if comm_data[alg]['halo_volumes']]
#         if halo_data:
#             ax8.boxplot(halo_data, labels=[alg for alg in self.algorithms if comm_data[alg]['halo_volumes']])
#         ax8.set_ylabel('Halo Exchange Volume')
#         ax8.set_title('Halo Exchange Distribution', fontweight='bold')
#         ax8.grid(True, alpha=0.3)
        
#         # Summary Statistics (bottom right)
#         ax9 = fig.add_subplot(3, 3, 9)
#         ax9.axis('off')
        
#         # Create summary table
#         summary_text = "PERFORMANCE SUMMARY\n" + "="*25 + "\n\n"
        
#         for i, alg in enumerate(self.algorithms):
#             eff = efficiencies[i]
#             cut = cut_fractions[i]
#             halo = halo_totals[i]
            
#             summary_text += f"{alg}:\n"
#             summary_text += f"  Load Balance: {eff:.3f}\n"
#             summary_text += f"  Cut Fraction: {cut:.3f}\n"
#             summary_text += f"  Halo Volume: {halo:,}\n\n"
        
#         # Determine best algorithm
#         best_balance = self.algorithms[np.argmax(efficiencies)]
#         best_locality = self.algorithms[np.argmin(cut_fractions)]
#         best_comm = self.algorithms[np.argmin(halo_totals)]
        
#         summary_text += f"BEST PERFORMERS:\n"
#         summary_text += f"Load Balance: {best_balance}\n"
#         summary_text += f"Spatial Locality: {best_locality}\n"
#         summary_text += f"Communication: {best_comm}"
        
#         ax9.text(0.05, 0.95, summary_text, transform=ax9.transAxes, 
#                 fontsize=10, verticalalignment='top', fontfamily='monospace',
#                 bbox=dict(boxstyle="round,pad=0.3", facecolor="lightgray", alpha=0.5))
        
#         plt.tight_layout()
#         plt.savefig('load_balance_comprehensive_analysis.png', dpi=300, bbox_inches='tight')
#         print("Saved comprehensive analysis to 'load_balance_comprehensive_analysis.png'")
#         plt.show()
        
#         return fig

# def main():
#     """Main execution function"""
#     print("Load Balance Analysis")
#     print("=" * 50)
    
#     analyzer = LoadBalanceAnalyzer()
#     analyzer.create_comprehensive_plot()
    
#     print("\nAnalysis complete!")
#     print("Check 'load_balance_comprehensive_analysis.png' for results")

# if __name__ == "__main__":
#     main()

# import os
# from collections import defaultdict
# from itertools import product
# import numpy as np
# import matplotlib.pyplot as plt
# import matplotlib as mpl
# from matplotlib.colors import ListedColormap, BoundaryNorm
# from mpl_toolkits.axes_grid1 import make_axes_locatable

# class SimData:
#     """
#     Structure for easy access to load costs reduced diagnostics
#     """

#     def __init__(self, filename, iterations, is_3D=False):
#         self.filename = filename
#         self.iterations = iterations
#         self.is_3D = is_3D
#         self.rank_arr = None
#         self.cost_arr = None
#         self.ranks = None
#         self.lb_efficiencies = None
#         self.lb_efficiency = None
#         self.lb_efficiency_max = None
#         self.lb_efficiency_min = None
#         self._load_data()

#     def __call__(self, iteration):
#         # For compatibility; does nothing as we load all data at once
#         pass

#     def _load_data(self):
#         try:
#             # Skip header lines starting with #
#             raw_data = np.loadtxt(self.filename, comments='#')
#             if self.is_3D:
#                 # box_id, rank, x, y, z, cost
#                 box_ids = raw_data[:, 0].astype(int)
#                 ranks = raw_data[:, 1].astype(int)
#                 x_coords = raw_data[:, 2].astype(int)
#                 y_coords = raw_data[:, 3].astype(int)
#                 z_coords = raw_data[:, 4].astype(int)
#                 costs = raw_data[:, 5]
#                 nx = x_coords.max() + 1
#                 ny = y_coords.max() + 1
#                 nz = z_coords.max() + 1
#                 self.rank_arr = np.full((nz, ny, nx), -1)
#                 self.cost_arr = np.zeros((nz, ny, nx))
#                 for i in range(len(box_ids)):
#                     self.rank_arr[z_coords[i], y_coords[i], x_coords[i]] = ranks[i]
#                     self.cost_arr[z_coords[i], y_coords[i], x_coords[i]] = costs[i]
#             else:
#                 # box_id, rank, x, y, cost
#                 box_ids = raw_data[:, 0].astype(int)
#                 ranks = raw_data[:, 1].astype(int)
#                 x_coords = raw_data[:, 2].astype(int)
#                 y_coords = raw_data[:, 3].astype(int)
#                 costs = raw_data[:, 4]
#                 nx = x_coords.max() + 1
#                 ny = y_coords.max() + 1
#                 self.rank_arr = np.full((ny, nx), -1)
#                 self.cost_arr = np.zeros((ny, nx))
#                 for i in range(len(box_ids)):
#                     self.rank_arr[y_coords[i], x_coords[i]] = ranks[i]
#                     self.cost_arr[y_coords[i], x_coords[i]] = costs[i]
#             # Compute load balance efficiency
#             rank_to_cost_map = {r: 0.0 for r in np.unique(ranks)}
#             for c, r in zip(costs, ranks):
#                 rank_to_cost_map[r] += c
#             efficiencies = np.array(list(rank_to_cost_map.values()))
#             efficiencies /= efficiencies.max() if efficiencies.max() > 0 else 1.0
#             self.ranks = np.array(list(rank_to_cost_map.keys()))
#             self.lb_efficiencies = efficiencies
#             self.lb_efficiency = efficiencies.mean()
#             self.lb_efficiency_max = efficiencies.max()
#             self.lb_efficiency_min = efficiencies.min()
#         except Exception as e:
#             print(f"Error loading data from {self.filename}: {str(e)}")
#             self.rank_arr = None
#             self.cost_arr = None

# def plot_2D_distribution(sim, title=""):
#     """Plot 2D distribution mapping"""
#     if not hasattr(sim, 'rank_arr') or sim.rank_arr is None:
#         print(f"No data available for {title}")
#         return None
    
#     # Get unique ranks (excluding -1)
#     unique_ranks = np.unique(sim.rank_arr[sim.rank_arr >= 0])
#     n_ranks = len(unique_ranks)
#     print(f"\nPlotting {title}")
#     print(f"Number of unique ranks: {n_ranks}")
    
#     # Create colormap with n_ranks + 1 colors (including one for unassigned cells)
#     colors = plt.cm.tab20(np.linspace(0, 1, n_ranks))
#     colors = np.vstack(([0.8, 0.8, 0.8, 1], colors))  # Add gray for unassigned
#     cmap = ListedColormap(colors)
    
#     # Create normalized rank array (0 for unassigned, 1+ for ranks)
#     rank_normalized = np.zeros_like(sim.rank_arr)
#     for i, rank in enumerate(unique_ranks):
#         rank_normalized[sim.rank_arr == rank] = i + 1
    
#     # Create plot
#     bounds = np.arange(-0.5, n_ranks + 1.5)
#     norm = BoundaryNorm(bounds, cmap.N)
#     im = plt.pcolormesh(rank_normalized, cmap=cmap, norm=norm)
    
#     # Add grid lines
#     plt.ylabel('j')
#     plt.xlabel('i')
#     plt.title(title)
#     plt.grid(True, which='both', color='gray', linewidth=0.1)
#     plt.gca().set_aspect('equal')
    
#     # Add rank labels
#     for j in range(sim.rank_arr.shape[0]):
#         for i in range(sim.rank_arr.shape[1]):
#             if sim.rank_arr[j, i] >= 0:  # Only label assigned cells
#                 plt.text(i + 0.5, j + 0.5, int(sim.rank_arr[j, i]),
#                         ha="center", va="center", color="black", fontsize=6)
    
#     # Add colorbar
#     divider = make_axes_locatable(plt.gca())
#     cax = divider.new_horizontal(size="5%", pad=0.05)
#     plt.gcf().add_axes(cax)
    
#     # Create custom colorbar labels
#     colorbar_labels = ['Unassigned'] + [f'Rank {int(r)}' for r in unique_ranks]
#     cb = plt.colorbar(im, cax=cax, orientation="vertical", 
#                      ticks=np.arange(len(colorbar_labels)))
#     cb.ax.set_yticklabels(colorbar_labels)
    
#     return im

# def plot_costs(sim, title=""):
#     """Plot computational costs"""
#     if not hasattr(sim, 'cost_arr') or sim.cost_arr is None:
#         print(f"No cost data available for {title}")
#         return
    
#     print(f"\nPlotting costs for {title}")
#     print(f"Cost range: {sim.cost_arr.min():.2e} to {sim.cost_arr.max():.2e}")
    
#     # Create mask for unassigned cells
#     mask = sim.rank_arr < 0
#     costs_masked = np.ma.array(sim.cost_arr, mask=mask)
    
#     # Plot with logarithmic scale
#     im = plt.pcolormesh(costs_masked, norm=mpl.colors.LogNorm())
#     plt.ylabel('j')
#     plt.xlabel('i')
#     plt.title(title)
#     plt.grid(True, which='both', color='gray', linewidth=0.1)
#     plt.gca().set_aspect('equal')
#     plt.colorbar(im, label='Cost (log scale)')

# def plot_3D_distribution(sim, title="", savefile=None):
#     """Plot 3D distribution mapping as a scatter plot colored by rank."""
#     if not hasattr(sim, 'rank_arr') or sim.rank_arr is None:
#         print(f"No data available for {title}")
#         return None
#     nz, ny, nx = sim.rank_arr.shape
#     z, y, x = np.indices((nz, ny, nx))
#     ranks = sim.rank_arr.flatten()
#     mask = ranks >= 0
#     x = x.flatten()[mask]
#     y = y.flatten()[mask]
#     z = z.flatten()[mask]
#     ranks = ranks[mask]
#     fig = plt.figure(figsize=(10, 8))
#     ax = fig.add_subplot(111, projection='3d')
#     p = ax.scatter(x, y, z, c=ranks, cmap='tab20', s=2)
#     ax.set_xlabel('i')
#     ax.set_ylabel('j')
#     ax.set_zlabel('k')
#     ax.set_title(title)
#     fig.colorbar(p, ax=ax, label='Rank')
#     plt.tight_layout()
#     if savefile:
#         plt.savefig(savefile, dpi=300, bbox_inches='tight')
#         print(f"Saved 3D distribution plot to {savefile}")
#     else:
#         plt.show()
#     plt.close(fig)

# def main():
#     print("Starting visualization...")
    
#     # Set is_3D to True for 3D data, False for 2D data
#     is_3D = True  # Change to False if using 2D files
    
#     # Load data for different algorithms
#     sim_knapsack = SimData('LBC_knapsack.txt', [1], is_3D=is_3D)
#     sim_sfc = SimData('LBC_sfc.txt', [1], is_3D=is_3D)
#     sim_hilbert = SimData('LBC_hilbert.txt', [1], is_3D=is_3D)
    
#     # Load the first iteration for each
#     for sim in [sim_knapsack, sim_sfc, sim_hilbert]:
#         sim(1)
    
#     # Compare distributions
#     if hasattr(sim_knapsack, 'rank_arr') and hasattr(sim_sfc, 'rank_arr') and hasattr(sim_hilbert, 'rank_arr'):
#         print("\nComparing distributions:")
#         print("Knapsack vs SFC identical:", np.array_equal(sim_knapsack.rank_arr, sim_sfc.rank_arr))
#         print("Knapsack vs Hilbert identical:", np.array_equal(sim_knapsack.rank_arr, sim_hilbert.rank_arr))
#         print("SFC vs Hilbert identical:", np.array_equal(sim_sfc.rank_arr, sim_hilbert.rank_arr))
    
#     if is_3D:
#         # Plot 3D distribution mappings
#         plot_3D_distribution(sim_knapsack, "Knapsack 3D Distribution", savefile='knapsack_3d_distribution_x.png')
#         plot_3D_distribution(sim_sfc, "SFC 3D Distribution", savefile='sfc_3d_distribution_x.png')
#         plot_3D_distribution(sim_hilbert, "Hilbert SFC 3D Distribution", savefile='hilbert_3d_distribution_x.png')
#     else:
#         # Plot 2D distribution mappings
#         plt.figure(figsize=(20, 6))
#         plt.subplot(131)
#         plot_2D_distribution(sim_knapsack, "Knapsack Distribution")
#         plt.subplot(132)
#         plot_2D_distribution(sim_sfc, "SFC Distribution")
#         plt.subplot(133)
#         plot_2D_distribution(sim_hilbert, "Hilbert SFC Distribution")
#         plt.tight_layout()
#         plt.savefig('distribution_mapping_comparison.png', dpi=300, bbox_inches='tight')
#         print("\nSaved distribution mapping comparison to distribution_mapping_comparison.png")
    
#     # Plot costs (2D only, or add 3D cost visualization if desired)
#     if not is_3D:
#         plt.figure(figsize=(20, 6))
#         plt.subplot(131)
#         plot_costs(sim_knapsack, "Knapsack Costs")
#         plt.subplot(132)
#         plot_costs(sim_sfc, "SFC Costs")
#         plt.subplot(133)
#         plot_costs(sim_hilbert, "Hilbert SFC Costs")
#         plt.tight_layout()
#         plt.savefig('cost_comparison.png', dpi=300, bbox_inches='tight')
#         print("Saved cost comparison to cost_comparison.png")
    
#     print("Knapsack LB efficiency:", sim_knapsack.lb_efficiency)
#     print("SFC LB efficiency:", sim_sfc.lb_efficiency)
#     print("Hilbert LB efficiency:", sim_hilbert.lb_efficiency)

# if __name__ == "__main__":
#     main()
