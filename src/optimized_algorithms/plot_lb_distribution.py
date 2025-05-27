import os
from collections import defaultdict
from itertools import product
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import ListedColormap, BoundaryNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable

class SimData:
    """
    Structure for easy access to load costs reduced diagnostics
    """

    def __init__(self, filename, iterations, is_3D=False):
        self.filename = filename
        self.iterations = iterations
        self.is_3D = is_3D
        self.rank_arr = None
        self.cost_arr = None
        self.ranks = None
        self.lb_efficiencies = None
        self.lb_efficiency = None
        self.lb_efficiency_max = None
        self.lb_efficiency_min = None
        self._load_data()

    def __call__(self, iteration):
        # For compatibility; does nothing as we load all data at once
        pass

    def _load_data(self):
        try:
            # Skip header lines starting with #
            raw_data = np.loadtxt(self.filename, comments='#')
            if self.is_3D:
                # box_id, rank, x, y, z, cost
                box_ids = raw_data[:, 0].astype(int)
                ranks = raw_data[:, 1].astype(int)
                x_coords = raw_data[:, 2].astype(int)
                y_coords = raw_data[:, 3].astype(int)
                z_coords = raw_data[:, 4].astype(int)
                costs = raw_data[:, 5]
                nx = x_coords.max() + 1
                ny = y_coords.max() + 1
                nz = z_coords.max() + 1
                self.rank_arr = np.full((nz, ny, nx), -1)
                self.cost_arr = np.zeros((nz, ny, nx))
                for i in range(len(box_ids)):
                    self.rank_arr[z_coords[i], y_coords[i], x_coords[i]] = ranks[i]
                    self.cost_arr[z_coords[i], y_coords[i], x_coords[i]] = costs[i]
            else:
                # box_id, rank, x, y, cost
                box_ids = raw_data[:, 0].astype(int)
                ranks = raw_data[:, 1].astype(int)
                x_coords = raw_data[:, 2].astype(int)
                y_coords = raw_data[:, 3].astype(int)
                costs = raw_data[:, 4]
                nx = x_coords.max() + 1
                ny = y_coords.max() + 1
                self.rank_arr = np.full((ny, nx), -1)
                self.cost_arr = np.zeros((ny, nx))
                for i in range(len(box_ids)):
                    self.rank_arr[y_coords[i], x_coords[i]] = ranks[i]
                    self.cost_arr[y_coords[i], x_coords[i]] = costs[i]
            # Compute load balance efficiency
            rank_to_cost_map = {r: 0.0 for r in np.unique(ranks)}
            for c, r in zip(costs, ranks):
                rank_to_cost_map[r] += c
            efficiencies = np.array(list(rank_to_cost_map.values()))
            efficiencies /= efficiencies.max() if efficiencies.max() > 0 else 1.0
            self.ranks = np.array(list(rank_to_cost_map.keys()))
            self.lb_efficiencies = efficiencies
            self.lb_efficiency = efficiencies.mean()
            self.lb_efficiency_max = efficiencies.max()
            self.lb_efficiency_min = efficiencies.min()
        except Exception as e:
            print(f"Error loading data from {self.filename}: {str(e)}")
            self.rank_arr = None
            self.cost_arr = None

def plot_2D_distribution(sim, title=""):
    """Plot 2D distribution mapping"""
    if not hasattr(sim, 'rank_arr') or sim.rank_arr is None:
        print(f"No data available for {title}")
        return None
    
    # Get unique ranks (excluding -1)
    unique_ranks = np.unique(sim.rank_arr[sim.rank_arr >= 0])
    n_ranks = len(unique_ranks)
    print(f"\nPlotting {title}")
    print(f"Number of unique ranks: {n_ranks}")
    
    # Create colormap with n_ranks + 1 colors (including one for unassigned cells)
    colors = plt.cm.tab20(np.linspace(0, 1, n_ranks))
    colors = np.vstack(([0.8, 0.8, 0.8, 1], colors))  # Add gray for unassigned
    cmap = ListedColormap(colors)
    
    # Create normalized rank array (0 for unassigned, 1+ for ranks)
    rank_normalized = np.zeros_like(sim.rank_arr)
    for i, rank in enumerate(unique_ranks):
        rank_normalized[sim.rank_arr == rank] = i + 1
    
    # Create plot
    bounds = np.arange(-0.5, n_ranks + 1.5)
    norm = BoundaryNorm(bounds, cmap.N)
    im = plt.pcolormesh(rank_normalized, cmap=cmap, norm=norm)
    
    # Add grid lines
    plt.ylabel('j')
    plt.xlabel('i')
    plt.title(title)
    plt.grid(True, which='both', color='gray', linewidth=0.1)
    plt.gca().set_aspect('equal')
    
    # Add rank labels
    for j in range(sim.rank_arr.shape[0]):
        for i in range(sim.rank_arr.shape[1]):
            if sim.rank_arr[j, i] >= 0:  # Only label assigned cells
                plt.text(i + 0.5, j + 0.5, int(sim.rank_arr[j, i]),
                        ha="center", va="center", color="black", fontsize=6)
    
    # Add colorbar
    divider = make_axes_locatable(plt.gca())
    cax = divider.new_horizontal(size="5%", pad=0.05)
    plt.gcf().add_axes(cax)
    
    # Create custom colorbar labels
    colorbar_labels = ['Unassigned'] + [f'Rank {int(r)}' for r in unique_ranks]
    cb = plt.colorbar(im, cax=cax, orientation="vertical", 
                     ticks=np.arange(len(colorbar_labels)))
    cb.ax.set_yticklabels(colorbar_labels)
    
    return im

def plot_costs(sim, title=""):
    """Plot computational costs"""
    if not hasattr(sim, 'cost_arr') or sim.cost_arr is None:
        print(f"No cost data available for {title}")
        return
    
    print(f"\nPlotting costs for {title}")
    print(f"Cost range: {sim.cost_arr.min():.2e} to {sim.cost_arr.max():.2e}")
    
    # Create mask for unassigned cells
    mask = sim.rank_arr < 0
    costs_masked = np.ma.array(sim.cost_arr, mask=mask)
    
    # Plot with logarithmic scale
    im = plt.pcolormesh(costs_masked, norm=mpl.colors.LogNorm())
    plt.ylabel('j')
    plt.xlabel('i')
    plt.title(title)
    plt.grid(True, which='both', color='gray', linewidth=0.1)
    plt.gca().set_aspect('equal')
    plt.colorbar(im, label='Cost (log scale)')

def plot_3D_distribution(sim, title="", savefile=None):
    """Plot 3D distribution mapping as a scatter plot colored by rank."""
    if not hasattr(sim, 'rank_arr') or sim.rank_arr is None:
        print(f"No data available for {title}")
        return None
    nz, ny, nx = sim.rank_arr.shape
    z, y, x = np.indices((nz, ny, nx))
    ranks = sim.rank_arr.flatten()
    mask = ranks >= 0
    x = x.flatten()[mask]
    y = y.flatten()[mask]
    z = z.flatten()[mask]
    ranks = ranks[mask]
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    p = ax.scatter(x, y, z, c=ranks, cmap='tab20', s=2)
    ax.set_xlabel('i')
    ax.set_ylabel('j')
    ax.set_zlabel('k')
    ax.set_title(title)
    fig.colorbar(p, ax=ax, label='Rank')
    plt.tight_layout()
    if savefile:
        plt.savefig(savefile, dpi=300, bbox_inches='tight')
        print(f"Saved 3D distribution plot to {savefile}")
    else:
        plt.show()
    plt.close(fig)

def main():
    print("Starting visualization...")
    
    # Set is_3D to True for 3D data, False for 2D data
    is_3D = True  # Change to False if using 2D files
    
    # Load data for different algorithms
    sim_knapsack = SimData('LBC_knapsack.txt', [1], is_3D=is_3D)
    sim_sfc = SimData('LBC_sfc.txt', [1], is_3D=is_3D)
    sim_hilbert = SimData('LBC_hilbert.txt', [1], is_3D=is_3D)
    
    # Load the first iteration for each
    for sim in [sim_knapsack, sim_sfc, sim_hilbert]:
        sim(1)
    
    # Compare distributions
    if hasattr(sim_knapsack, 'rank_arr') and hasattr(sim_sfc, 'rank_arr') and hasattr(sim_hilbert, 'rank_arr'):
        print("\nComparing distributions:")
        print("Knapsack vs SFC identical:", np.array_equal(sim_knapsack.rank_arr, sim_sfc.rank_arr))
        print("Knapsack vs Hilbert identical:", np.array_equal(sim_knapsack.rank_arr, sim_hilbert.rank_arr))
        print("SFC vs Hilbert identical:", np.array_equal(sim_sfc.rank_arr, sim_hilbert.rank_arr))
    
    if is_3D:
        # Plot 3D distribution mappings
        plot_3D_distribution(sim_knapsack, "Knapsack 3D Distribution", savefile='knapsack_3d_distribution.png')
        plot_3D_distribution(sim_sfc, "SFC 3D Distribution", savefile='sfc_3d_distribution.png')
        plot_3D_distribution(sim_hilbert, "Hilbert SFC 3D Distribution", savefile='hilbert_3d_distribution.png')
    else:
        # Plot 2D distribution mappings
        plt.figure(figsize=(20, 6))
        plt.subplot(131)
        plot_2D_distribution(sim_knapsack, "Knapsack Distribution")
        plt.subplot(132)
        plot_2D_distribution(sim_sfc, "SFC Distribution")
        plt.subplot(133)
        plot_2D_distribution(sim_hilbert, "Hilbert SFC Distribution")
        plt.tight_layout()
        plt.savefig('distribution_mapping_comparison.png', dpi=300, bbox_inches='tight')
        print("\nSaved distribution mapping comparison to distribution_mapping_comparison.png")
    
    # Plot costs (2D only, or add 3D cost visualization if desired)
    if not is_3D:
        plt.figure(figsize=(20, 6))
        plt.subplot(131)
        plot_costs(sim_knapsack, "Knapsack Costs")
        plt.subplot(132)
        plot_costs(sim_sfc, "SFC Costs")
        plt.subplot(133)
        plot_costs(sim_hilbert, "Hilbert SFC Costs")
        plt.tight_layout()
        plt.savefig('cost_comparison.png', dpi=300, bbox_inches='tight')
        print("Saved cost comparison to cost_comparison.png")
    
    print("Knapsack LB efficiency:", sim_knapsack.lb_efficiency)
    print("SFC LB efficiency:", sim_sfc.lb_efficiency)
    print("Hilbert LB efficiency:", sim_hilbert.lb_efficiency)

if __name__ == "__main__":
    main()