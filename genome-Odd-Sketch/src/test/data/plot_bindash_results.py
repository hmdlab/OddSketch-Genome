#!/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

def plot_bindash_results(input_file):
    """
    Create plots for BinDash Jaccard coefficient results.
    """
    # Read the data
    data = []
    with open(input_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('#') or not line or line.startswith('pair_id'):
                continue
            
            parts = line.split('\t')
            if len(parts) >= 6:
                pair_id = int(parts[0])
                mutation_dist = float(parts[3])
                p_value = float(parts[4])
                jaccard = float(parts[5])
                
                data.append({
                    'pair_id': pair_id,
                    'mutation_distance': mutation_dist,
                    'p_value': p_value,
                    'jaccard': jaccard
                })
    
    df = pd.DataFrame(data)
    
    # Set up the plotting style
    plt.style.use('default')
    sns.set_palette("husl")
    
    # Create a figure with multiple subplots
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    fig.suptitle('BinDash Jaccard Coefficient Analysis for 500 Genome Pairs', fontsize=16, fontweight='bold')
    
    # 1. Histogram of Jaccard coefficients
    axes[0, 0].hist(df['jaccard'], bins=30, alpha=0.7, color='skyblue', edgecolor='black')
    axes[0, 0].set_xlabel('Jaccard Coefficient')
    axes[0, 0].set_ylabel('Frequency')
    axes[0, 0].set_title('Distribution of Jaccard Coefficients')
    axes[0, 0].grid(True, alpha=0.3)
    
    # Add statistics text
    mean_jaccard = df['jaccard'].mean()
    std_jaccard = df['jaccard'].std()
    axes[0, 0].axvline(mean_jaccard, color='red', linestyle='--', linewidth=2, label=f'Mean: {mean_jaccard:.3f}')
    axes[0, 0].legend()
    
    # 2. Jaccard vs Mutation Distance
    axes[0, 1].scatter(df['mutation_distance'], df['jaccard'], alpha=0.6, s=30, color='coral')
    axes[0, 1].set_xlabel('Mutation Distance')
    axes[0, 1].set_ylabel('Jaccard Coefficient')
    axes[0, 1].set_title('Jaccard Coefficient vs Mutation Distance')
    axes[0, 1].grid(True, alpha=0.3)
    
    # Add trend line
    z = np.polyfit(df['mutation_distance'], df['jaccard'], 1)
    p = np.poly1d(z)
    axes[0, 1].plot(df['mutation_distance'], p(df['mutation_distance']), "r--", alpha=0.8, linewidth=2)
    
    # Calculate correlation
    correlation = df['mutation_distance'].corr(df['jaccard'])
    axes[0, 1].text(0.05, 0.95, f'Correlation: {correlation:.3f}', transform=axes[0, 1].transAxes, 
                   bbox=dict(boxstyle="round", facecolor='white', alpha=0.8))
    
    # 3. Jaccard coefficients by pair ID (time series like)
    axes[1, 0].plot(df['pair_id'], df['jaccard'], linewidth=1, alpha=0.7, color='green')
    axes[1, 0].scatter(df['pair_id'], df['jaccard'], s=20, alpha=0.6, color='darkgreen')
    axes[1, 0].set_xlabel('Pair ID')
    axes[1, 0].set_ylabel('Jaccard Coefficient')
    axes[1, 0].set_title('Jaccard Coefficients Across Genome Pairs')
    axes[1, 0].grid(True, alpha=0.3)
    
    # 4. Box plot and summary statistics
    axes[1, 1].boxplot(df['jaccard'], patch_artist=True, 
                      boxprops=dict(facecolor='lightblue', alpha=0.7),
                      medianprops=dict(color='red', linewidth=2))
    axes[1, 1].set_ylabel('Jaccard Coefficient')
    axes[1, 1].set_title('Jaccard Coefficient Distribution Summary')
    axes[1, 1].grid(True, alpha=0.3)
    
    # Add summary statistics text
    stats_text = f"""Statistics:
Mean: {df['jaccard'].mean():.4f}
Median: {df['jaccard'].median():.4f}
Std: {df['jaccard'].std():.4f}
Min: {df['jaccard'].min():.4f}
Max: {df['jaccard'].max():.4f}
Q1: {df['jaccard'].quantile(0.25):.4f}
Q3: {df['jaccard'].quantile(0.75):.4f}"""
    
    axes[1, 1].text(1.1, 0.5, stats_text, transform=axes[1, 1].transAxes, 
                   bbox=dict(boxstyle="round", facecolor='lightyellow', alpha=0.8),
                   verticalalignment='center', fontsize=10)
    
    plt.tight_layout()
    
    # Save the plot
    output_file = 'bindash_jaccard_analysis.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Plot saved as: {output_file}")
    
    # Create a separate detailed histogram
    plt.figure(figsize=(10, 6))
    plt.hist(df['jaccard'], bins=50, alpha=0.7, color='steelblue', edgecolor='black')
    plt.xlabel('Jaccard Coefficient', fontsize=12)
    plt.ylabel('Frequency', fontsize=12)
    plt.title('Detailed Distribution of Jaccard Coefficients (BinDash Results)', fontsize=14, fontweight='bold')
    plt.grid(True, alpha=0.3)
    
    # Add vertical lines for statistics
    plt.axvline(df['jaccard'].mean(), color='red', linestyle='--', linewidth=2, label=f'Mean: {df["jaccard"].mean():.3f}')
    plt.axvline(df['jaccard'].median(), color='orange', linestyle='--', linewidth=2, label=f'Median: {df["jaccard"].median():.3f}')
    plt.legend()
    
    plt.tight_layout()
    plt.savefig('bindash_jaccard_histogram_detailed.png', dpi=300, bbox_inches='tight')
    print("Detailed histogram saved as: bindash_jaccard_histogram_detailed.png")
    
    # Print summary statistics
    print("\n=== BinDash Jaccard Coefficient Summary ===")
    print(f"Total pairs: {len(df)}")
    print(f"Mean Jaccard coefficient: {df['jaccard'].mean():.6f}")
    print(f"Standard deviation: {df['jaccard'].std():.6f}")
    print(f"Minimum: {df['jaccard'].min():.6f}")
    print(f"Maximum: {df['jaccard'].max():.6f}")
    print(f"Median: {df['jaccard'].median():.6f}")
    print(f"25th percentile: {df['jaccard'].quantile(0.25):.6f}")
    print(f"75th percentile: {df['jaccard'].quantile(0.75):.6f}")
    print(f"Correlation with mutation distance: {df['mutation_distance'].corr(df['jaccard']):.6f}")

if __name__ == "__main__":
    input_file = "bindash_500pairs_jaccard.txt"
    
    try:
        plot_bindash_results(input_file)
        plt.show()
    except Exception as e:
        print(f"Error creating plots: {e}")
        print("Make sure matplotlib, pandas, and seaborn are installed:")
        print("pip install matplotlib pandas seaborn")