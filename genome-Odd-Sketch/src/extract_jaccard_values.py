#!/usr/bin/env python3
import sys
import os

def extract_jaccard_values(input_file, output_file):
    """
    Extract Jaccard coefficients from bindash output and create a summary.
    
    Bindash output format:
    query_genome  target_genome  mutation_distance  p_value  jaccard_index
    """
    jaccard_values = []
    pair_count = 0
    
    with open(input_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            
            parts = line.split('\t')
            if len(parts) >= 5:
                query = parts[0]
                target = parts[1]
                mutation_dist = parts[2]
                p_value = parts[3]
                jaccard_str = parts[4]
                
                # Parse jaccard index (format: numerator/denominator)
                if '/' in jaccard_str:
                    numerator, denominator = jaccard_str.split('/')
                    jaccard_index = float(numerator) / float(denominator)
                else:
                    jaccard_index = float(jaccard_str)
                
                jaccard_values.append({
                    'query': query,
                    'target': target,
                    'mutation_dist': float(mutation_dist),
                    'p_value': float(p_value),
                    'jaccard': jaccard_index
                })
                pair_count += 1
    
    # Write summary to output file
    with open(output_file, 'w') as f:
        f.write("# BinDash Jaccard Coefficient Results Summary\n")
        f.write(f"# Total genome pairs analyzed: {pair_count}\n")
        f.write(f"# Average Jaccard coefficient: {sum(v['jaccard'] for v in jaccard_values) / len(jaccard_values):.6f}\n")
        f.write(f"# Min Jaccard coefficient: {min(v['jaccard'] for v in jaccard_values):.6f}\n")
        f.write(f"# Max Jaccard coefficient: {max(v['jaccard'] for v in jaccard_values):.6f}\n")
        f.write("\n")
        f.write("query_genome\ttarget_genome\tmutation_distance\tp_value\tjaccard_coefficient\n")
        
        for v in jaccard_values:
            f.write(f"{v['query']}\t{v['target']}\t{v['mutation_dist']:.6f}\t{v['p_value']:.6e}\t{v['jaccard']:.6f}\n")
    
    return len(jaccard_values)

if __name__ == "__main__":
    input_file = "bindash_jaccard_results.txt"
    output_file = "bindash_jaccard_summary.txt"
    
    if os.path.exists(input_file):
        count = extract_jaccard_values(input_file, output_file)
        print(f"Processed {count} genome pairs")
        print(f"Results saved to: {output_file}")
    else:
        print(f"Input file {input_file} not found")