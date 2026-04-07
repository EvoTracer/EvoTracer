import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import numpy as np

def main():
    # File configuration
    input_file = 'Analysis_Large_Indels_1kb/Large_Indels_Details.csv'
    output_dir = 'Analysis_Large_Indels_1kb'
    output_img = os.path.join(output_dir, 'Indel_Length_Distribution.png')

    if not os.path.exists(input_file):
        print(f"Error: {input_file} not found.")
        return

    # 1. Read Data
    df = pd.read_csv(input_file)
    
    # Ensure Actual_Seq_Length is numeric
    df['Actual_Seq_Length'] = pd.to_numeric(df['Actual_Seq_Length'], errors='coerce')
    df = df.dropna(subset=['Actual_Seq_Length'])
    
    # Filter out indels smaller than 1KB
    df = df[df['Actual_Seq_Length'] >= 1000]

    lengths = df['Actual_Seq_Length']
    if len(lengths) == 0:
        print("No indels >= 1KB found.")
        return

    max_len = lengths.max()
    print(f"Max Indel Length: {max_len}")

    # Print top 6 largest indels
    print("\nTop 6 Largest Indels:")
    top_6 = df.nlargest(6, 'Actual_Seq_Length')
    # Print all columns
    print(top_6.to_string(index=False))
    print("\n")

    # 2. Define Bins (50KB steps)
    # The user requested intervals like <50KB, 50-100KB, etc.
    step = 50000
    # Determine the upper bound for bins
    upper_bound = int(np.ceil(max_len / step)) * step
    if upper_bound == 0: upper_bound = step # Handle case with no data or small data

    bins = range(0, upper_bound + step + 1, step)
    
    # 3. Create Labels
    labels = []
    for i in range(len(bins)-1):
        start = bins[i]
        end = bins[i+1]
        
        # Convert to KB
        start_kb = start // 1000
        end_kb = end // 1000
        
        if i == 0:
            labels.append(f"<{end_kb}KB")
        else:
            labels.append(f"{start_kb}-{end_kb}KB")

    # 4. Binning
    # right=True/False? usually [start, end). large indels usually discrete?
    # Let's use right=False: [0, 50000), [50000, 100000)
    df['Bin'] = pd.cut(df['Actual_Seq_Length'], bins=bins, labels=labels, right=False)
    
    # 5. Aggregate
    # Use value_counts and reindex to ensure all bins are present (or at least up to max)
    counts = df['Bin'].value_counts().sort_index()

    print("Counts per bin:")
    print(counts)

    # Convert to DataFrame for easier plotting with seaborn
    plot_df = counts.reset_index()
    plot_df.columns = ['Range', 'Count']

    # Filter out empty bins at the tail if desired, but here we show all up to max
    # If the tail is extremely sparse, the plot might look wide.
    # We will plot all existing bins.


    # 6. Plot
    plt.figure(figsize=(10, 6))
    sns.set_style("whitegrid")
    ax = sns.barplot(data=plot_df, x='Range', y='Count', color='steelblue', edgecolor='black')

    # Add text on top of bars
    for p in ax.patches:
        height = p.get_height()
        if height > 0:
            ax.text(p.get_x() + p.get_width()/2., height + 0.5,
                    f'{int(height)}', ha='center', va='bottom')

    plt.title('Distribution of Large Indel Lengths', fontsize=16)
    plt.xlabel('Length Range', fontsize=12)
    plt.ylabel('Count', fontsize=12)
    plt.xticks(rotation=45)
    plt.tight_layout()

    # Save
    plt.savefig(output_img)
    print(f"\nHistogram saved to: {output_img}")

    # Print top 6 largest indels AT THE END
    print("\n" + "="*50)
    print("TOP 6 LARGEST INDELS (>1KB)")
    print("="*50)
    top_6 = df.nlargest(6, 'Actual_Seq_Length')
    # Print all columns
    print(top_6.to_string(index=False))
    print("="*50 + "\n")


if __name__ == "__main__":
    main()
