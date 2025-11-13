import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import seaborn as sns
import numpy as np
import argparse
import os

# Define a color palette inspired by the R script but with better contrast/aesthetics
# Using a dictionary for easy lookup
CATEGORY_COLORS = {
    "Metabolic": "#D55E00", 
    "Inflammatory": "#0072B2"
}

# Define markers for subtypes
SUBTYPE_MARKERS = {
    "Lung adenocarcinoma": "o",  # Circle
    "Squamous cell lung cancer": "s"  # Square
}

def format_pvalue(p):
    """Formats a p-value for display."""
    if pd.isna(p):
        return "N/A"
    if p < 0.001:
        return "<0.001"
    return f"{p:.3f}"

def parse_arguments():
    """Parses command-line arguments."""
    parser = argparse.ArgumentParser(description="Generate a publication-quality subtype-specific forest plot from MR results.")
    parser.add_argument(
        "--input",
        type=str,
        required=True,
        help="Path to the input CSV file (e.g., 'step05_mr_results_summary.csv')."
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        default="results/figures/python_generated",
        help="Directory to save the output figures."
    )
    parser.add_argument(
        "--output_name",
        type=str,
        default="FigureS2_SubtypeSpecific_ForestPlot_Optimized",
        help="Base name for the output PDF and PNG files."
    )
    return parser.parse_args()

def prepare_plot_data(df):
    """
    Filters and prepares the data for plotting, identifying significant exposures
    and structuring the data for the forest plot.
    """
    # Identify significant exposures based on any of the three outcomes
    # An exposure is significant if it has a significant FDR or nominal p-value for any subtype
    significant_exposures = df[
        (df['significant_fdr'] == True) | (df['pval'] < 0.05)
    ]['exposure'].unique()

    # Filter for the significant exposures and the two main subtypes
    plot_df = df[
        df['exposure'].isin(significant_exposures) &
        df['outcome'].isin(['Lung adenocarcinoma', 'Squamous cell lung cancer'])
    ].copy()

    if plot_df.empty:
        print("Warning: No significant data found to plot.")
        return pd.DataFrame()

    # Get the overall lung cancer p-values for sorting
    overall_lc_pvals = df[df['outcome'] == 'Lung cancer'].set_index('exposure')['pval']
    plot_df['overall_lc_pval'] = plot_df['exposure'].map(overall_lc_pvals)

    # Sort exposures: by category, then by overall lung cancer p-value
    # This mimics the sorting from the R script to group related exposures
    plot_df = plot_df.sort_values(
        by=['category', 'overall_lc_pval'],
        ascending=[True, True]
    )
    
    # Create a stable categorical order for exposures for plotting
    exposure_order = plot_df['exposure_label'].unique()
    plot_df['exposure_label_cat'] = pd.Categorical(plot_df['exposure_label'], categories=exposure_order, ordered=True)

    return plot_df

def create_forest_plot(df, output_dir, output_name):
    """
    Generates and saves a publication-quality forest plot using matplotlib and seaborn.
    """
    if df.empty:
        return

    # Set plot style
    sns.set_style("whitegrid")
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['axes.labelweight'] = 'bold'

    # Determine figure size based on number of exposures
    num_exposures = len(df['exposure_label'].unique())
    fig_height = max(6, num_exposures * 0.5)
    fig_width = 10

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    # Get the ordered list of exposures for the y-axis
    y_labels = df['exposure_label_cat'].cat.categories
    y_ticks = np.arange(len(y_labels))

    # Create a mapping from label to y-position
    y_map = {label: i for i, label in enumerate(y_labels)}

    # Define the vertical offset for dodging
    dodge = 0.2

    # Plot each point and its error bar
    for _, row in df.iterrows():
        y_pos = y_map[row['exposure_label']]
        subtype = row['outcome']
        
        # Determine vertical position with dodging
        y_plot = y_pos + dodge if subtype == 'Lung adenocarcinoma' else y_pos - dodge

        # Determine color by category
        color = CATEGORY_COLORS.get(row['category'], 'grey')
        
        # Plot error bar
        ax.errorbar(
            x=[row['or_lci'], row['or_uci']],
            y=[y_plot, y_plot],
            color=color,
            linewidth=1,
            capsize=3
        )
        
        # Plot point estimate
        # Use filled markers for FDR significant, open for nominal, smaller for non-sig
        marker = SUBTYPE_MARKERS[subtype]
        if row['significant_fdr']:
            ax.plot(row['or'], y_plot, marker=marker, color=color, markersize=7, linestyle='None')
        elif row['pval'] < 0.05:
            ax.plot(row['or'], y_plot, marker=marker, color=color, markersize=7, linestyle='None', mfc='white')
        else:
            # Optionally, don't plot non-significant points or make them very subtle
            ax.plot(row['or'], y_plot, marker=marker, color='lightgrey', markersize=5, linestyle='None')


    # --- Plot Aesthetics ---
    
    # Vertical line at OR=1
    ax.axvline(x=1, linestyle='--', color='grey', linewidth=0.8)

    # X-axis (log scale)
    ax.set_xscale('log')
    ax.xaxis.set_major_formatter(mticker.FormatStrFormatter('%.2f'))
    x_ticks = [0.7, 0.8, 1.0, 1.2, 1.5]
    ax.set_xticks(x_ticks)
    ax.set_xlabel("Odds Ratio (95% CI)")
    
    # Y-axis
    ax.set_yticks(y_ticks)
    ax.set_yticklabels(y_labels, fontsize=10)
    ax.set_ylabel("")
    ax.invert_yaxis() # To match typical forest plot order

    # Title and Subtitle
    fig.suptitle(
        "Subtype-Specific Associations with Lung Cancer",
        fontsize=16,
        fontweight='bold',
        y=0.98
    )
    ax.set_title(
        "Mendelian Randomization Analysis of Significant Exposures",
        fontsize=12,
        color='dimgray',
        pad=20
    )

    # Grid and Spines
    ax.grid(axis='x', linestyle=':', color='lightgrey')
    ax.grid(axis='y', which='major', linestyle='',) # No horizontal grid lines
    sns.despine(left=True, bottom=False) # Remove left and top spines

    # --- Legend ---
    legend_elements = []
    # Category legend
    for category, color in CATEGORY_COLORS.items():
        legend_elements.append(plt.Line2D([0], [0], color=color, lw=4, label=category))
    
    # Subtype legend
    for subtype, marker in SUBTYPE_MARKERS.items():
        label = "Adenocarcinoma" if "adeno" in subtype else "Squamous"
        legend_elements.append(plt.Line2D([0], [0], marker=marker, color='grey', label=label, linestyle='None', markersize=8))

    # Significance legend
    legend_elements.append(plt.Line2D([0], [0], marker='o', color='black', label='FDR < 0.05', linestyle='None', markersize=8))
    legend_elements.append(plt.Line2D([0], [0], marker='o', color='black', mfc='white', label='P < 0.05', linestyle='None', markersize=8))

    ax.legend(handles=legend_elements, loc='lower center', bbox_to_anchor=(0.5, -0.15), ncol=4, frameon=False)

    plt.tight_layout(rect=[0, 0.05, 1, 0.95]) # Adjust layout to make space for legend and title

    # Save figure
    os.makedirs(output_dir, exist_ok=True)
    pdf_path = os.path.join(output_dir, f"{output_name}.pdf")
    png_path = os.path.join(output_dir, f"{output_name}.png")
    
    fig.savefig(pdf_path, dpi=300, bbox_inches='tight')
    fig.savefig(png_path, dpi=300, bbox_inches='tight')
    
    print(f"✓ Optimized forest plot saved to:")
    print(f"  - PDF: {pdf_path}")
    print(f"  - PNG: {png_path}")


def main():
    """Main function to run the script."""
    # In a real scenario, we would use parse_arguments()
    # args = parse_arguments()
    # For this environment, let's hardcode the paths.
    
    input_file = "c:\\Users\\surob\\Documents\\GWAS\\1.主MR分析\\results\\data\\step05_mr_results_summary.csv"
    output_dir = "c:\\Users\\surob\\Documents\\GWAS\\1.主MR分析\\results\\figures\\python_generated"
    output_name = "FigureS2_SubtypeSpecific_ForestPlot_Optimized"

    if not os.path.exists(input_file):
        print(f"Error: Input file not found at {input_file}")
        return

    try:
        df = pd.read_csv(input_file)
        plot_data = prepare_plot_data(df)
        create_forest_plot(plot_data, output_dir, output_name)
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    main()