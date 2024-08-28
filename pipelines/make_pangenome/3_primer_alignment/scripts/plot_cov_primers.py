import os
import pandas as pd
import matplotlib.pyplot as plt
import sys

# Get the input and output file paths from the command line arguments
input_file_path, output_file_path = sys.argv[1], sys.argv[2]

# Extract the sample name from the input file name
sample_name = os.path.basename(input_file_path).split('.')[0]

# Initialize the plot
fig, ax = plt.subplots(figsize=(15,10))

# Check if the input file exists before trying to read it
if os.path.exists(input_file_path):
    df = pd.read_table(input_file_path, sep="\t", names=["Ref", "Pos", "Depth"])
    ax.fill_between(df["Pos"], df["Depth"], label=sample_name, color='blue', alpha=0.5)

    # Set the title and legend
    ax.set_title(f"Coverage Depth for {sample_name}")
    ax.legend()

    # Set the y-axis limits
    ax.set_ylim(0, df["Depth"].max())  # Adjust these values as needed

    # Save the plot
    plt.savefig(output_file_path)
else:
    print(f"Input file {input_file_path} does not exist.")