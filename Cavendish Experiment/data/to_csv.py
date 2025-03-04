import os
import pandas as pd

# List of files to convert
file_paths = [
    "Clockwise_Feb7th.txt",
    "Counterclockwise_Feb7th.txt",
    "Empty_Feb7th.txt"
]

# Directory to save the CSV files
output_dir = "data_IN_csv"
os.makedirs(output_dir, exist_ok=True)


# Function to process and save files as CSV
def convert_txt_to_csv(file_path, output_dir):
    # Read the file into a dataframe, skipping the first two lines for the header
    df = pd.read_csv(file_path, delimiter='\t', skiprows=2,
                     names=["Time (s)", "Position (m)"])
    # Drop incomplete or invalid rows
    df.dropna(inplace=True)

    # Save to CSV
    output_file_name = os.path.basename(file_path).replace('.txt', '.csv')
    output_path = os.path.join(output_dir, output_file_name)
    df.to_csv(output_path, index=False)
    print(f"Saved: {output_path}")


# Convert each file
for file in file_paths:
    convert_txt_to_csv(file, output_dir)
