import os
import pandas as pd
import matplotlib.pyplot as plt

# Path to your folder
folder = os.path.expanduser("~/Documents/dev/ESB-Workshop-tutorial/pancreatic_tumor_model/data-export")

# Define colors for consistency
colors = {
    "C": "blue",
    "P": "red",
    "E": "orange",
    "N": "purple",
    "H": "green",
    "R": "cyan"
}


# Create subplots: 2 rows, 1 column
fig, axes = plt.subplots(2, 1, figsize=(10, 10), sharex=True)

# -------------------------
# First plot: individual CSV files
# -------------------------
csv_files = ["C-Cells.csv", "E-Cells.csv", "H-Cells.csv", "N-Cells.csv", "P-Cells.csv", "R-Cells.csv"]

for csv_file in csv_files:
    file_path = os.path.join(folder, csv_file)
    df = pd.read_csv(file_path, header=None, names=["Days", "Cells Population"])
    key = csv_file.split("-")[0]  # e.g., "C" from "C-Cells.csv"
    axes[0].plot(df["Days"], df["Cells Population"], label=key, color=colors[key])

axes[0].set_ylabel("Cells Population (log scale)")
axes[0].set_title("Cells Population Growth Over Time (Paper Actual Data)")
axes[0].set_yscale("log")
axes[0].legend()
axes[0].grid(True, which="both", linestyle="--", linewidth=0.5)

# -------------------------
# Second plot: populations.csv
# -------------------------
pop_file = os.path.join(folder, "populations.csv")
df_pop = pd.read_csv(pop_file)

for col in ["C", "P", "E", "N", "H", "R"]:
    axes[1].plot(df_pop["days"], df_pop[col], label=col, color=colors[col])

axes[1].set_xlabel("Days")
axes[1].set_ylabel("Cells Population (log scale)")
axes[1].set_title("Cells Population Growth Over Time (Simulation Data using BioDynaMo)")
axes[1].set_yscale("log")
axes[1].legend()
axes[1].grid(True, which="both", linestyle="--", linewidth=0.5)

# Adjust layout
plt.tight_layout()
plt.show()
