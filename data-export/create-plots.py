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

# -------------------------
# First plot: individual CSV files
# -------------------------
csv_files = ["C-Cells.csv", "E-Cells.csv", "H-Cells.csv", "N-Cells.csv", "P-Cells.csv", "R-Cells.csv"]

plt.figure(figsize=(10, 6))

for csv_file in csv_files:
    file_path = os.path.join(folder, csv_file)
    df = pd.read_csv(file_path, header=None, names=["Days", "Cells Population"])
    key = csv_file.split("-")[0]  # e.g., "C" from "C-Cells.csv"
    plt.plot(df["Days"], df["Cells Population"], label=key, color=colors[key])

plt.xlabel("Days")
plt.ylabel("Cells Population (log scale)")
plt.title("Cells Population Growth Over Time (from separate CSVs)")
plt.yscale("log")
plt.legend()
plt.grid(True, which="both", linestyle="--", linewidth=0.5)
plt.tight_layout()
plt.show()

# -------------------------
# Second plot: populations.csv
# -------------------------
pop_file = os.path.join(folder, "populations.csv")
df_pop = pd.read_csv(pop_file)

plt.figure(figsize=(10, 6))

for col in ["C", "P", "E", "N", "H", "R"]:
    plt.plot(df_pop["days"], df_pop[col], label=col, color=colors[col])

plt.xlabel("Days")
plt.ylabel("Cells Population (log scale)")
plt.title("Cells Population Growth Over Time (from populations.csv)")
plt.yscale("log")
plt.legend()
plt.grid(True, which="both", linestyle="--", linewidth=0.5)
plt.tight_layout()
plt.show()
