import os
import pandas as pd
import matplotlib.pyplot as plt

# Path to your folder
folder = os.path.expanduser(
    "~/Documents/dev/ESB-Workshop-tutorial/pancreatic_tumor_model/data-export"
)

# Fixed colors per population (matches your reference figure)
colors = {
    "C": "blue",
    "P": "red",
    "E": "orange",
    "N": "green",
    "H": "cyan",
    "R": "purple",
}

# Original and scaled CSVs
original_csvs = ["C-Cells.csv", "E-Cells.csv", "H-Cells.csv", "N-Cells.csv", "P-Cells.csv", "R-Cells.csv"]
scaled_csvs   = ["C-Cells_scaled_global.csv", "E-Cells_scaled_global.csv", "H-Cells_scaled_global.csv",
                 "N-Cells_scaled_global.csv", "P-Cells_scaled_global.csv", "R-Cells_scaled_global.csv"]

# Create subplots: 2 rows, 1 column
fig, axes = plt.subplots(2, 1, figsize=(11, 10), sharex=True)

# -------------------------
# Plot 1: separate CSVs (original + scaled)
# -------------------------
for csv_file in original_csvs:
    path = os.path.join(folder, csv_file)
    df = pd.read_csv(path, header=None, names=["Days", "Cells Population"])
    key = csv_file.split("-")[0]  # "C" from "C-Cells.csv"
    axes[0].plot(df["Days"], df["Cells Population"], label=key, color=colors[key], linewidth=2)

for csv_file in scaled_csvs:
    path = os.path.join(folder, csv_file)
    if os.path.exists(path):
        df = pd.read_csv(path, header=None, names=["Days", "Cells Population"])
        key = csv_file.split("-")[0]
        axes[1].plot(
            df["Days"], df["Cells Population"],
            label=f"{key} (scaled)", color=colors[key],
            linestyle="--", linewidth=1.8, alpha=0.9
        )

axes[0].set_ylabel("Cells Population (log scale)")
axes[0].set_title("Cells Population Over Time (separate CSVs: original & scaled)")
axes[0].set_yscale("log")
axes[0].legend(ncol=3)
axes[0].grid(True, which="both", linestyle="--", linewidth=0.5)

# -------------------------
# Plot 2: populations.csv (C,P,E,N,H,R vs days)
# -------------------------
pop_file = os.path.join(folder, "populations.csv")
df_pop = pd.read_csv(pop_file)

for col in ["C", "P", "E", "N", "H", "R"]:
    axes[1].plot(df_pop["days"], df_pop[col], label=col, color=colors[col], linewidth=2)

axes[1].set_xlabel("Days")
axes[1].set_ylabel("Cells Population (log scale)")
axes[1].set_title("Cells Population Over Time (from populations.csv)")
axes[1].set_yscale("log")
axes[1].legend(ncol=6)
axes[1].grid(True, which="both", linestyle="--", linewidth=0.5)

plt.tight_layout()
plt.show()
