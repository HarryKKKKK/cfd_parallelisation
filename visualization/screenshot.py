import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("t=0.csv")

plt.tricontourf(df["x"], df["y"], df["rho"], levels=100)
plt.colorbar(label="Density")
plt.axis("equal")
plt.xlabel("x")
plt.ylabel("y")
plt.title("rho")
plt.savefig("fig.png")
