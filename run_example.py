from gong_beta import run_cell
import matplotlib.pyplot as plt

results = {}
for iso in [0.0, 1.0]:
    print(f"Running ISO = {iso} uM ...", flush=True)
    results[iso] = run_cell(iso=iso, beats=5, pcl=1000)
    print("Last-beat APD90:", round(results[iso]["history"][-1, 1], 2), "ms")

fig, axes = plt.subplots(1, 2, figsize=(11, 4))
for iso, result in results.items():
    axes[0].plot(result["time"], result["voltage"], label=f"ISO={iso} uM")
    axes[1].plot(result["time"], result["calcium_uM"], label=f"ISO={iso} uM")
axes[0].set(xlabel="Time (ms)", ylabel="Voltage (mV)", title="Last AP of a five-beat check")
axes[1].set(xlabel="Time (ms)", ylabel="Calcium (uM)", title="Calcium transient")
for ax in axes:
    ax.legend()
    ax.grid(alpha=0.3)
plt.tight_layout()
plt.show()
