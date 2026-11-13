import matplotlib.pyplot as plt
import pandas as pd


def plot_dose_comparison(trace_paths, dosages):
    fig, ax = plt.subplots(figsize=(10, 6))

    for i, (trace_path, dosage) in enumerate(zip(trace_paths, dosages)):
        trace_data = pd.read_csv(trace_path)
        ax.plot(trace_data['time'][:4500], trace_data['v'][:4500], label=f'Dose {dosage} nM')

    ax.set_xlabel('Time (ms)')
    ax.set_ylabel('Voltage (mV)')
    ax.set_title('Voltage Across Doses')
    ax.legend()
    plt.show()


if __name__ == "__main__":
    trace_paths = [
        "results\\ap\\fixed\\control\\dose_0\\trace.csv",
        "results\\ap\\fixed\\quinidine\\dose_1\\trace.csv",
        "results\\ap\\fixed\\quinidine\\dose_2\\trace.csv",
        "results\\ap\\fixed\\quinidine\\dose_3\\trace.csv",
        "results\\ap\\fixed\\quinidine\\dose_4\\trace.csv"
    ]

    dosages = [0, 3237, 6474, 9711, 12948]
    plot_dose_comparison(trace_paths, dosages)