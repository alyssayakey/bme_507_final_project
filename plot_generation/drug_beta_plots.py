import matplotlib.pyplot as plt
import pandas as pd

def plot_drug_beta_tests(trace_paths, dosages, iso_concs):
    fig, ax = plt.subplots(2, 1,figsize=(10, 6))
    
    for i, (trace_path, dosage) in enumerate(zip(trace_paths, dosages)):
        trace_data = pd.read_csv(trace_path)

        time = trace_data.time
        voltage = trace_data.v
    
        INa = trace_data.INa
        INaL = trace_data.INaL
        Ito = trace_data.Ito
        ICaL = trace_data.ICaL
        IKr = trace_data.IKr
        IKs = trace_data.IKs
        IK1 = trace_data.IK1
        dv = trace_data.dv

        ax[0].plot(time[10:1000], voltage[10:1000], label=f'Quinidine: {dosage} nM, Iso: {iso_concs[i]} uM')
        ax[1].plot(time[10:1000], dv[10:1000])

    ax[0].set_xlabel('Time (ms)')
    ax[0].set_ylabel('Voltage (mV)')
    ax[1].set_xlabel('Time (ms)')
    ax[1].set_ylabel('dv/dt')
    fig.suptitle('dv/dt Across Experimental Conditions (without first 10ms)')
    fig.legend()
    plt.show()



if __name__ == "__main__":
    trace_paths = [
        "results\\cipa_beta\\iso_0_uM\\fixed\\quinidine\\dose_0\\trace.csv",
        "results\\cipa_beta\\iso_0_uM\\fixed\\quinidine\\dose_1\\trace.csv",
        "results\\cipa_beta\\iso_1_uM\\fixed\\quinidine\\dose_0\\trace.csv",
        "results\\cipa_beta\\iso_1_uM\\fixed\\quinidine\\dose_1\\trace.csv"
    ]

    drug_dosages = [0, 3237, 0, 3237]
    iso_concs = [0, 0, 1, 1]
    plot_drug_beta_tests(trace_paths, drug_dosages, iso_concs)