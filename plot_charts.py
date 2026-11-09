import matplotlib
# matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd



def plot_currents(trace_path):
    trace = pd.read_csv(trace_path)

    time = trace.time
    voltage = trace.v

    INa = trace.INa
    INaL = trace.INaL
    Ito = trace.Ito
    ICaL = trace.ICaL
    IKr = trace.IKr
    IKs = trace.IKs
    IK1 = trace.IK1
    dv = trace.dv

    fig, axes = plt.subplots(2, 4, figsize=(12, 10))
    axes[0][0].plot(time, voltage)
    axes[0][0].set_xlabel("Time (ms)")
    axes[0][0].set_ylabel("Voltage (mV)")
    axes[0][0].set_title("Voltage Trace")

    axes[0][1].plot(time, INa, label="INa")
    axes[0][1].set_xlabel("Time (ms)")
    axes[0][1].set_title("INa")

    axes[0][2].plot(time, INaL, label="INaL")
    axes[0][2].set_xlabel("Time (ms)")
    axes[0][2].set_ylim((-2, 2))
    axes[0][2].set_title("INaL")

    axes[0][3].plot(time, Ito, label="Ito")
    axes[0][3].set_xlabel("Time (ms)")
    axes[0][3].set_ylim((-2, 2))
    axes[0][3].set_title("Ito")

    axes[1][0].plot(time, ICaL, label="ICaL")
    axes[1][0].set_xlabel("Time (ms)")
    axes[1][0].set_ylim((-2, 2))
    axes[1][0].set_title("ICaL")

    axes[1][1].plot(time, IKr, label="IKr")
    axes[1][1].set_xlabel("Time (ms)")
    axes[1][1].set_ylim((-2, 2))
    axes[1][1].set_title("IKr")

    axes[1][2].plot(time, IKs, label="IKs")
    axes[1][2].set_xlabel("Time (ms)")
    axes[1][2].set_ylim((-2, 2))
    axes[1][2].set_title("IKs")

    axes[1][3].plot(time, IK1, label="IK1")
    axes[1][3].set_xlabel("Time (ms)")
    axes[1][3].set_ylim((-2, 2))
    axes[1][3].set_title("IK1")

    plt.tight_layout()
    plt.show()


    fig2, axes2 = plt.subplots(3, 1, figsize=(8, 6))
    axes2[0].plot(time[:800], voltage[:800])
    axes2[0].set_ylabel("Voltage (mV)")
    axes2[0].set_title("Voltage Trace")

    axes2[1].plot(time[:800], INaL[:800], label="INaL")
    axes2[1].plot(time[:800], Ito[:800], label="Ito")
    axes2[1].plot(time[:800], ICaL[:800], label="ICaL")
    axes2[1].plot(time[:800], IKr[:800], label="IKr")
    axes2[1].plot(time[:800], IKs[:800], label="IKs")
    axes2[1].plot(time[:800], IK1[:800], label="IK1")
    axes2[1].set(xlabel="Time (ms)", ylabel="Ionic Current")
    axes2[1].legend(loc="right")
    axes2[1].set_title("Ionic Currents")

    axes2[2].plot(time[:800], INa[:800], label="INa")
    axes2[2].set_title("INa")
    axes2[2].set_xlabel("Time (ms)")
    axes2[2].set_ylabel("Current (pA/pF)")

    plt.tight_layout()
    plt.show()

    fig3, axes3 = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
    axes3[0].plot(time[:500], voltage[:500])
    axes3[0].set_ylabel("Voltage (mV)")
    axes3[1].plot(time[:500], dv[:500])
    axes3[1].set(xlabel="Time (ms)", ylabel="dv/dt (mV/ms)")

    plt.tight_layout()
    plt.show()
    
    

# keys within trace dataframe
    # ['time', 'v', 'nai', 'nass', 'ki', 'kss', 
    #  'cai', 'cass', 'cansr', 'cajsr', 'm', 'hf', 
    #  'hs', 'j', 'hsp', 'jp', 'mL', 'hL', 'hLp', 
    #  'a', 'iF', 'iS', 'ap', 'iFp', 'iSp', 'd', 
    #  'ff', 'fs', 'fcaf', 'fcas', 'jca', 'nca', 
    #  'ffp', 'fcafp', 'xs1', 'xs2', 'xk1', 'Jrelnp', 
    #  'Jrelp', 'CaMKt', 'IC1', 'IC2', 'C1', 'C2', 
    #  'O', 'IO', 'IObound', 'Obound', 'Cbound', 'D', 
    #  'qNet', 'INa', 'INaL', 'Ito', 'ICaL', 'IKr', 
    #  'IKs', 'IK1', 'dv']


if __name__ == "__main__":
    trace_path = "results\\ap\\fixed\\control\\dose_0\\trace.csv"
    plot_currents(trace_path)