import pickle
from scipy.stats import nbinom
import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression

def predict_nz(PCR_cyc, MNZ, N_reads):
    ans = np.zeros((2, 6))
    with open("sPCR_sim_1.RData", "rb") as f:
        pcr_sim = pickle.load(f)
    
    min_diff = np.min(np.abs(pcr_sim[:, 6] - MNZ))
    ans[0, :] = pcr_sim[np.abs(pcr_sim[:, 6] - MNZ) == min_diff, :][0]
    ans[1, :] = pcr_sim[np.abs(pcr_sim[:, 6] - MNZ) == min_diff, :][0]
    return ans


def predict_nz1(PCR_cyc, MNZ, N_reads):
    ans = np.zeros((2, 6))
    with open("sPCR_sim_1.RData", "rb") as f:
        pcr_sim = pickle.load(f)
    
    inds = (np.abs(pcr_sim[:, 3] - N_reads) <= 5000) & (pcr_sim[:, 4] == PCR_cyc) & (np.abs(pcr_sim[:, 6] - MNZ) < 0.025)
    if np.sum(inds) == 0:
        inds = (np.abs(pcr_sim[:, 3] - N_reads) <= 5000) & (pcr_sim[:, 4] == PCR_cyc) & (np.abs(pcr_sim[:, 6] - MNZ) < 0.05)
        print("No appropriate simulation")
    
    if np.sum(inds) == 0:
        raise ValueError("No appropriate simulation")
    
    ans[0, :] = pcr_sim[pcr_sim[:, 5] == np.max(pcr_sim[inds, 5]), :][0]
    ans[1, :] = pcr_sim[pcr_sim[:, 5] == np.min(pcr_sim[inds, 5]), :][0]
    return ans

def sc_saturation_sub_sampling(r1, scale_by=np.arange(0.1, 1.1, 0.1)):
    n_captured = []
    f = [0] * (max(r1) + 1)
    for i in r1:
        f[i] += 1
    
    for scale in scale_by:
        captured = 0
        for x in range(len(f)):
            captured += f[x] * (1 - (1 - scale)**x)
        n_captured.append(captured)
    
    return n_captured

def sc_saturation(r1, pcr_cyc, scale_by=np.arange(1, 3.1, 0.1), verbose=False):
    if scale_by[0] != 1:
        raise ValueError("scale_by should start with 1")
    
    r2 = r1[:min(len(r1), 15000)]
    x = np.arange(1, max(r2) + 1)
    nUMI = len(r2)
    y = np.zeros(max(r2) + 1)
    
    for i in r2:
        y[i] += 1
    
    fit = LinearRegression().fit(x[:2].reshape(-1, 1), y[:2])
    y0_lin = fit.predict(np.array([[0]]))[0]
    
    if verbose:
        print(f"Running predict with {pcr_cyc}, {np.mean(r2)}, {np.sum(r2)}")
    
    p = predict_nz(pcr_cyc, np.mean(r2), np.sum(r2))
    print(p)
    y0_upp = nUMI * (p[0, 5]) / (1 - p[0, 5])
    y0_low = nUMI * (p[1, 5]) / (1 - p[1, 5])
    
    if verbose:
        plt.plot(np.append([0], x), np.append([y0_lin], y), ylim=[0, max(y0_lin, y0_upp, y0_low, max(y))])
        plt.plot(x[:2], fit.predict(x[:2].reshape(-1, 1)), color='red')
        plt.scatter(0, y0_upp, color='green')
        plt.scatter(0, y0_low, color='blue')
        plt.show()
    
    nb = nbinom.fit(np.concatenate([np.zeros(int(y0_upp)), r2]), floc=0)
    if verbose:
        plt.figure()
        plt.hist(np.concatenate([np.zeros(int(y0_upp)), r2]), bins=50, density=True, alpha=0.6, color='g')
        plt.show()
    
    s, m = nb[0], nb[2]
    perc_zero = [nbinom.pmf(0, s, s / (s + scale * m)) for scale in scale_by]
    N = round(len(r1) / (1 - perc_zero[0]))
    n_captured = [N * (1 - perc_zero[i]) for i in range(len(scale_by))]
    
    return n_captured
