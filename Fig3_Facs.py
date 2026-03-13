import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
import scipy.optimize as sci
from scipy.integrate import solve_ivp
from scipy.integrate import odeint
from scipy.optimize import curve_fit

import _pat_pmpc as pat 

"""
This script has the fitting for the PMPC Psome binidng with a549 cells.
binding_odeint() is a function of bining. 
get_odeint_solution() is a function of resolution we can get form the binding. 
You need to the experiment time (t_exp), the data to fit (y_exp). 
Then, use curve_fit(fun_name, t_exp, y_exp, bounds=boundary) to fit the data. 
You can also choose other algrithm to fit your parameters. 
"""

# 1. Read the data. 
fp_read = "./Data/TB_Plus.xlsx"
NA = 6.022*10**23  
A549_CellType = np.array([61.2, 0.1, 334.1])

def binding_odeint(ft, t, nu, kon, koff, ke):
    """
    Here is the binding function which can be solved by the 
    scipy.integrate.odeint. 

    Parameters: 
    ----------
    ft: solution for the function. 

    t: time. 

    kon: on rate. 

    koff: off rate.

    ke: endocytosis rate. 
    """
    drhoPdt = -kon*ft[0]*ft[1]**nu + koff * ft[2]
    drhoRdt = -kon*ft[0]*ft[1]**nu + koff * ft[2]
    drhoPRdt = kon*ft[0]*ft[1]**nu - koff * ft[2] - ke*ft[2]

    drhoEdt = ke*ft[2]

    return [drhoPdt, drhoRdt, drhoPRdt, drhoEdt]

def binding_odeint_tube(ft, t, nu, kon, koff, ke):
    """
    Here is the binding function which can be solved by the 
    scipy.integrate.odeint. 

    Parameters: 
    ----------
    ft: solution for the function. 

    t: time. 

    kon: on rate. 

    koff: off rate.

    ke: endocytosis rate. 
    """
    kon = kon/rho0_tube**nu

    drhoPdt = -kon*ft[0]*ft[1]**nu + koff * ft[2]
    drhoRdt = -kon*ft[0]*ft[1]**nu + koff * ft[2]
    drhoPRdt = kon*ft[0]*ft[1]**nu - koff * ft[2] - ke*ft[2]

    drhoEdt = ke*ft[2]
    
    return [drhoPdt, drhoRdt, drhoPRdt, drhoEdt]

def get_odeint_solution(t, rhoR, nu, kon1, koff1, ke1, kon2, ke2, koff2):
    """
    Here is the fiting function. 

    Parameters:
    ----------
    fun_name: function name. 

    t_exp: experiment time point. 

    kon: on rate. 

    koff: off rate 

    ke: endocytosis rate. 
    """
    # 1. Initial values.
    # rhoP_ves, rhoR0, rhoPR0, rhoE0 =10**-8, 10**-5, 0, 0 # Giving at before. 
    rhoPR0, rhoE0 = 0, 0 
    rhoR_ves = rhoR/nu

    """PAT: get koff."""
    y_init = [rho0_ves, rhoR_ves, rhoPR0, rhoE0]
    args1 = (nu, kon1/rhoR_ves**nu, koff1, ke1)
    # 2. Solve function and get the result. 
    sol1 = odeint(binding_odeint, y_init, t_exp_ves, args=args1)
    y1 = sol1[:, 2]/np.max(sol1[:2])

    nu2 = nu/(40**2 / 23**2)
    rhoR_mic = rhoR / nu2 
    y_init = [rho0_mic, rhoR_mic, rhoPR0, rhoE0]
    args2 = (nu2, kon2/rhoR_mic**nu2, koff2, ke2)
    # 2. Solve function and get the result. 
    sol2 = odeint(binding_odeint, y_init, t_exp_mic, args=args2)
    y2 = sol2[:, 2]/np.max(sol2[:2])
    y = np.hstack((y1, y2))
    
    return y

def get_odeint_solution_tube(t, nu, kon, koff, ke):
    # 1. Initial values.
    # rhoP_ves, rhoR0, rhoPR0, rhoE0 =10**-8, 10**-5, 0, 0 # Giving at before. 
    rhoPR0, rhoE0 = 0, 0 
    y_init = [rho0_tube, rhoR_tube, rhoPR0, rhoE0]
    args1 = (nu, kon, koff, ke)
    # 2. Solve function and get the result. 
    sol1 = odeint(binding_odeint_tube, y_init, t_exp_tube_fit, args=args1)
    y = sol1[:, 2]/np.max(sol1[:2]) * 0.5

    return y


def pltfig(texp, yexp, yerror, tsim, ysim, colori, par_type, dashx=None, 
           dashy=None):
    plt.figure(figsize=(10,8))
    plt.errorbar(texp, yexp, yerr=yerror,     
                 fmt='o', # 'o' for scatter, no lines
                 markersize=8,  # Adjust scatter point size
                 capsize=5,     # Caps on error bars
                 capthick=2,    # Thickness of cap lines
                 elinewidth=2,  # Thickness of error bar lines
                 color=colori,   # Set color
                 alpha=0.8  # Transparency
                 )
    plt.grid(axis='y', linestyle='-', linewidth=0.7, alpha=0.7)
    plt.plot(tsim, ysim, color=colori)
    # Micelle 
    # plt.plot(t_sim[:4000], y_sim[:4000], color=color_i)
    # plt.plot(t_sim[4000:6000], y_sim[4000:6000], color=color_i, linestyle="--")
    fp_save = "/Users/ttomis9651/Downloads/{0}.svg".format(par_type)
    plt.xlim([0, 60])
    # plt.ylim([0, 1])
    if dashx is None:
        pass 
    else: 
        plt.plot(dashx, dashy, linestyle="--")
    plt.savefig(fp_save, dpi=800)
    # plt.show()

    # # Cover the kon into real kon. 
    # kon = popt[1]
    # # rhoR0 = 10**6/valence
    # print(valence)
    # kon_real = kon/(rhoR0/(valence*NA*10**-9))
    # print("{:.2e}".format(kon_real))


# 2. Initial parameters. 
rhoR0 = 0

ves_rho = 1.43*10**-8
mic_rho = 2.35*10**-7
tube_rho = 1.63*10**-9    

rho0_ves = ves_rho * NA * 10 **-9 / 15  # turn to NP per cell. 
rho0_tube = tube_rho * NA * 10 **-9 / 15  # turn to NP per cell. 
rho0_mic = mic_rho * NA * 10 **-9 / 15  # turn to NP per cell. 

boundary = np.array([[1.25*10**4, 2.5,  0.01*10**-1, 10**-9, 4*10**-2,  0.1, 10**-9, 1*10**-2],
                     [25*10**4, 4,  0.05*10**-1, 10**-2,     5*10**-2,   0.2, 10**-2, 2*10**-2]]) 
fitting_function = "curve_fit"

data_ves = pd.read_excel(fp_read, sheet_name="Psome")
data_mic = pd.read_excel(fp_read, sheet_name="Micelle")
data_tube = pd.read_excel(fp_read, sheet_name="Tube")

t_exp_ves = data_ves.values[:, 0]
y_ves = data_ves.values[:, 1]
error_ves = data_ves.values[:, 2]
y_max_ves = np.max(y_ves)
theta_ves = y_ves/y_max_ves 

t_exp_mic = data_mic.values[:, 0]
y_mic = data_mic.values[:, 1]
error_mic = data_mic.values[:, 2]
y_max_mic = np.max(y_mic)
theta_mic = y_mic/y_max_mic

t_exp_tube_fit = data_tube.values[:, 0]
y_tube_fit = data_tube.values[:, 1]
y_max_tube = np.max(y_tube_fit)
theta_tube_fit = y_tube_fit/y_max_tube

t_exp = np.hstack((t_exp_ves, t_exp_mic))

colors = [[0.8, 0.1, 0.1], [0, 0.392, 0], [0.804, 0.522, 0.247]]
plt.figure(figsize=(10,10))
if fitting_function=="curve_fit":
    # Fite the data. 
    theta = np.hstack((theta_ves, theta_mic))
    popt, *pcov = curve_fit(get_odeint_solution, t_exp, theta, bounds=boundary)
    print(popt)
    nu = popt[1]
    rhoR_ves = popt[0]/nu
    popt_ves = (popt[1], popt[2]/rhoR_ves**nu, popt[3],popt[4])

    y_init = [rho0_ves, rhoR_ves, 0, 0]
    t_sim = np.linspace(0, 65, 6501)
    y_sol = odeint(binding_odeint, y_init, t_sim, args=popt_ves)
    print("Occupied Receptor", np.max(y_sol[:, 2])/rhoR_ves)

    y_sol = y_sol[:, 2]/ np.max(y_sol[:,2]) * np.max(y_ves)

    plt.errorbar(t_exp_ves, y_ves, error_ves,   
                fmt='o', # 'o' for scatter, no lines
                markersize=14,  # Adjust scatter point size
                capsize=8,     # Caps on error bars
                capthick=3,    # Thickness of cap lines
                elinewidth=3,  # Thickness of error bar lines
                color=colors[0],   # Set color
                alpha=0.8  # Transparency
                )
    plt.grid(axis='y', linestyle='-', linewidth=0.7, alpha=0.7)
    plt.plot(t_sim, y_sol, color=colors[0], linewidth=3)

    t_exp_mic = data_mic.values[:, 0]
    y_mic = data_mic.values[:, 1]
    error_mic = data_mic.values[:, 2]
    y_max_mic = np.max(y_mic)
    theta_mic = y_mic/y_max_mic
    """Micelle"""
    nu_mic = nu/(40**2/28**2)
    rhoR_mic = popt[0]/nu_mic
    popt_mic = (nu_mic, popt[5]/rho0_mic**nu_mic, popt[6], popt[7])
    # Simulation. 
    y_init = [rho0_mic, rhoR_mic, 0, 0]
    y_sol = odeint(binding_odeint, y_init, t_sim, args=popt_mic)
    print("Occupied Receptor", np.max(y_sol[:,2])/rhoR_mic)
    y_sol = y_sol[:, 2]/ np.max(y_sol[:,2]) * y_max_mic
    """Micelle plot."""
    plt.errorbar(t_exp_mic, y_mic, error_mic,   
                fmt='o', # 'o' for scatter, no lines
                markersize=14,  # Adjust scatter point size
                capsize=8,     # Caps on error bars
                capthick=3,    # Thickness of cap lines
                elinewidth=3,  # Thickness of error bar lines
                color=colors[1],   # Set color
                alpha=0.8  # Transparency
                )
    plt.grid(axis='y', linestyle='-', linewidth=0.7, alpha=0.7)
    plt.plot(t_sim, y_sol, color=colors[1], linewidth=3)

    """Use PAT to modify the koff."""
    dsh = 26
    ves_list = nu*A549_CellType/ np.sum(A549_CellType)
    ka_Psome = pat.pmpc_ka(ves_list, dsh, 40, 25, phi=1)

    kon_Psome = popt[2]/(rhoR_ves*15/(10**-9*NA))**nu
    koff_Psome = kon_Psome / ka_Psome
    print("nu_Psomes=", "{:.3e}".format(popt[1]), 
          "; kon_Psomes=", "{:.3e}".format(ka_Psome), 
          "; kon_Psome=", "{:.3e}".format(kon_Psome), 
          "; koff_Psome=", "{:.3e}".format(koff_Psome))

    mic_list = nu_mic*A549_CellType/ np.sum(A549_CellType)
    ka_mic = pat.pmpc_ka(mic_list, dsh, 28, 25, phi=1)

    kon_mic = popt[5]/(rhoR_mic*15/(10**-9*NA)) ** nu_mic
    koff_mic = kon_mic / ka_mic
    print("nu_mic=","{:.3e}".format(nu_mic), 
          "; ka_mic=","{:.3e}".format(ka_mic), 
          "; kon_mic=","{:.3e}".format(kon_mic), 
          "; koff_mic=","{:.3e}".format(koff_mic))
    
    """The density of the NPs."""
    c_R = nu / (40**2*np.pi) 
    nu_tube = 20* 200*c_R 
    tube_list = nu_tube*A549_CellType/ np.sum(A549_CellType)
    ka_mic = pat.pmpc_ka(tube_list, dsh, 28, 25, phi=1)
    rhoR_tube = popt[0]/nu_tube 
 
    y_init = [rho0_tube, rhoR_tube, 0, 0]
    boundary = np.array([[1.9, 0.01, 10**-9, 10**-6], 
                        [2,  0.02, 10**-2, 10**-5]])

    popt, *pcov = curve_fit(get_odeint_solution_tube, t_exp_tube_fit, theta_tube_fit, bounds=boundary)

    args_tube = (popt[0], popt[1], popt[2], popt[3])
    y_sol = odeint(binding_odeint_tube, y_init, t_sim, args=args_tube)
    print("Occupied Receptor", np.max(y_sol[:,2])/rhoR_tube)
    y_sol = y_sol[:,2]/np.max(y_sol[:, 2]) *  y_max_tube
    t_exp_tube = data_tube.values[:, 0]
    y_tube = data_tube.values[:, 1]
    error_tube = data_tube.values[:, 2]
    y_max_tube = np.max(y_tube)
    theta_tube = y_tube/y_max_tube

    plt.errorbar(t_exp_tube, y_tube, error_tube,   
            fmt='o', # 'o' for scatter, no lines
            markersize=14,  # Adjust scatter point size
            capsize=8,     # Caps on error bars
            capthick=3,    # Thickness of cap lines
            elinewidth=3,  # Thickness of error bar lines
            color=colors[2],   # Set color
            alpha=0.8  # Transparency
            )
    plt.grid(axis='y', linestyle='-', linewidth=0.7, alpha=0.7)
    plt.plot(t_sim, y_sol, color=colors[2], linewidth=3)


# Micelle 
# plt.plot(t_sim[:4000], y_sim[:4000], color=color_i)
# plt.plot(t_sim[4000:6000], y_sim[4000:6000], color=color_i, linestyle="--")

plt.xlim([0, 65])
plt.ylim([0, 80000])

plt.yticks([0, 2*10**4, 4*10**4, 6*10**4, 8*10**4])
# plt.legend()

fp_save = "/Users/ttomis9651/Downloads/{0}.svg".format("Facs")
plt.savefig(fp_save, dpi=800)
# plt.show()
    