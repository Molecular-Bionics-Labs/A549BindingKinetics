import numpy as np 
import matplotlib.pyplot as plt 
import pandas as pd 

import _pat_pmpc as pat 
import contour_plot as cplt 
from constant import C_rho 

import constant as const

A549_CellType = [61.2, 0.1, 334.1] # nPTM from humman protein atlas.  
kon_ves = 8.43*10**6
kon_mic= 6.92*10**5
dsh = const.dsh_A549  # The Glycocaly per cell.
ft = [32, 30, 28]
adjust = [0.9, 0.15, 0.15, 0.95]

nu_psome, nu_disk, nu_tube = 2.5, 1.37, 10.98 

# 2. PMPC particle data. 
R, polymerization = np.arange(10, 110, 5), np.arange(1, 106, 5) # Spherical

def fig2_avidity_theta_spherical():
    """Fig1(c) binding avidity."""
    # 1. Cell data. 
    cols = len(polymerization)
    rows = len(R)
    ka_arra = np.zeros([rows, cols])
    kd_arra = np.zeros([rows, cols])
    # 2. get ka and kd
    d0 = pat.pmpc_chain_length(25)
    R0 = 40
    all_rec = 61.2 + 0.1 + 334.1
    rec_list0 = [61.2*nu_psome/all_rec, 0.1*nu_psome/all_rec, 
                 334.1*nu_psome/all_rec]

    for col in range(cols):
        for row in range(rows):
            NPC = polymerization[col]
            rec_list = pat.pmpc_zeta_in_vb(rec_list0, R0, R[row])
            ka_i = pat.pmpc_ka(rec_list, dsh, R[row], NPC, phi=1)
            ka_arra[row, col] = ka_i 
            # Limitations of the ka 
            if ka_i == 0:
                kd_arra[row, col] = 1
            else:
                kd = 1/ka_i 
                if kd > 1: 
                    kd = 1 
                kd_arra[row, col] = kd
    # 4. Plot kd
    con_line = [0, 28, 8]

    if con_line[2]!=0: 
        arr = np.linspace(con_line[0], con_line[1], con_line[2])
        ticks = []
        tick_labels = [] 
        for i in range(len(arr)):
            ticks.append(arr[i])
            tick_labels.append(f'$10^{{{-arr[i]:.1f}}}$')

        cticks = {"cticks": ticks,
                "clables": tick_labels}
        
    label = [r"$N_{PC}$", "R(nm)", r"Spherical shap avidity $K_{D}(M)$", r"$K_{D}(M)$"]
    kd_log = -np.log10(kd_arra)
    cplt.contour_plot(polymerization, R, kd_log, MyColor="Gradient", label=label, 
                      cticks=cticks, con_line=con_line, pltfig=[12, 8], pltcbar=1,
                      ContourLable="10^-val")
    # plt.show()
    plt.savefig("./Figure/Fig2{0}_vesical_ka.png".format("A549"), dpi=800)


    # 5. Plot energy 
    # kd_log = np.log(1/kd_arra)
    # print(np.max(kd_log))
    # cticks = {"cticks": [0, 20, 40, 60, 80, 100, 120, 140],
    #           "clables":[0, -20, -40, -60, -80, -100, -120, -140]}
    # label = [r"$N_{PC}$", "R(nm)", "Binding energy (kT)", "Binding energy (kT)"]
    # cplt.contour_plot(polymerization, R, kd_log, MyColor="Gradient", 
    #                   pltfig=[10,7], cticks=cticks, label=label, con_line=[0, 140, 8] )
    # # plt.show()
    # plt.savefig("/Users/ttomis9651/Downloads/{0}Posome_ka.png".format("A549"), 
    #             dpi=800)

    # 6. Plot theta. 
    theta = ka_arra * C_rho / (1 + ka_arra * C_rho)
    label = [r"$N_{PC}$", "R(nm)", None, r"Binding fraction $\theta$"]
    # label = [r"$N_{PC}$", "R(nm)", r"Binding fraction $\theta$, $\rho=10^{-10}M$", r"Binding fraction $\theta$"]
    cticks = {"cticks": [0, 0.25, 0.5,0.75, 1],
              "clables": [0, 0.25, 0.5,0.75, 1]}
    cplt.contour_plot(polymerization, R, theta, MyColor=1, 
                      pltfig=[10,8],cticks=cticks, ft=ft, label=label, adjust=adjust, pltcbar=None)
    plt.savefig("./Figure/Fig2{0}_vesicle_theta.png".format("A549"), dpi=800)


R_disk = np.arange(8, 42, 2) # Spherical
def fig2_avidity_theta_disk():
    """Fig1(c) binding avidity."""
    # 1. Cell data. 
    cols = len(polymerization)
    rows = len(R_disk)
    ka_arra = np.zeros([rows, cols])
    kd_arra = np.zeros([rows, cols])
    # 2. get ka and kd
    vs = 2.5
    d0 = pat.pmpc_chain_length(25)
    R0 = 40
    all_rec = 61.2 + 0.1 + 334.1
    rec_list0 = [61.2*vs/all_rec, 0.1*vs/all_rec, 334.1*vs/all_rec]

    for col in range(cols):
        for row in range(rows):

            NPC = polymerization[col]
            rec_list = pat.pmpc_zeta_in_vb(rec_list0, R0, R_disk[row])
            
            ka_i = pat.pmpc_ka(rec_list, dsh, R_disk[row], NPC, phi=1, L=5, ParType="Disk")
            if ka_i > 10**56:
                ka_i = 10**56

            ka_arra[row, col] = ka_i 
            if ka_i == 0:
                kd_arra[row, col] = 1
            else:
                kd = 1/ka_i 
                if kd > 1: 
                    kd = 1 
                kd_arra[row, col] = kd


    # 4. Plot kd
       # 4. Plot kd
    con_line = [0, 35, 8]

    if con_line[2]!=0: 
        arr = np.linspace(con_line[0], con_line[1], con_line[2])
        ticks = []
        tick_labels = [] 
        for i in range(len(arr)):
            ticks.append(arr[i])
            tick_labels.append(f'$10^{{{-arr[i]:.1f}}}$')

        cticks = {"cticks": ticks,
                "clables": tick_labels}
    label = [r"$N_{PC}$", "R(nm)", r"Disk shap avidity $K_{D}(M)$", r"$K_{D}(M)$"]
    kd_log = -np.log10(kd_arra)
    cplt.contour_plot(polymerization, R_disk, kd_log, MyColor="Gradient", label=label, 
                      cticks=cticks, con_line=con_line, pltfig=[12, 8], pltcbar=1,
                      ContourLable="10^-val")
    
    plt.savefig("./Figure/Fig2{0}_Disk_ka.png".format("A549"), dpi=800)
    
    # 6. Plot theta. 
    theta = ka_arra * C_rho / (1 + ka_arra * C_rho)
    label = [r"$N_{PC}$", "R(nm)", None, r"Binding fraction $\theta$"]
    # label = [r"$N_{PC}$", "R(nm)", r"Binding fraction $\theta$, $\rho=10^{-10}M$", r"Binding fraction $\theta$"]
    cticks = {"cticks": [0, 0.25, 0.5,0.75, 1],
              "clables": [0, 0.25, 0.5,0.75, 1]}
    cplt.contour_plot(polymerization, R_disk, theta, MyColor=1, 
                      pltfig=[10,8],cticks=cticks, ft=ft, adjust=adjust, label=label, pltcbar=None)
    plt.savefig("./Figure/Fig2{0}_disk_theta.png".format("A549"), dpi=800)


R_cross = 20 # 
length = np.linspace(40, 1000, 96)
def figs3_avidity_theta_tube():
    """Plot the figure s2. The avidyty change with the tube legnth"""
    cols = len(polymerization)
    rows = len(length)
    kd_arra = np.zeros([rows, cols])
    ka_arra = np.zeros([rows, cols])

    for row in range(rows):
        for col in range(cols):
            all_rec = 61.2 + 0.1 + 334.1
            vs = nu_tube/300 * length[row]
            rec_list = [61.2*vs/all_rec, 0.1*vs/all_rec, 334.1*vs/all_rec]
            ka_i = pat.pmpc_ka(rec_list, dsh, 40, polymerization[col], 
                               phi=1, ParType="Tube", L=length[row])

            ka_arra[row, col] = ka_i 
            if ka_i == 0:
                kd_arra[row, col] = 1
            else:
                kd = 1/ka_i 
                if kd > 1: 
                    kd = 1 
                kd_arra[row, col] = kd

    con_line = [0, 175, 8]
    # 4. Plot kd
    if con_line[2]!=0: 
        arr = np.linspace(con_line[0], con_line[1], con_line[2])
        ticks = []
        tick_labels = [] 
        for i in range(len(arr)):
            ticks.append(arr[i])
            tick_labels.append(f'$10^{{{-arr[i]:.1f}}}$')

        cticks = {"cticks": ticks,
                "clables": tick_labels}
    label = [r"$N_{PC}$", "Length(nm)", r"Tube shap avidity $K_{D}(M)$", r"$K_{D}(M)$"]
    kd_log = -np.log10(kd_arra)
    cplt.contour_plot(polymerization, length, kd_log, MyColor="Gradient", label=label, 
                      cticks=cticks, con_line=con_line, pltfig=[12, 8], pltcbar=1,
                      ContourLable="10^-val")
    plt.savefig("./Figure/Fig2{0}_Tube_ka.png".format("A549"), dpi=800)

    # 6. Plot theta. 
    theta = ka_arra * C_rho / (1 + ka_arra * C_rho)
    label = [r"$N_{PC}$", "Length(nm)", None, r"Binding fraction $\theta$"]
    # label = [r"$N_{PC}$", "R(nm)", r"Binding fraction $\theta$, $\rho=10^{-10}M$", r"Binding fraction $\theta$"]
    cticks = {"cticks": [0, 0.25, 0.5,0.75, 1],
              "clables": [0, 0.25, 0.5,0.75, 1]}
    cplt.contour_plot(polymerization, length, theta, MyColor=1, 
                      pltfig=[10,8],cticks=cticks, ft=ft, adjust=adjust, label=label, pltcbar=None)
    plt.savefig("./Figure/Fig2{0}_Tube_theta.png".format("A549"), dpi=800)


fig2_avidity_theta_spherical()
fig2_avidity_theta_disk()
figs3_avidity_theta_tube()
