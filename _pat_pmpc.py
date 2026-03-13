import numpy as np
NA = 6.02214076 * 10 ** 23

"""
Functions: 
d = pmpc_chain_length(NPC)
reclist = pmpc_zeta_in_vb(rec_list, NPC)
"""

def pmpc_chain_length(NPC):
    """
    Calculate the pmpc chain lenght.  return: chain length. 

    Parameter:
    ----------
    NPC: polymerization of the PMPC chain. 
    
    """
    aPC = 0.285
    d = NPC * aPC/3  # Devied by 3 condidering the rotation of PC units. 
    return d 


def pmpc_zeta_in_vb(rec_list0, R0, R):
    """
    Renew receptor number in the binding volume. Return: receptor list. 

    Nanoparticle are assumned to bind at a distacne d=Rg from the surface.
    The interact area can be calculated as: S = [(R+Rg)**2 -R**2]*pi, 
    which is approximate to be 2*pi*R*Rg. 
    You can find the method in SI: 10.1038/s41467-020-18603-5

    Parameters:
    -----------
    rec_list0: receptor list for particle with parameter [R0, d0].

    R0: Radius

    d0: Chain lenfth

    R: Radius for the new particle. 

    NPC: the polymerization. 
    """
    rec_list = [] 
    for i in range(len(rec_list0)):
        rec_i = rec_list0[i]*R/R0
        rec_list.append(rec_i)

    return rec_list


def pmpc_ka(rec_list, dSH, R, polymerization, L=200, phi=None, 
            ParType="Vesicle", Print=False, retrun_type="avidity"):
    """ 
    This fucntion is simplified function for PMPC pat. 

    Parameter:
    --------------
    Receptor: list of float. The receptor expression from the protein atlas and 
    fluorescence.

    dSH: float. Distance between two glycosaminoglycans(GAGs).

    R: float. Radius of the pmpc particle.

    polymerization: int. Polymerization of a single pmpc chain.

    L: float. An argument in the length of the  tube pmpc polymersome.

    Phi: float. Transfer the receptor expression into receptor number within vB.

    ParType: str. Morphology of the NP: micelle, vesicle, and tube.

    Print: Print the information of the particle data.
    
    Example:
    --------
    configure: = {}

    Attension: 
    ----------

    Don't use (R-d) this will casuse an error when d>R (NPC>100, R<15).

    """
    NPC = round(polymerization * 0.34)  # Roation of the PMPC chain. 
    # Receptor expression
    SRB1 = rec_list[0]
    CD36 = rec_list[1]
    CD81 = rec_list[2]
    sigmaGAG = np.pi * dSH ** 2

    # Get the receptor number in the binding volume.
    if phi is None:
        phi = 0.116
        n_zeta1, n_zeta2, n_zeta3 = SRB1 * phi, CD36 * phi, CD81 * phi
    elif type(phi) == float:
        n_zeta1, n_zeta2, n_zeta3 = SRB1 * phi, CD36 * phi, CD81 * phi
    elif type(phi) == list:
        n_zeta1, n_zeta2, n_zeta3 = SRB1 * phi[0], CD36 * phi[1], CD81 * phi[2]
    else:
        n_zeta1, n_zeta2, n_zeta3 = SRB1 * phi, CD36 * phi, CD81 * phi


    v_zeta1, v_zeta2, v_zeta3 = 49.61, 49.61, 11.53  # the volume of the receptor
    v_zeta = [v_zeta1, v_zeta2, v_zeta3]  # the volume of the receptor

    # The steric potential of juxtaposition
    u_zeta1 = 0.38169024 / 0.5992
    u_zeta2 = 0.36030411 / 0.5992
    u_zeta3 = 0.13455454 / 0.5992
    u_zeta = [u_zeta1, u_zeta2, u_zeta3]  # The steric potential of juxtaposition
    b_zeta1, b_zeta2, b_zeta3 = 0.07846136, 0.05445527, -0.02074469
    b_zeta = [b_zeta1, b_zeta2, b_zeta3]  # the juxtaposing bind sites
    h_zeta1, h_zeta2, h_zeta3 = 6.86, 6.86, 2.25
    h_zeta = [h_zeta1, h_zeta2, h_zeta3]  # the height of the receptor

    # Binding energy for a single PMPC unit.
    EB_SRB1 = -4.46050141 / 0.5992 + 4   # SRB1 = -4.46050141 / 0.5992 + 4
    EB_CD36 = -4.95330502 / 0.5992 + 7  # CD36 = -4.95330502 / 0.5992 + 7
    EB_CD81 = -3.53280759 / 0.5992 + 2 # CD81 = -3.53280759 / 0.5992 + 0
    EB = [EB_SRB1, EB_CD36, EB_CD81]

    # Calculate the number of the PMPC chain.
    Detailed_Calculation = 0  # Calculate the PMPC chain in a precision way. 
    aPC = 0.285
    aDPA = 0.285
    poly_DPA = 70
    DPA = 213.32
    MPDPA = DPA * poly_DPA
    rho_PDPA = 0.9
    # PC = 259.27
    d = aPC * polymerization ** 0.9
    t = aDPA * poly_DPA ** (2 / 3)
    v_PDPA = MPDPA / (rho_PDPA * NA) * 10 ** 21
    N_PMPC = 0
    sigma_PC = 6.17

    if Detailed_Calculation == 1:
        if ParType == "Vesicle":
            V_ves = 4/3*np.pi*((R - d)**3-(R - d - t)**3)
            Sex = 4 * np.pi * (R-d)**2
            Sin = 4 * np.pi * (R - d - t)**2
            xout = (1-np.sqrt(Sin/Sex))/(1-Sin/Sex)
            N_PMPC = V_ves / v_PDPA * xout

        elif ParType == "Micelle":
            p_mi = 1/3
            V_mi = 4/3 * np.pi * (R-d) ** 3
            N_PMPC = V_mi / v_PDPA * p_mi

        elif ParType == "Tube":
            p_tu = 0.774
            V_tu = np.pi * ((R-d)**2 - (R-d-t)**2) * L
            N_PMPC = V_tu / v_PDPA * p_tu

        elif ParType == "Disk":
            Sur_disk = 2*(np.pi*R**2)+2*np.pi*R*L
            N_PMPC = V_tu / v_PDPA * p_tu
        else:
            print("Please put the correct particle type.")

    else:
        if ParType == "Vesicle":
            Sur_ves = 4 * np.pi * R ** 2
            N_PMPC = Sur_ves / sigma_PC

        elif ParType == "Micelle":
            Sur_mi = 4 * np.pi * R ** 2
            N_PMPC = Sur_mi / sigma_PC

        elif ParType == "Tube":
            Sur_tu = np.pi * R ** 2 * L
            N_PMPC = Sur_tu / sigma_PC
        elif ParType == "Disk":
            Sur_disk = 2*(np.pi*R**2)+2*np.pi*R*L
            N_PMPC = Sur_disk / sigma_PC
        else:
            print("Please put the correct particle type.")

    N_PMPC = np.floor(N_PMPC)

    # 3. Binding energy for single chain
    q_zeta = np.zeros([9, NPC])  # Check the binding energy, q and single bound. 
    E_zeta = []
    for i in range(0, 3, 1):
        # 3 different receptors.
        e_zeta = 0
        for j in range(1, NPC+1, 1):
            n1 = 1 + b_zeta[i] * (j - 1)
            # Here is no degeneracy.
            q1 = -np.log((1 + np.exp(- EB[i] - u_zeta[i] * (j - 1))) ** n1 - 1)
            GC = aPC * j / R

            if GC <= 3 ** 0.5 - 1:
                gamma = (GC + 1) ** 2
            else:
                gamma = 3

            numerator = (NPC - j + 1) * aPC * v_zeta[i] * (1 - ((NPC - j + 1) / NPC) ** 2) ** (9 / 4)
            denominator = h_zeta[i] * (sigma_PC * (1 + aPC * (NPC - j) / R) ** (gamma - 1)) ** (3 / 2)
            q2 = numerator / denominator

            e_zeta = e_zeta + q1 + q2
            q_zeta[3*i, j-1] = q1
            q_zeta[3*i+1, j-1] = q2
            q_zeta[3*i+2, j-1] = e_zeta
        # Revise the e_zeta
        if e_zeta >= EB[i]:
            e_zeta = EB[i]
        E_zeta.append(e_zeta)

    # 4. Total binding energy
    x1, x2, x3 = N_PMPC * np.exp(-E_zeta[0]), N_PMPC * np.exp(-E_zeta[1]), N_PMPC * np.exp(-E_zeta[2])

    ka1 = (1 + x1) ** n_zeta1 - 1
    ka2 = (1 + x2) ** n_zeta2 - 1
    ka3 = (1 + x3) ** n_zeta3 - 1

    E_1, E_2, E_3 = np.log1p(ka1), np.log(ka2), np.log(ka3)
    E_total = E_1 + E_2 + E_3

    # 5. the steric potential of the glycosaminoglycan.
    v_Par = 0
    if ParType == "Vesicle":
        v_Par = 4 / 3 * np.pi * R ** 3
    elif ParType == "Tube":
        v_Par = (np.pi*R**2) * L
    elif ParType == "Micelle":
        v_Par = 4 / 3 * np.pi * R ** 3
    elif ParType == "Disk":
        v_Par = (np.pi*R**2) * L # Here L means the thickness of the disk.
    else:
        v_Par = 4 / 3 * np.pi * R ** 3

    u_gag = v_Par * (sigmaGAG ** (-3 / 2))
    # if ParType =="Tube":
    #     print("Here is the Ugag", u_gag)
    # 6. Total binding energy.
    E_F = E_total - u_gag

    if np.issubdtype(type(E_F), np.float64) and E_F < -199:
        E_F = -100
    elif np.issubdtype(type(E_F), np.ndarray):
        E_F[E_F <= -199] = -100

    if np.issubdtype(type(E_F), np.float64) and E_F > 400:
        E_F = 400
    elif np.issubdtype(type(E_F), np.ndarray):
        E_F[E_F >= 200] = 400

    qpsome = np.exp(E_F)

    # The binding proportion.
    equal_unit = 10 ** -8
    R1 = R * equal_unit
    d = (aPC * equal_unit * NPC)

    if ParType == "Vesicle":
        v_B = np.pi * (3 * (R1+d) ** 3 - R1 ** 3) / 3 * NA
    elif ParType == "Micelle":
        v_B = np.pi * (3 * (R1+d) ** 3 - R1 ** 3) / 3 * NA
    elif ParType == "Tube":
        v_B = np.pi*(3*((R1+d)**2*L*equal_unit)-R1**2*L*equal_unit)/3*NA
    elif ParType == "Disk":
        v_B = np.pi * (3 * (R1+d) **2*L*equal_unit - R1**2*L*equal_unit)/3*NA
    else:
        print("Please input the correct particle type.")
    # numerator_theta_right = NA * v_B * qpsome
    avidity = v_B * qpsome

    if retrun_type == "Energy":
        return E_F
    elif retrun_type == "Repel_Binding":
        return u_gag, -E_total
    elif retrun_type == "avidity":
        if Print is True:
            """Print settings"""
            print("PC unit EB with receptor(kT):", "\n", "(SRB1,CD36,CD81)         ", EB)
            print("Single chain E_zeta with receptor(kT): ", "\n", "(SRB1,CD36,CD81)         ", E_zeta)
            if ParType == "Tube":
                print("Average length of tube is(nm):", L)
            print(ParType, ": R(nm)=", R, ", N_PMPC=", N_PMPC, ", sigma_PC(nm**2):", sigma_PC)
            print("U_GAG(kT):", u_gag, "Qsome(ka):", qpsome, "     Binding Energy(kT):", -np.log(qpsome))
            print("qsome", qpsome, "V_B", v_B)
            print("Single chain molar weight: ")
            print("ka (number per liter):", avidity)
            print("The hydrophilic tether length d: ", d, "The hydrophobic length t:", t, "\n")

        return avidity
    