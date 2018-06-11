import numpy as np 
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import ctm_single as ctm


def simulate_openloop(rho_0, P, T, Res):
    '''
    simulate the open-loop ctm for a single initial condition
    rho_0: n x 1 [doulbe] initial condition of rho
    P: (d, C_d, C, epsilon, vf, w, wbar, rho_j, rho_jbar)
    return
    T: Termination time
    Res: number of simulation steps in each minute
    rho: n x T [double] trajectory of rho
    q:   (n+1) x T [double] flow rate
    '''

    # time points
    t = np.linspace(0, T, T*60*Res + 1)
    rho = np.zeros((1, len(t)))
    q = np.zeros((2, len(t)))
    rho[:,0] = rho_0
    q[:,0] = ctm.q_openloop(rho_0, P)

    for i in range(1, len(t)):
        # time span for the simulation step
        tspan = [t[i-1], t[i]]
        temp = odeint(ctm.ctm_single_openloop, rho[:,i-1], tspan, args=(P,))
        rho[:, i] = temp[-1,:]
        q[:,i] = ctm.q_openloop(rho[:, i], P)

    return rho, q

def simulate_closedloop(rho_0, P, delta, T, Res):
    '''
    simulate the open-loop ctm for a single initial condition
    rho_0: n x 1 [doulbe] initial condition of rho
    P: (d, C_d, C, epsilon, vf, w, wbar, rho_j, rho_jbar)
    delta: (delta_1, delta_2) control constants
        T: Termination time
    Res: number of simulation steps in each minute
    return
    rho: n x T [double] trajectory of rho
    q:   (n+1) x T [double] flow rate
    v:   n x T [double] VSL
    '''
    d, C_d, C, epsilon, vf, w, wbar, rho_j, rho_jbar = P
    delta_1, delta_2 = delta
    # time points
    t = np.linspace(0, T, T*60*Res+1)
    # density profile
    rho = np.zeros((1, len(t)))
    # flow rate profile
    q = np.zeros((2, len(t)))
    # VSL profile
    v = np.zeros((1, len(t)))

    # initial condition
    rho[:, 0] = rho_0
    # Control Mode
    if rho_0 > C_d / vf:
        cMode = 1
    else:
        cMode = 2
    q[:, 0], v[:, 0] = ctm.q_closedloop(rho[:,0], P, delta, cMode)

    for i in range(1, len(t)):

        # time span for the simulation step
        tspan = [t[i-1], t[i]]
        temp = odeint(ctm.ctm_single_closedloop, rho[:, i-1], tspan, args=(P, delta, cMode))
        rho[:,i] = temp[-1, :]
        q[:, i], v[:, i] = ctm.q_closedloop(rho[:, i], P, delta, cMode)

        # Check control mode
        if cMode == 1 and rho[:, i][-1] < C_d/vf - delta_2:
            cMode = 2

    return rho, q, v


# def phase_portrait(n, P, delta, T, Res, ctrl):
#     '''
#     Find the phase portrait of the system
#     n: Number of road sections
#     P: (d, C_d, C, epsilon, vf, w, wbar, rho_j, rho_jbar)
#     delta: (delta_1, delta_2)
#     T: Termination time
#     Res: number of simulation steps in each minute
#     ctrl: control flag. If ctrl == 0: open-loop. If ctrl == 1: closed-loop 
#     '''
#     d, C_d, C, epsilon, vf, w, wbar, rho_j, rho_jbar = P
#     N = 10
#     # Set of initial conditions
#     ICs = np.vstack((np.zeros(N), np.linspace(0, rho_j, N)))
#     ICs = np.hstack((ICs, np.vstack((rho_j*np.ones(N), np.linspace(0, rho_j, N)))))
#     ICs = np.hstack((ICs, np.vstack((np.linspace(0, rho_j, N), np.zeros(N)))))
#     ICs = np.hstack((ICs, np.vstack((np.linspace(0, rho_j, N), rho_j*np.ones(N)))))

#     Rho = np.zeros((n, int(T*60*Res + 1), N * 4))
#     phaseFig = plt.figure(figsize=(8,6.5))
#     for i in range(4*N):
#     #for i in range(2):
#         if ctrl==0: 
#             Rho[:, :, i] = simulate_openloop(ICs[:, i], P, T, Res)[0]
#         if ctrl==1:
#             Rho[:, :, i] = simulate_closedloop(ICs[:, i], P, delta, T, Res)[0]
#         plt.plot(Rho[0, :, i], Rho[1, :, i], 'k-')


#         arrow_base = (Rho[0, int((T*60*Res + 1)/40), i], Rho[1, int((T*60*Res + 1)/40), i])
#         arrow_head = (Rho[0, int((T*60*Res + 1)/40)+5, i], Rho[1, int((T*60*Res + 1)/40)+5, i])
#         # plt.plot(arrow_base[0], arrow_base[1], 'r*', markersize = 10)
#         # plt.plot(arrow_head[0], arrow_head[1], 'r*', markersize = 10)


#         plt.arrow(arrow_base[0], arrow_base[1], arrow_head[0] - arrow_base[0], arrow_head[1] - arrow_base[1], head_length=10, head_width=5, fc='k', ec='k')

#     # Mark the equilibrium points
#     if ctrl == 0:
#         if C_d < C:
#             if d < (1-epsilon)*C_d:
#                 plt.plot([d/vf], [d/vf], 'ro', markersize=8)
#             if d == (1-epsilon)*C_d:
#                 plt.plot([d/vf], [d/vf], 'ro', markersize=8)
#                 plt.plot([d/vf, d/vf], [C_d/vf, rho_j-d/w], 'r-', linewidth=3)
#                 plt.plot([d/vf, rho_j-d/w], [rho_j-d/w, rho_j-d/w], 'r-', linewidth=3)
#             if (1-epsilon)*C_d < d <= C_d:
#                 plt.plot([d/vf], [d/vf], 'ro', markersize=8)
#                 plt.plot([rho_j-(1-epsilon)*C_d/w], [rho_j-(1-epsilon)*C_d/w], 'r*', markersize=12)
#             if d > C_d:
#                 plt.plot([rho_j-(1-epsilon)*C_d/w], [rho_j-(1-epsilon)*C_d/w], 'ro', markersize=8)
#         else:
#             plt.plot([min(d, C)/vf], [min(d, C)/vf], 'ro', markersize=8)
#     else:
#         plt.plot([min(d, C_d)/vf], [min(d, C_d)/vf], 'ro', markersize=8)


#     plt.xlim((0, rho_j))
#     plt.ylim((0, rho_j))
#     plt.xlabel(r'$\rho_1$', fontsize=24)
#     plt.ylabel(r'$\rho_2$', fontsize=24)
#     subFig = phaseFig.add_subplot(111)
#     subFig.tick_params(axis='both', which='major', labelsize=18)
#     plt.show()


def curve_time(rho_0, n, P, delta, T, Res):
    '''
    Find the curver of rho, q and v vs. time of the system
    rho_0: n x 1 [double] initial condition
    n: Number of road sections
    P: (d, C_d, C, epsilon, vf, w, wbar, rho_j, rho_jbar)
    delta: (delta_1, delta_2)
    T: Termination time
    Res: number of simulation steps in each minute

    '''
    rho_ol, q_ol = simulate_openloop(rho_0, P, T, Res)
    rho_cl, q_cl, v = simulate_closedloop(rho_0, P, delta, T, Res)
    t = np.linspace(0, T, T*60*Res+1)
    rhoFig = plt.figure(figsize=(8.5,6.5))
    plt.plot(t*60, rho_cl[0, :], 'k-', linewidth = 2.5)
    plt.plot(t*60, rho_ol[0, :], 'k--', linewidth = 2.5)
    #plt.plot(t*60, rho[1, :], 'k--', linewidth = 2.5)
    plt.xlabel(r't (min)', fontsize = 24)
    plt.ylabel(r'density (veh/mi)', fontsize = 24)
    plt.legend([r'Closed-loop', r'Open-loop'], fontsize = 24, loc = 'best')
    subFig = rhoFig.add_subplot(111)
    subFig.tick_params(axis='both', which='major', labelsize=18)
    #plt.show()
    # rhoFig_cl = plt.figure(0, figsize=(8.5,6.5))
    # plt.plot(t*60, rho_cl[0, :], 'k-', linewidth = 2.5)
    # plt.xlabel(r't (min)', fontsize = 24)
    # plt.ylabel(r'density (veh/mi)', fontsize = 24)
    # #plt.legend([r'Closed-loop', r'Open-loop'], fontsize = 24, loc = 'best')
    # subFig = rhoFig_cl.add_subplot(111)
    # subFig.tick_params(axis='both', which='major', labelsize=18)

    qFig = plt.figure(figsize=(8.5,6.5))
    plt.plot(t*60, q_ol[0, :], 'k-', linewidth = 2.5)
    plt.plot(t*60, q_ol[1, :], 'k--', linewidth = 2.5)
    #plt.plot(t*60, q[2, :], 'k-.', linewidth = 2.5)
    plt.xlabel(r't (min)', fontsize = 24)
    plt.ylabel(r'Flow Rate (veh/hr)', fontsize = 24)
    plt.legend([r'$q_1$', r'$q_2$'], fontsize = 24, loc = 'best')
    plt.axis([0, T*60, 4200, 6100])
    subFig = qFig.add_subplot(111)
    subFig.tick_params(axis='both', which='major', labelsize=18)

    qFig_cl = plt.figure(figsize=(8.5,6.5))
    plt.plot(t*60, q_cl[0, :], 'k-', linewidth = 2.5)
    plt.plot(t*60, q_cl[1, :], 'k--', linewidth = 2.5)
    #plt.plot(t*60, q[2, :], 'k-.', linewidth = 2.5)
    plt.xlabel(r't (min)', fontsize = 24)
    plt.ylabel(r'Flow Rate (veh/hr)', fontsize = 24)
    plt.legend([r'$q_1$', r'$q_2$'], fontsize = 24, loc = 'best')
    subFig = qFig_cl.add_subplot(111)
    subFig.tick_params(axis='both', which='major', labelsize=18)   
   

    vFig = plt.figure(figsize=(8.5,6.5))
    plt.plot(t*60, v[0, :], 'k-', linewidth = 2.5)
    plt.xlabel(r't (min)', fontsize = 24)
    plt.ylabel(r'VSL (mi/hr)', fontsize = 24)
    #plt.legend([r'$v$'], fontsize = 24, loc='best')
    subFig = vFig.add_subplot(111)
    subFig.tick_params(axis='both', which='major', labelsize=18)
    plt.show()



    

    



def main():
    #d = 4000
    #d = 4420
    #d = 5000
    d = 6000
    C_d = 5200
    #C_d = 7000
    C = 6500
    epsilon = 0.15
    vf = 65
    w = 20
    wbar = 10
    rho_j = 425
    rho_jbar = 750
    P = (d, C_d, C, epsilon, vf, w, wbar, rho_j, rho_jbar)
    n = 2
    ctrl = 1
    if ctrl == 0:
        T = 0.5
    else:
        T = 16/60.0
    Res = 20
    delta = (20, 10)

    # phase_portrait(n, P, delta, T, Res)
    rho_0 = 110
    curve_time(rho_0, n, P, delta, T, Res)
    pass



if __name__ == '__main__':
    main()


    