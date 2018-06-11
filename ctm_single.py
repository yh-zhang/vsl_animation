# formulation of the cell transmission model and the VSL controller for the multiple section network.
import numpy as np


def ctm_single_openloop(rho, t, P):
    '''
    The signle-section ctm open-loop model
    rho: double, vehicle density in the section
    t:   double, time
    P:   (d, C_d, C, epsilon, vf, w, wbar, rho_j, rho_jbar)
         d: upstream demand veh/hr
         C_d: downstream capacity veh/hr
         C: capacity of the section veh/hr
         epsilon: capacity drop coefficient
         vf: free flow speed mi/hr
         w: back propagation speed mi/hr
         wbar: outflow decreasing rate
         rho_j: jam density
         rho_jbar: jam density with outflow
    return 
    dRho:double, the time derivation of rho 
    '''

    q = q_openloop(rho, P)
    dRho = q[0] - q[1]
    return dRho


def ctm_single_closedloop(rho, t, P, delta, cMode):
    '''
    The ctm closed-loop model
    rho: double, vehicle density in the section
    t:   double, time
    P:   (d, C_d, C, epsilon, vf, w, wbar, rho_j, rho_jbar)
         d: upstream demand veh/hr
         C_d: downstream capacity veh/hr
         C: capacity of the section veh/hr
         epsilon: capacity drop coefficient
         vf: free flow speed mi/hr
         w: back propagation speed mi/hr
         wbar: outflow decreasing rate
         rho_j: jam density
         rho_jbar: jam density with outflow
    delta: (delta_1, delta_2) control constants
    return 
    dRho:double, the time derivation of rho 
    '''
    q = q_closedloop(rho, P, delta, cMode)[0]
    dRho = q[0] - q[1]

    return dRho

def q_openloop(rho, P):
    '''
    compute the flow rate for open loop

    '''
    d, C_d, C, epsilon, vf, w, wbar, rho_j, rho_jbar = P
    q = np.zeros(2)
    q[0] = min(d, C, w*(rho_j - rho))
    if C_d < C:
        if rho <= C_d/vf:
            q[1] = vf * rho
        else:
            q[1] = min(wbar*(rho_jbar - rho), (1-epsilon)*C_d)
    else:
        q[1] = min(vf * rho, wbar*(rho_jbar - rho))
    return q

def q_closedloop(rho, P, delta, cMode):
    '''
    Compute the flow rate and VSL control command for the closed-loop system
    rho: double, vehicle density in each section
    P:   (d, C_d, C, epsilon, vf, w, wbar, rho_j, rho_jbar)
    cMode: control mode of the switching logic. If cMode ==1, v[n-1] = v_{n-1,1}. If cMode==2, v[n-1] = v_{n-1,2}
    return 
    q: 2 x 1 [double], the flow rate in each section
    v: double, the speed limit in each section
    '''
    d, C_d, C, epsilon, vf, w, wbar, rho_j, rho_jbar = P
    # Controller constants
    delta_1, delta_2 = delta
    L = 40

    # flow rates
    q = np.zeros(2)

    # Error state
    x = rho - 0.999 * C_d / vf
     
    if rho <= C_d/vf:
        q[1] = vf * rho
    else:
        q[1] = min(wbar*(rho_jbar - rho), (1-epsilon)*C_d)

    if cMode == 1:
        v = w * (q[1] - L * (x + delta_1)) / (w*rho_j - (q[1] - L * (x + delta_1)))
    else:
        v = w * (q[1] - L * x) / (w*rho_j - (q[1] - L * x))
    v = np.median([0.1, 65, v])
    q[0] = min(d, v*w*rho_j/(v+w), C, w*(rho_j - rho))
    return q, v

# def vsl(rho, cMode, P):
#     '''
#     compute the VSL control command of the multiple section system
#     rho: n x 1 [double], vehicle density in each section
#     cMode: control mode of the switching logic. If cMode ==1, v[n-1] = v_{n-1,1}. If cMode==2, v[n-1] = v_{n-1,2}
#     P:   (d, C_d, C, epsilon, vf, w, wbar, rho_j, rho_jbar)
#     return
#     v:   (n-1) x 1 [double], the speed limit in each section
#     '''
#     n = len(rho)
#     d, C_d, C, epsilon, vf, w, wbar, rho_j, rho_jbar = P
#     x = rho - C_d / vf * np.ones(n)
#     q = 



def main():
    d = 3000
    C_d = 5200
    C = 6500
    epsilon = 0.15
    vf = 65
    w = 20
    wbar = 10
    rho_j = 425
    rho_jbar = 750
    P = (d, C_d, C, epsilon, vf, w, wbar, rho_j, rho_jbar)
    rho = 70
    dRho = ctm_single_openloop(rho, 0, P)
    print(dRho)
    

if __name__ == '__main__':
    main()
    
