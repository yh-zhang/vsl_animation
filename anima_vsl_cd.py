# Create the animation that shows points moving on fundamental diagram

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
import ctm_single as ctm
import simulation as sim

plt.rcParams['animation.ffmpeg_path'] = 'C:/ffmpeg/bin/ffmpeg'

def init_ol(P):
    qFun = ctm.q_openloop
    d, C_d, C, epsilon, vf, w, wbar, rho_j, rho_jbar = P
    rho = np.linspace(0, rho_j, 500)
    q_in = np.zeros((len(rho),))
    rho_1 = np.linspace(0, C_d/vf, 200)
    rho_2 = np.linspace(C_d/vf+1, rho_j, 300)
    q_out_1 = np.zeros((len(rho_1),))
    q_out_2 = np.zeros((len(rho_2),))
    for i in range(len(rho)):
        q_in[i] = qFun(rho[i], P)[0]
    for i in range(len(rho_1)):
        q_out_1[i] = qFun(rho_1[i], P)[-1]
    for i in range(len(rho_2)):
        q_out_2[i] = qFun(rho_2[i], P)[-1]

    # q_in = ctm.q_openloop(rho, P)[0]
    # q_out = ctm.q_openloop(rho, P)[-1]
    fig = plt.figure()
    ax = plt.axes(xlim = (0, rho_j+10), ylim = (0, C+500))
    line3, line4, line5, = ax.plot(rho_1, q_out_1, 'b-', rho_2, q_out_2, 'b-', rho, q_in, 'y-', linewidth=2)
    line1, = ax.plot([], [], 'ro', ms = 10)
    line2, = ax.plot([], [], 'ro', ms = 10)
    ax.legend((line1, line5, line3), ('Density', 'Inflow', 'Outflow'))
    return line1, line2, fig

def init_vsl(P):
    d, C_d, C, epsilon, vf, w, wbar, rho_j, rho_jbar = P
    rho = np.linspace(0, rho_j, 500)
    #q_in = np.zeros((len(rho),))
    # q_out = np.zeros((len(rho),))
    # for i in range(len(rho)):
    #     q_in[i] = min(d, C, w*(rho_j - rho[i]))
    #     q_out[i] = min(vf*rho[i], C_d - wbar*(rho[i] - C_d/vf))

    # fig = plt.figure()
    # ax = plt.axes(xlim = (0, rho_j+10), ylim = (0, C+100))
    # ax.plot(rho, q_out, 'b-', rho, q_in, 'y-', linewidth=2)
    qFun = ctm.q_openloop
    rho_1 = np.linspace(0, C_d/vf, 200)
    rho_2 = np.linspace(C_d/vf+1, rho_j, 300)
    q_out_1 = np.zeros((len(rho_1),))
    q_out_2 = np.zeros((len(rho_2),))
    # for i in range(len(rho)):
    #     q_in[i] = qFun(rho[i], P)[0]
    for i in range(len(rho_1)):
        q_out_1[i] = qFun(rho_1[i], P)[-1]
    for i in range(len(rho_2)):
        q_out_2[i] = qFun(rho_2[i], P)[-1]

    # q_in = ctm.q_openloop(rho, P)[0]
    # q_out = ctm.q_openloop(rho, P)[-1]
    fig = plt.figure()
    ax = plt.axes(xlim = (0, rho_j+10), ylim = (0, C+500))
    ax.plot(rho_1, q_out_1, 'b-', rho_2, q_out_2, 'b-', linewidth=2)
    # line1: 2 density dots; line2: curve of q_in
    line1, = ax.plot([], [], 'ro', ms = 10)
    line2, = ax.plot([], [], 'y-', linewidth=2)
    # line3, = ax.plot([], [], 'b--', linewidth=2)
    # line4, = ax.plot([], [], 'k--', linewidth=2)
    # ax.legend((line1, line2, line4, line3), ('Density', 'VSL', 'Inflow', 'Outflow'), loc = 1)

    #return line1, line2, line3, line4, fig
    return line1, line2, fig

def animate_ol(rho, line1, line2, P):
    qFun = ctm.q_openloop
    rho0 = rho[0][0]
    rho1 = rho[0][1]
    line1.set_data([rho0, rho0], [qFun(rho0, P)[-2], qFun(rho0, P)[-1]])
    line2.set_data([rho1, rho1], [qFun(rho1, P)[-2], qFun(rho1, P)[-1]])
    #line.set_data([rho], [ctm.q_openloop([rho], P)[-1]])
    return line1, line2

def animate_vsl(frame, line1, line2, P):
    rho = frame[0]
    q = frame[1:3]
    v = frame[3]
    d, C_d, C, epsilon, vf, w, wbar, rho_j, rho_jbar = P
    #q, v = ctm.q_closedloop(rho, P)

    # r1 = np.linspace(0, min(w*rho_j/(w+v), d/v), 500)
    # qr = np.zeros((len(r1),))
    # for i in range(len(r1)):
    #     #qr[i] = min(v*r[i], w*(rho_j - r[i]))
    #     qr[i] = v*r1[i]
    # r2 = np.linspace(0, rho_j, 500)
    # line1.set_data(rho, q[-1])
    # line2.set_data(r1, qr)
    # line3.set_data(r2, q[1])
    # line4.set_data(r2, q[0])
    r1 = np.linspace(0, rho_j, 1000)
    q_in = np.minimum(min(d, v*w*rho_j/(v+w)), w*(rho_j - r1))
    line1.set_data([rho, rho], [q[0], q[1]])
    line2.set_data(r1, q_in)

    return line1, line2

def animation_vsl(P, rho_0, T, delta, Res):
    d, C_d, C, epsilon, vf, w, wbar, rho_j, rho_jbar = P
    rho, q, v = sim.simulate_closedloop(rho_0, P, delta, T, Res)
    line1, line2, fig = init_vsl(P)
    if d < (1-epsilon)*C_d:
        omega = 1
    elif d == (1-epsilon)*C_d:
        omega = 2
    elif d < C_d:
        omega = 3
    else:
        omega = 4
    if C_d >= C:
        omega = 5
    fileName = 'vsl_%d.mp4'%omega
    #frames = np.transpose(rho).tolist()
    frames = np.vstack((rho, q, v))
    frames = np.transpose(frames).tolist()
    anim = animation.FuncAnimation(fig, animate_vsl, frames = frames, fargs = (line1, line2, P), interval = 20, blit = True, repeat=True)

    FFwriter = animation.FFMpegWriter(fps=30, extra_args=['-vcodec', 'libx264'])
    anim.save(fileName, writer=FFwriter)
    plt.show()

def animation_nc(P, rho_0, T, res):
    
    d, C_d, C, epsilon, vf, w, wbar, rho_j, rho_jbar = P

    line1, line2, fig = init_ol(P)
    if d < (1-epsilon)*C_d:
        omega = 1
    elif d == (1-epsilon)*C_d:
        omega = 2
    elif d < C_d:
        omega = 3
    else:
        omega = 4
    if C_d >= C:
        omega = 5
    fileName = 'nc_%d.mp4' % omega
    rho0, q0 = sim.simulate_openloop(0, P, T, res)
    rho1, q1 = sim.simulate_openloop(425, P, T, res)

    rho = [rho0, rho1]
    frames = np.transpose(rho).tolist()
    anim = animation.FuncAnimation(fig, animate_ol, frames = frames, fargs = (line1, line2, P), interval = 20, blit = True, repeat=False)

    #FFwriter = animation.FFMpegWriter(fps=60, extra_args=['-vcodec', 'libx264'])
    FFwriter = animation.FFMpegWriter(fps=30, extra_args=['-vcodec', 'libx264'])
    anim.save(fileName, writer=FFwriter)
    plt.show()


def main():
    ctrl = 'nc'
    #ctrl = 'vsl'
    #d = 3000
    #d = 4160
    #d = 5000
    d = 6000
    #C_d = 5200
    C_d = 7000
    C = 6500
    epsilon = 0.20
    vf = 65
    w = 20
    wbar = 3
    rho_j = 425
    rho_jbar = C/vf + C/wbar
    P = (d, C_d, C, epsilon, vf, w, wbar, rho_j, rho_jbar)
    rho_0 = rho_j
    #rho_0 = [0, 0]
    T = 0.3
    Res = 20
    delta = (20, 10)

    if ctrl == 'nc':
        animation_nc(P, rho_0, T, Res)
    else:
        animation_vsl(P, rho_0, T, delta, Res)
if __name__ == '__main__':
    main()