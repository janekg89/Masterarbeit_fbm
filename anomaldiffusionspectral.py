# -*- coding: utf-8 -*-
import numpy as np  # module for scientific computing
import matplotlib.pyplot as plt  # module for plotting "a la" matlab

D = 4
K_alpha = D
particles = 1000
n = 1000  # len of signal
alpha = 0.5

r_t_allparticles = []
v_ano_frq_all_particles = []
v_frq_all_particles = []


def Z(frq, alpha, K_alpha):
    Z = (((-1j * frq) ** (1 - alpha)) * K_alpha * np.math.gamma(1 + alpha))
    return Z


def compute_v_ano_t(v_t):
    v_t = np.array(v_t)
    v_frq = np.fft.fft(v_t)
    frq = np.fft.fftfreq(n) * 2 * np.pi
    # plt.plot(sin(frq))




    v_ano_frq = np.sqrt(Z(frq, alpha, K_alpha).real * 2) * v_frq
    # v_ano_frq=ones_like(v_ano_frq)
    # plt.plot(Z(frq,alpha,K_alpha).real)
    v_ano_frq[0] = np.random.normal(0, np.sqrt(2 * K_alpha * n ** alpha))
    # v_ano_frq[1]=np.random.normal(0,np.sqrt(2*K_alpha*(n-1)**alpha))

    v_ano_t = np.fft.ifft(v_ano_frq)

    return v_ano_t, v_ano_frq, v_frq, frq


def autocorr(x):
    result = np.correlate(x, x, mode='full')
    return result[result.size / 2:]


for particle in range(particles):
    v_t = []
    for i in range(n):
        v_t.append(np.random.normal(0, 1))


    # t = np.arange(n)
    # v_t=np.sin(t/(2*np.pi) )


    v_ano_t, v_ano_frq, v_frq, frq = compute_v_ano_t(v_t)

    r_t = np.cumsum(v_ano_t)  #Ort bei anomaler diffusion in abhängigkeit von der zeit
    r_t_allparticles.append(r_t)  #Ort bei anomaler diffusion fuer alle teilchen
    v_ano_frq_all_particles.append(v_ano_frq)
    v_frq_all_particles.append(v_frq)


def compute_msd(pathx):
    totalsize = len(pathx)
    msd = [None] * len(pathx)
    for j in range(0, totalsize):
        msdink = []
        for i in range(0, totalsize - j):
            msdink.append(((pathx[i] - pathx[j + i]) ** 2))
        msd[j] = sum(msdink) / float(totalsize - j)
    return msd


def msdanalyt(t, K_alpha, alpha):
    return 2 * ((t) ** alpha) * K_alpha


msd = [None] * particles
for particle in range(15):
    msd[particle] = compute_msd(r_t_allparticles[particle])
print len(msd[0])


def msd_ensemble(r_t_allparticles):
    r_t_allparticles = np.array(r_t_allparticles)
    r_t_allparticles_squared = r_t_allparticles ** 2
    msd = r_t_allparticles_squared.mean(axis=0)
    return msd


msd_ensemble = msd_ensemble(r_t_allparticles)


def caluculate_s_w(v_ano_frq_all_particles):
    v_ano_frq_all_particles = np.array(v_ano_frq_all_particles)
    v_ano_frq_all_particles_squared = abs(v_ano_frq_all_particles) ** 2
    # v_ano_frq_all_particles_squared=abs(v_ano_frq_all_particles_squared)
    s_w = v_ano_frq_all_particles_squared.mean(axis=0)
    return s_w


s_w_ano = caluculate_s_w(v_ano_frq_all_particles)
s_w1 = caluculate_s_w(v_frq_all_particles * np.sqrt(2 * D))

plt.figure(figsize=(20, 10))
# plt.xlim(0,100)
plt.plot(frq, s_w1[:] * np.sqrt(K_alpha * 2))
plt.plot(frq, s_w_ano[:])

colors = ['r', 'b', 'g', 'k', 'c', 'w', 'b', 'r', 'g', 'b', 'k', 'c', 'w', 'b', 'r', 'g', 'b', 'k', 'c', 'w', 'bo',
          'ro', 'go', 'bo', 'ko', 'co', 'wo', 'bo']
plt.figure(figsize=(20, 10))
plt.xlim(msd_ensemble.size - 40, msd_ensemble.size)

for particle in range(15):
    plt.plot(range(n), msd[particle], colors[particle], label="a. d. simulation")
plt.plot(range(msd_ensemble.size), msd_ensemble, colors[2], label="ensemble msd")
plt.plot(range(n), msdanalyt(np.array(range(n)), D, alpha), ":", color=colors[1], label="analytisch")
plt.legend(loc=3)
plt.xlabel('Steps', fontsize=14)
plt.ylabel('MSD', fontsize=14)
plt.figure(figsize=(50, 20))
