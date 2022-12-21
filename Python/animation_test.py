import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time

def plot(network, artist):
    artist[0].set_data(network[:, 0], network[:, 1])

def next_step(pos, spd, rest, k, dt):
    # calculate pairwise distances
    dist = np.linalg.norm(pos[None, :, :] - pos[:, None, :], ord=2, axis=-1)
    dist[rest == 0] = 0 # remove distances where no spring attached
    mask = dist != 0 # mask to ensure nonzero division
    # calculate pairwise differences (along x and y)
    dist_x = pos[None, :, 0] - pos[:, None, 0]
    dist_y = pos[None, :, 1] - pos[:, None, 1]
    # forces matrix
    f_spring = (dist - rest)
    f_damp = np.zeros(f_spring.shape)
    f_tot = k * (f_damp + f_spring)
    # brodcast f on x and y
    dist_x[mask] /= dist[mask] # cos = adj / hyp
    dist_y[mask] /= dist[mask] # sin = opp / hyp
    fx = np.diag(f_tot.dot(dist_x.transpose()))
    fy = np.diag(f_tot.dot(dist_y.transpose()))
    f = np.array((fx, fy)).transpose()
    
    new_spd = spd.copy() + f * dt
    new_pos = pos.copy() + new_spd * dt

    return new_pos, new_spd

plt.ion()
# initial positions
net_pos = np.array((
    (0.25, 0.25),
    (0.25, 0.75),
    (0.75, 0.75),
    (0.75, 0.25)
), np.float32)
# initial speeds
net_spd = np.array((
    (0, 0),
    (0, 0),
    (0, 0),
    (0, 0)
), np.float32)
# rest length of springs
l0 = 0.6
rest_l = np.array((
    (0, l0, 0, l0),
    (l0, 0, l0, 0),
    (0, l0, 0, l0),
    (l0, 0, l0, 0),
), np.float32)
# stiffness
k = 0.1

fig, ax = plt.subplots()
artist = ax.plot(net_pos[:, 0], net_pos[:, 1], '.')
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.set_xticks([])
ax.set_yticks([])
ax.set_aspect('equal')

fps = 30
t_tot = 5
dt = 0.2
for i in range(t_tot * fps): # t_tot * fps
    net_pos, net_spd = next_step(net_pos, net_spd, rest_l, k, dt)
    plot(net_pos, artist)
    plt.pause(1/fps)

print('---- done ----')
plt.pause(2)