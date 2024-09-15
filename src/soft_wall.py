import numpy as np
import matplotlib.pyplot as plt

class particle:
    def __init__(self):
        self.x=0
        self.y=0
        self.vx=0
        self.vy=0
        self.Fx=0
        self.Fy=0
        self.m=1



wallx=30
wally=0
wall_theta=30
wall_theta*=np.pi/180
k=1e-5

def wall_draw(wall_x, wall_y, theta):

    th=np.linspace(-10, 100, 100)
    
    x=wall_x+np.cos(theta)*th
    y=wall_y+np.sin(theta)*th

    plt.plot(x, y, color='k')
p1=particle()
p1.Fx=10
p1.Fy=0


dt=0.05
Nt=100

traj=np.zeros((Nt, 2))

for i in range(Nt):

    traj[i, 0]=p1.x
    traj[i, 1]=p1.y

    p1.vx+=dt*p1.Fx/p1.m
    p1.vy+=dt*p1.Fy/p1.m

    p1.x+=dt*p1.vx
    p1.y+=dt*p1.vy

    

plt.plot(traj[:, 0], traj[:, 1], lw=1, color='g') 
print (traj[-1, 0], traj[-1, 1])

p2=particle()
p2.Fx=p1.Fx
p2.Fy=p1.Fy

Nt=65000
traj=np.zeros((Nt, 2))
for i  in range(Nt):
    traj[i, 0]=p2.x
    traj[i, 1]=p2.y
    r=0.03
    th=np.random.uniform(0, 2*np.pi)
    dx=r*np.cos(th)
    dy=r*np.sin(th)

    d=(p2.y+dy-wally)*np.cos(wall_theta)-(p2.x+dx-wallx)*np.sin(wall_theta)
    energy_wall=0
    if d<0:
        continue
    else:
        energy_wall=k*(dx)**2+k*dy**2/2.
        p2.Fx=p2.Fx-k*(dx)
        p2.Fy=p2.Fy-k*dy
    


    energy=-p2.Fx*(dx)-p2.Fy*dy+energy_wall
    
    if energy<0:
        p2.x+=dx
        p2.y+=dy
    else:
        prob=np.exp(-energy)
        if prob > np.random.uniform(0, 1):
            p2.x+=dx
            p2.y+=dy
plt.plot(traj[:, 0], traj[:, 1], color='r', lw=0.5) 
#plt.xlim(0, 140)
#plt.ylim(-65, 65)

wall_draw(wallx, wally, wall_theta)
print (traj[-1, 0], traj[-1, 1])


plt.show()