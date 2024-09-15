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


def find_distance(x1, y1, wth):
    r=np.sin(wth)*(x1-wx)-np.cos(wth)*(y1-wy)
    if r>0:
        return -r*np.sin(wth), r*np.cos(wth) 
    else:
        return 0, 0



wx=0
wy=0
wth=85*np.pi/180. 

t=np.linspace(-10, 10, 100)

wxx=wx+t*np.cos(wth)
wyy=wy+t*np.sin(wth)

WFx=np.zeros(100)
WFy=np.zeros(100)

p=particle()
p.x=-10
p.y=0

p.Fx=10
p.Fy=0

r=0.4
Nt=400
k=10
traj=np.zeros((Nt, 2))
for i in range(Nt):
    traj[i, 0]=p.x
    traj[i, 1]=p.y

    th=np.random.uniform(-np.pi, np.pi)
    dist1=find_distance(p.x, p.y, wth)
    dx=r*np.cos(th)
    dy=r*np.sin(th)

    p.x+=dx 
    p.y+=dy
    wFx, wFy=find_distance(p.x, p.y, wth)
    d=np.cos(wth)*(p.y-wy)-np.sin(wth)*(p.x-wx)
    #print ("wall_energy", 0.5*k*(dist2**2))
    print ("wall_energy", wFx*dx+wFy*dy)
    energy=-(p.Fx*dx+p.Fy*dy)-k*(wFx*dx+wFy*dy)
    if energy>0:
        if np.exp(-energy)<=np.random.uniform(0,1):
            p.x-=dx
            p.y-=dy


        

    plt.plot(traj[:i+1, 0], traj[:i+1, 1], 'r')
    plt.scatter(traj[i, 0], traj[i, 1], marker='+', color='b')
    plt.plot(wxx, wyy, 'k', marker='.')
    plt.pause(0.1)
    plt.cla()


'''
x1=-0.67
y1=-0.32

r=np.sin(wth)*(x1-wx)-np.cos(wth)*(y1-wy)
x0=x1-r*np.sin(wth)
y0=y1+r*np.cos(wth)
print (r)
plt.plot(wxx, wyy, 'k', marker='.')

plt.scatter(x1, y1, color='b')
plt.scatter(x0, y0, color='g')'''
plt.axis('equal')
plt.show()