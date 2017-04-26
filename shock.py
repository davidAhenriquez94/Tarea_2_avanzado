import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as sciopt

N = 1001
infile = np.genfromtxt("datos.txt")

density = infile[:,0]
pression = infile[:,1]
velocity = infile[:,2]

#Exact solution                                                                                                                                                                                                                              
density_left = 1.0/0.1
density_right = 1.0
pressure_left = 1.0/(1.4*0.1)
pressure_right = 1.0/1.4
velocity_left = 0.0
velocity_right = 0.0

def a(pressure, density):
    return np.sqrt((1.4)*(pressure/density))
def compatibility_equation(M):
    return M-(1.0/M)-a(pressure_left,density_left)*6*(1-((pressure_right/pressure_left)*((7.0/6)*(M**2)-(1.0/6)))**(1.0/7))
def Match_number_calculator(f,x0):
    return sciopt.newton_krylov(f,x0)[0]

Ms = Match_number_calculator(compatibility_equation,[2.0])
pressure_1 = pressure_right*((7.0/6)*(Ms**2)-(1.0/6))
density_1 = ((1.0/density_right)*((5.0/6)*(1.0/(Ms**2))+(1.0/6)))**(-1)
velocity_1 = (5.0/6)*(Ms - 1.0/Ms)
pressure_2  = pressure_1
velocity_2  = velocity_1
density_2  = density_left*((pressure_2/pressure_left)**(5.0/7))

def x1(t):
    return 0.5 - a(pressure_left,density_left)*t
def x2(t):
    return 0.5 + (velocity_2 - a(pressure_2,density_2))*t
def x3(t):
    return 0.5 + velocity_2*t
def x4(t):
    return 0.5 + (Ms)*t

def velocity_expansion_fan(x,t):
    return (5.0/6)*(a(pressure_left, density_left)+((x-0.5)/t))
def a_expasion_fan(x,t):
    return a(pressure_left, density_left) - (0.2)*velocity_expansion_fan(x,t)
def pressure_expansion_fan(x,t):
    return pressure_left*((a_expasion_fan(x,t)/a(pressure_left,density_left))**7)
def density_expansion_fan(x,t):
    return ((a_expasion_fan(x,t)**2)/(1.4*pressure_expansion_fan(x,t)))**(-1)

t = 0.25
z1 = np.linspace(0,x1(t),25)
z2 = np.linspace(x1(t),x2(t),25)
z3  = np.linspace(x2(t),x3(t),25)
z4  = np.linspace(x3(t),x4(t),25)
z5 = np.linspace(x4(t),1,25)
y1 = np.ones(25)

y_pressure = np.append(pressure_left*y1,[pressure_expansion_fan(z2,t),pressure_2*y1,pressure_1*y1,pressure_right*y1])
y_density = np.append(density_left*y1,[density_expansion_fan(z2,t),density_2*y1,density_1*y1,density_right*y1])
y_velocity = np.append(velocity_left*y1,[velocity_expansion_fan(z2,t),velocity_2*y1,velocity_1*y1,velocity_right*y1])

z = np.append(z1,[z2,z3,z4,z5])

x = np.linspace(0,1,N) 


fig = plt.figure()
ax = plt.axes()
ax.set_xlabel("x/L")
ax.set_ylabel("Pressure (p/gamma*p_r)")
ax.set_title("P vs x")
plt.plot(x,pression,'-o',color= "red", label = "Lax-Wendroff Solution")
plt.plot(z,y_pressure,label = "Analytical Solution")
plt.plot(())
ax.legend()
filename = 'pressure' 
plt.savefig(filename + '.pdf',format = 'pdf')
plt.close()

fig = plt.figure()
ax = plt.axes()
ax.set_xlabel("x/L")
ax.set_ylabel("density (d/d_r)")
ax.set_title("D vs x")
plt.plot(x,density,'-o',color= "red",label = "Lax-Wendroff Solution")
plt.plot(z,y_density,label = "Analytical Solution")
ax.legend()
filename = 'density'
plt.savefig(filename + '.pdf',format = 'pdf')
plt.close()

fig = plt.figure()
ax = plt.axes()
ax.set_xlabel("x/L")
ax.set_ylabel("Velocity (v/a_r)")
ax.set_title("V vs x")
plt.plot(x,velocity,'-o',color= "red",label = "Lax-Wendroff Solution")
plt.plot(z,y_velocity,label = "Analytical Solution")
ax.legend()
filename = 'velocity'
plt.savefig(filename + '.pdf',format = 'pdf')
plt.close()
