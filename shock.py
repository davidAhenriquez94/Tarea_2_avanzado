import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as sciopt

N = 101
infile = np.genfromtxt("datos.txt")

density = infile[:,0]
pression = infile[:,1]
velocity = infile[:,2]

#Exact solution                                                                                                                                                                                                                              
density_left = 1.0
density_right = 0.1
pressure_left = 1.0
pressure_right =0.1
velocity_left = 0.0
velocity_right = 0.0

def a(pressure, density):
    return np.sqrt((5.0/3)*(pressure/density))
def compatibility_equation(M):
    return M-(1.0/M)-a(pressure_left,density_left)*4*(1-((pressure_right/pressure_left)*((5.0/4)*(M**2)-(1.0/4)))**(1.0/5))
def Match_number_calculator(f,x0):
    return sciopt.newton_krylov(f,x0)[0]

pressure_1 = pressure_right*((5.0/4)*(Match_number_calculator(compatibility_equation,[1.0])**2))
density_1 = ((1.0/density_right)*((3.0/4)*(1.0/(Match_number_calculator(compatibility_equation,[1.0]))**2)+(1.0/4)))**(-1)
velocity_1 = (3.0/4)*(Match_number_calculator(compatibility_equation,[1.0])-(1.0/Match_number_calculator(compatibility_equation,[1.0])))
pressure_2  = pressure_1
velocity_2  = velocity_1
density_2  = density_left*((pressure_2/pressure_left)**(3.0/5))

def x1(t):
    return 0.5 - a(pressure_left,density_left)*t
def x2(t):
    return 0.5 - (velocity_2 - a(pressure_2,density_2))*t
def x3(t):
    return 0.5 + velocity_2*t
def x4(t):
    return 0.5 + ((Match_number_calculator(compatibility_equation,[1.0]))**2)*t

def velocity_expansion_fan(x,t):
    return (3.0/4)*(a(pressure_left, density_left)+(x-0.5)/t)
def a_expasion_fan(x,t):
    return a(pressure_left, density_left) - (1.0/3)*velocity_expansion_fan(x,t)
def pressure_expansion_fan(x,t):
    return pressure_left*((a_expasion_fan(x,t)/a(pressure_left,density_left))**5)
def density_expansion_fan(x,t):
    return ((density_2-density_left)/(x2(t)-x1(t))*x + (density_left*x2(t)-density_2*x1(t))/(x2(t)-x1(t)))

z1 = np.linspace(0,x1(0.13),25)
z2 = np.linspace(x1(0.13),x2(0.13),25)
z3  = np.linspace(x2(0.13),x3(0.13),25)
z4  = np.linspace(x3(0.13),x4(0.13),25)
z5 = np.linspace(x4(0.13),1,25)
y1 = np.ones(25)

y_pressure = np.append(pressure_left*y1,[pressure_expansion_fan(z1,2),pressure_2*y1,pressure_1*y1,pressure_right*y1])
y_density = np.append(density_left*y1,[density_expansion_fan(z1,2),density_2*y1,density_1*y1,density_right*y1])
y_velocity = np.append(velocity_left*y1,[velocity_expansion_fan(z1,2),velocity_2*y1,velocity_1*y1,velocity_right*y1])

z = np.append(z1,[z2,z3,z4,z5])

x = np.linspace(0,1,N) 


fig = plt.figure()
ax = plt.axes()
ax.set_xlabel("x")
ax.set_ylabel("Pressure")
ax.set_title("P vs x")
plt.plot(x,pression,'-o',color= "red", label = "Lax-Wendroff ")
plt.plot(z,y_pressure,label = "Exact solution")
plt.plot(())
ax.legend()
filename = 'pressure' 
plt.savefig(filename + '.pdf',format = 'pdf')
plt.close()

fig = plt.figure()
ax = plt.axes()
ax.set_xlabel("x")
ax.set_ylabel("Density")
ax.set_title("D vs x")
plt.plot(x,density,'-o',color= "red",label = "Lax-Wendroff")
plt.plot(z,y_density,label = "Exact solution")
ax.legend()
filename = 'Density'
plt.savefig(filename + '.pdf',format = 'pdf')
plt.close()

fig = plt.figure()
ax = plt.axes()
ax.set_xlabel("x")
ax.set_ylabel("Velocity")
ax.set_title("V vs x")
plt.plot(x,velocity,'-o',color= "red",label = "Lax-Wendroff")
plt.plot(z,y_velocity,label = "Exact solution")
ax.legend()
filename = 'velocity'
plt.savefig(filename + '.pdf',format = 'pdf')
plt.close()
