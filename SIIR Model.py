#!/usr/bin/env python
# coding: utf-8

# In[15]:


from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# In[16]:


#Mean Vaccination % per day
df = pd.read_csv("daily-covid-19-vaccination-doses.csv")
v1 = df[df["Entity"]=="World"][14:]["new_vaccinations_smoothed"].mean()/8000000000
v2 = 1/365


# In[17]:


#SIIR Model for COVID-19 (d1 = 0.001)
def ChemRxnNet_ODE(t, y, N, beta, t1, t2):

    S,I1,I2,R1,R2,R3 = y

    b1 = 0.3/t1
    b2 = (1/t1)-b1
    c1 = 1/t2
    c2 = 0.8/t2
    c3 = (1/t2)-c2
    d1 = 0.001

    # differential equations
    dSdt = -beta*S*(I1+I2) + d1*(N-S-I1-I2-R1-R3)
    dI1dt = beta*S*(I1+I2) - b1*I1 - b2*I1
    dI2dt = b2*I1 - c1*I2
    dR1dt = b1*I1 - c2*R1 - c3*R1
    dR2dt = c1*I2 + c2*R1 - d1*R2
    dR3dt = c3*R1

    dydt = [dSdt,dI1dt,dI2dt,dR1dt,dR2dt,dR3dt]

    return dydt

# the time interval of the simulation (days)
tspan = [0,200]

# set initial conditions
y0 = [0.9996,0.0004,0,0,0,0] # S, I1, I2, R1, R2, R3

#Disease Modeling Parameters
N = 1
beta = 0.2
t1 = 5
t2 = 17

# integrate the ODE
sol = solve_ivp(ChemRxnNet_ODE, tspan, y0, args=(N, beta, t1, t2))

# plot the results
plt.plot(sol.t, sol.y[0], label='S')
plt.plot(sol.t, sol.y[1], label='I1')
plt.plot(sol.t, sol.y[2], label='I2')
plt.plot(sol.t, sol.y[3], label='R1')
plt.plot(sol.t, sol.y[4], label='R2')
plt.plot(sol.t, sol.y[5], label='R3')
plt.legend(loc="upper right")
plt.title("SIIR Model For COVID-19 (top row)")
plt.xlabel('Days')
plt.ylabel('% of Population')
plt.show()

print("% Dead = " + str(round(sol.y[5][-1]*100, 3)))


# In[18]:


#SIIR Model for COVID-19 (d = 0.003)
def ChemRxnNet_ODE(t, y, N, beta, t1, t2):

    S,I1,I2,R1,R2,R3 = y # unpack y

    b1 = 0.3/t1
    b2 = (1/t1)-b1
    c1 = 1/t2
    c2 = 0.8/t2
    c3 = (1/t2)-c2
    d1 = 0.003

    # differential equations
    dSdt = -beta*S*(I1+I2) + d1*(N-S-I1-I2-R1-R3)
    dI1dt = beta*S*(I1+I2) - b1*I1 - b2*I1
    dI2dt = b2*I1 - c1*I2
    dR1dt = b1*I1 - c2*R1 - c3*R1
    dR2dt = c1*I2 + c2*R1 - d1*R2
    dR3dt = c3*R1

    dydt = [dSdt,dI1dt,dI2dt,dR1dt,dR2dt,dR3dt] # repack dydt

    return dydt

# the time interval of the simulation
tspan = [0,500]

# set initial conditions
y0 = [0.9996,0.0004,0,0,0,0] # S, I1, I2, R1, R2, R3

#Disease Modeling Parameters
N = 1
beta = 0.2
t1 = 5
t2 = 17

# integrate the ODE
sol = solve_ivp(ChemRxnNet_ODE, tspan, y0, args=(N, beta, t1, t2))

# plot the results
plt.plot(sol.t, sol.y[0], label='S')
plt.plot(sol.t, sol.y[1], label='I1')
plt.plot(sol.t, sol.y[2], label='I2')
plt.plot(sol.t, sol.y[3], label='R1')
plt.plot(sol.t, sol.y[4], label='R2')
plt.plot(sol.t, sol.y[5], label='R3')
plt.legend(loc="upper right")
plt.title("SIIR Model For COVID-19 (bottom row)")
plt.xlabel('Days')
plt.ylabel('% of Population')
plt.show()

print("% Dead = " + str(round(sol.y[5][-1]*100, 3)))


# In[19]:


#SIIR Model for COVID-19 w/ Average Immunization Rate
def ChemRxnNet_ODE(t, y, N, beta, t1, t2):

    S,I1,I2,R1,R2,R3 = y # unpack y

    b1 = 0.3/t1
    b2 = (1/t1)-b1
    c1 = 1/t2
    c2 = 0.8/t2
    c3 = (1/t2)-c2
    d1 = 0.001

    # differential equations
    dSdt = -beta*S*(I1+I2) + d1*(N-S-I1-I2-R1-R3) - v1*S
    dI1dt = beta*S*(I1+I2) - b1*I1 - b2*I1
    dI2dt = b2*I1 - c1*I2
    dR1dt = b1*I1 - c2*R1 - c3*R1
    dR2dt = c1*I2 + c2*R1-d1*R2 + v1*S
    dR3dt = c3*R1

    dydt = [dSdt,dI1dt,dI2dt,dR1dt,dR2dt,dR3dt] # repack dydt

    return dydt

# the time interval of the simulation
tspan = [0,200]

# set initial conditions
y0 = [0.9996,0.0004,0,0,0,0] # S, I1, I2, R1, R2, R3

#Disease Modeling Parameters
N = 1
beta = 0.2
t1 = 5
t2 = 17

# integrate the ODE
sol = solve_ivp(ChemRxnNet_ODE, tspan, y0, args=(N, beta, t1, t2))

# plot the results
plt.plot(sol.t, sol.y[0], label='S')
plt.plot(sol.t, sol.y[1], label='I1')
plt.plot(sol.t, sol.y[2], label='I2')
plt.plot(sol.t, sol.y[3], label='R1')
plt.plot(sol.t, sol.y[4], label='R2')
plt.plot(sol.t, sol.y[5], label='R3')
plt.legend(loc="upper right")
plt.title("Vaccinated SIIR Model For COVID-19 (Average Rate)")
plt.xlabel('Days')
plt.ylabel('% of Population')
plt.show()

print("% Dead = " + str(round(sol.y[5][-1]*100, 3)))


# In[20]:


#SIIR Model for COVID-19 w/ Optimal Rate
def ChemRxnNet_ODE(t, y, N, beta, t1, t2):

    S,I1,I2,R1,R2,R3 = y # unpack y

    b1 = 0.3/t1
    b2 = (1/t1)-b1
    c1 = 1/t2
    c2 = 0.8/t2
    c3 = (1/t2)-c2
    d1 = 0.001

    # differential equations
    dSdt = -beta*S*(I1+I2) + d1*(N-S-I1-I2-R1-R3) - v2*S
    dI1dt = beta*S*(I1+I2) - b1*I1 - b2*I1
    dI2dt = b2*I1 - c1*I2
    dR1dt = b1*I1 - c2*R1 - c3*R1
    dR2dt = c1*I2 + c2*R1 - d1*R2 + v2*S
    dR3dt = c3*R1

    dydt = [dSdt,dI1dt,dI2dt,dR1dt,dR2dt,dR3dt] # repack dydt

    return dydt

# the time interval of the simulation
tspan = [0,200]

# set initial conditions
y0 = [0.9996,0.0004,0,0,0,0] # S, I1, I2, R1, R2, R3

#Disease Modeling Parameters
N = 1
beta = 0.2
t1 = 5
t2 = 17

# integrate the ODE
sol = solve_ivp(ChemRxnNet_ODE, tspan, y0, args=(N, beta, t1, t2))

# plot the results
plt.plot(sol.t, sol.y[0], label='S')
plt.plot(sol.t, sol.y[1], label='I1')
plt.plot(sol.t, sol.y[2], label='I2')
plt.plot(sol.t, sol.y[3], label='R1')
plt.plot(sol.t, sol.y[4], label='R2')
plt.plot(sol.t, sol.y[5], label='R3')
plt.legend(loc="upper right")
plt.title("Vaccinated SIIR Model For COVID-19 (Optimal Rate)")
plt.xlabel('Days')
plt.ylabel('% of Population')
plt.show()

print("% Dead = " + str(round(sol.y[5][-1]*100, 3)))


# In[ ]:




