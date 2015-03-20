# Assignment-6
Radioactivity Decay and Monte Carlo Integration Homework 

#Elizabeth Rodriguez
#Computational Physics
#Assignment 6


# 1: Radioactivity Decay
# This exercise looks at a more advanced version of the simple radioactive decay
#simulation in Example 10.1 from the book. The isotope 213Bi decays to stable 209Bi
#via one of two different routes, with probabilities and half-lives thus:
# 213Bi (46 min)----(2.09%)---> 209Ti (2.2 min)----> 209Pb (3.3 min)----> 209Bi
# 213Bi (46 min)----(97.91%)---> 209Pb (3.3 min)----> 209Bi
#(Technically, 209Bi isn't really stable, but it has a half-life of more than 10^19
#years, a billion times the age of the universe, so it might as well be). Starting 
#with a sample consisting of 10,000 atoms of 213Bi, simulate the decay of atoms as
#in Example 10.1 by dividing time into slices of length delta t = 1 s each and on 
#each step doing the following:


# (A): For each atom of 209Pb in turn, decide at random, with the appropriate probability,
#whether it decays or not. (The probability can be calculated from Eq. (10.3) in the book).
# Count the total number that decay, subtract it from the number of 209Pb atoms, and add 
#it to the number of 209Bi atoms.

import pylab as pl
import numpy as np

lam_Pb = 0.0035 # Decay Rate
dt = 1.0 # Time Step
s0 = s = 10000 # Initial number of nuclei
t0 = 0.0 # Initial Time
N = 2000 # Number of steps
 

# Initialize lists of results s, t, and e lists
slist = []
tlist = []
elist = []

for i in range(N): # Perform N Euler Steps
    e = s0*pl.exp(-lam_Pb*(t0)) # Compute Exact Solution
    d = 100.0*(abs(e-s)/e) # Get Percent Difference
    f = 4*'%7.3f' # Output Format
    slist.append(s) # Keep the s Value
    tlist.append(t0) # Keep the Time
    elist.append(e) # Keep Exact Solution
    s += -lam_Pb*s*dt # Compute the Change in s Using the Euler Algorithm
    t0 += dt # Increment the Time 

# Approximate Point at Which the Equation Goes From Linear to Exponential Decay  
e_1 = s0*pl.exp(-lam_Pb*(500))
print 'The approximate average decayed nuclei for Pb is:', round(e_1,3)
print 'The aprroximate number of undecayed nuclei for Pb at 500s is:', round(e_1-s,3) 

# Plotting and Showing Results    
pl.plot(tlist, slist, 'g.', tlist, elist, 'r-')
pl.xlabel('Time(s)')
pl.ylabel('Number of remaining nuclei')
pl.title('Radioactive Decay')
pl.legend(['Euler Algorithm', 'Exact Solution'])
pl.show()


# (B): Now do the same for 209Ti, except that decaying atoms are subtracted from the 
#total for 209Ti and added to the total for 209Pb.

import pylab as pl
import numpy as np

lam_Tl = 0.00525 # Decay Rate
dt = 4.0 # Time Step
s0 = s = 10000 # Initial number of nuclei
t0 = 0.0 # Initial Time
N = 500 # Number of steps
 

# Initialize lists of results s, t, and e lists
slist = []
tlist = []
elist = []

for i in range(N): # Perform N Euler Steps
    e = s0*pl.exp(-lam_Tl*(t0)) # Compute Exact Solution
    d = 100.0*(abs(e-s)/e) # Get Percent Difference
    f = 4*'%7.3f' # Output Format
    slist.append(s) # Keep the s Value
    tlist.append(t0) # Keep the Time
    elist.append(e) # Keep Exact Solution
    s += -lam_Tl*s*dt # Compute the Change in s Using the Euler Algorithm
    t0 += dt # Increment the Time 

# Approximate Point at Which the Equation Goes From Linear to Exponential Decay  
Pb_decayed = 1737.739
Pb_undecayed = 1728.732
e_1 = s0*pl.exp(-lam_Tl*(500))
c = Pb_decayed+e_1
print 'The approximate average decayed nuclei for Tl is:', round(e_1,3)
print 'The aprroximate number of undecayed nuclei for Tl at 500s is:', round(e_1-s,3)
print 'The approximate total of decayed nuclei is:', round(c,3)

# Plotting and Showing Results    
pl.plot(tlist, slist, 'g.', tlist, elist, 'r-')
pl.xlabel('Time(s)')
pl.ylabel('Number of remaining nuclei')
pl.title('Radioactive Decay')
pl.legend(['Euler Algorithm', 'Exact Solution'])
pl.show()


# (C): For 213Bi the situation is more complicated: when a 213Bi atom decays you have
#to decide at random with the appropriate probability the route by which it decays.
# Count the numbers that decay by each route and add and subtract accordingly.

import pylab as pl
import numpy as np

lam_213 = 2.19e-27 # Decay Rate
dt = 4.47e26 # Time Step
s0 = s = 10000 # Initial number of nuclei
t0 = 0.0 # Initial Time
N = 10000 # Number of steps
 

# Initialize lists of results s, t, and e lists
slist = []
tlist = []
elist = []

for i in range(N): # Perform N Euler Steps
    e = s0*pl.exp(-lam_213*(t0)) # Compute Exact Solution
    d = 100.0*(abs(e-s)/e) # Get Percent Difference
    f = 4*'%7.3f' # Output Format
    slist.append(s) # Keep the s Value
    tlist.append(t0) # Keep the Time
    elist.append(e) # Keep Exact Solution
    s += -lam_213*s*dt # Compute the Change in s Using the Euler Algorithm
    t0 += dt # Increment the Time 

# Approximate Point at Which the Equation Goes From Linear to Exponential Decay  
Pb_decayed = 1737.739
Pb_undecayed = 1728.732
Tl_decayed = 724.398
Tl_undecayed = 724.151
Total_Tl_Pb = 2462.137
e_1 = s0*pl.exp(-lam_213*(500))
c = e_1-Pb_decayed-Tl_decayed
print 'The approximate average decayed nuclei for Tl is:', round(e_1,3)
print 'The aprroximate number of undecayed nuclei for Tl at 500s is:', round(e_1-s,3)
print 'The approximate total of decayed nuclei is:', round(c,3)

# Plotting and Showing Results    
pl.plot(tlist, slist, 'g.', tlist, elist, 'r-')
pl.xlabel('Time(s)')
pl.ylabel('Number of remaining nuclei')
pl.title('Radioactive Decay')
pl.legend(['Euler Algorithm', 'Exact Solution'])
pl.show()
    
    
# Note that you have to work up the chain from the bottom like this, not down from the 
#top, to avoid inadvertently making the same atom decay twice on a single step. Keep
#track of the number of atoms each of the four isotopes at all times for 20,000 seconds
#and make a graph showing the four numbers as a function of time on the same axes.


# For full credit: Upload your program to you Github (make project for HW6) along with
#a text file (.txt) and a image file (.png or .eps) showing the graph it produces.



# 2: Monte Carlo Integration
# Calculate a value for the integral:
# I = int from 0-1 ((x^-1/2)/(e^x +1))dx 
#using the importance sampling formula, Eq. (10.45), with w(x) = x^1/2, as follows.


# (A): Show that the probability distribution p(x) from which the sample points should
#be drawn is given by:
# p(x) = 1/(2(sqrt(x)))
#and derive a transformation formula for generating random numbers between zero and 
#one from this distribution.

# In a Picture file

# (B): Using your formula, sample N = 1,000,000 random points and hence evaluate the
#integral. You should get a value around 0.84.

from scipy import integrate
from numpy import exp
from math import sqrt

a = 0.0
b = 1.0

def f(x):
    return((x**(-0.5))/(exp(x)+1))
    
c = integrate.quad(f, a, b)
print 'The estimated integral with its variance is:', c

# For full credit: Turn your derivations for part (a) and upload your program to your
#Github (in project HW6) along with a text file (.txt) showing the answers it calculates.

