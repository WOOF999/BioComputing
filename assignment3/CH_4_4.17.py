import sympy as sy
import math 

nA_0=float(input("nA_O: "))
nB_0=float(input("nB_0: "))
nT_0=nA_0+nB_0
yA_0=nA_0/nT_0
yB_0=nB_0/nT_0

P=float(input("Pressure: "))
T=float(input("Temperature: "))

print(yA_0,yB_0)

Ke=4.79*10**(-13)*math.exp(11458/T)
x=sy.Symbol('x')

equation = sy.Eq((x*(nT_0-2*x)**2)/((nA_0-x)*(nB_0-2*x)**2)*((1/P)**2), Ke)

answer=sy.solve(equation)
print(answer)
