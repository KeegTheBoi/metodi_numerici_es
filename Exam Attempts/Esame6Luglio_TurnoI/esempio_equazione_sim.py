####### FILE NON DATO ALL'ESEAME ########


import sympy as sym
import sympy.utilities.lambdify as lambdify
x = sym.symbols('x')
b = sym.symbols('b')
c = sym.symbols('c')
# Delta = b**2 - 4ac => -b +- sqrt(b**2 - 2c)/2a => -b +- sqrt(b**2 - 2*c)

x1 = -b + sym.sign(b)*sym.sqrt(b**2 - 2*c)
pr = c*2
xs2 = pr/c

x1d = sym.diff(x1,c)
x2d = sym.diff(x2,c)

x1_b = x1.subs(b,10**8)
x1d_b = x1d.subs(b,10**8)
x2_b = x2.subs(b,10**8)
x2_b = x2d.subs(b,10**8)

i = np.arange(-5,10,1)
vet_c = 2.0**(-i)
bnum = 10**8

x1_num = lambdify(c,x1_b,np)
x1d_num = lambdify(c,x1d_b,np)

x2_num = lambdify(c,x2_b,np)
x2d_num = lambdify(x,x2d_b,np)

sol1 = x1_num(vet_c) 
sol2 = x2_num(vet_c)
cond1 = np.abs(x1_num(vet_c)*vet_c/x1d_num(vet_c))

print("Sol 1: ",sol1)
print("Condizionamento 1: ", cond1)
print("Sol 2: ",sol2)

print(sol1,sol2)
print(vet_c)
plt.plot(sol1,'o-',sol2,'o-')
plt.show()
