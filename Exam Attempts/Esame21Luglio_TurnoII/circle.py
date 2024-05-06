####### FILE NON DATO ALL'ESEAME ########



import numpy as np
import scipy.linalg as spl
import RisolviSis
import matplotlib.pyplot as plt
import matplotlib.patches as patches

x1,y1 = 0,0
x2,y2 = 5,0
x3,y3 = 0,5
x4, y4 = 6,7

A = np.array([[x1,y1,1],[x2,y2,1],[x3,y3,1],[x4,y4,1]])
b = np.array([[-x1**2-y1**2],[-x2**2-y2**2],[-x3**2-y3**2],[-x4**2-y4**2]])

print(A)
print(b)

#Soluzione di un sistema sovradeterminato facendo uso della fattorizzazione QR    
def qrLS(A,b):
    n=A.shape[1]  # numero di colonne di A
    Q,R=spl.qr(A)
    h=Q.T@b
    x,flag=RisolviSis.Usolve(R[0:n,:],h[0:n])
    residuo=np.linalg.norm(h[n:])**2
    return x,residuo

x,res = qrLS(A,b)
print(x)
print(res)

r = np.sqrt(x[0]**2/4 + x[1]**2/4 - x[2])

fig, ax = plt.subplots()

# Crea una circonferenza con centro (h, k) e raggio r
circonferenza_sol = patches.Circle((-x[0]/2, -x[1]/2), r, fill=False)

# Aggiungi la circonferenza al grafico
ax.add_patch(circonferenza_sol)
ax.scatter([x1, x2, x3, x4], [y1, y2, y3, y4], color='red')


# Imposta gli assi per avere la stessa scala
ax.axis('equal')

# Mostra il grafico
plt.show()
