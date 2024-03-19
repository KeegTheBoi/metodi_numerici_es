import matplotlib.pyplot as plt
import numpy as np

relativerr = lambda f, g: np.abs(f - g) / np.abs(g)

class DomainCreator:
    
    def __init__(self, interval, step=1, BETA=10.0):
        self.beta = BETA    
        self.k = np.arange(interval[0],interval[1] + 1, step) 

    def compute(self, func):
        self.n = func(self.beta)


class FunctionEvaluator:
    
    def __init__(self, func, var, obsolete=None):
        self.func = func
        self.x = var
        self.y = self.func(var)
        self.yobs = obsolete(var) if obsolete != None else self.func(var)

    def plot_err(self, reals, compare=False, what=None, log=False):
        self.reals = reals
        x = self.x if what == None else what
        plt.loglog(x, relativerr(self.y, reals), 'r')
        if compare:
            plt.loglog(self.x, relativerr(self.yobs, reals), 'b:')
            plt.legend(["New", "Obsolete"])
        plt.title("Grafico dell'errore relativo di funzioni approssimate")
        plt.show()    

    def plot_compares(self, x, log=False):
        plt.plot(x, self.reals, 'b-', x, self.y, 'y-o')
        if log:
            plt.xscale('log')
        plt.legend(["Reals", "Machine"])
        plt.title("Confronto tra funzioni calcolate in aritmetica reale e macchina")