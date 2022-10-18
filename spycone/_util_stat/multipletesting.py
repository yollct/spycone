import numpy as np 
import gc

class pvalue_correction():
    def __init__(self, pvals, method, alpha = 0.05):
        self.pvals= np.sort(pvals) # sort pvals in ascending order
        self.method = method
        self.alpha = alpha
        self.n = len(pvals)
        self.run()
        

    def run(self):
        METHODS = {'bonf': self.bonferroni, 'holm_bonf': self.holm_bonferroni, 'fdr_bh':self.fdr_bh}
        func = METHODS[self.method]
        func()
    
    def bonferroni(self):
        self.corrected_pvals  = self.pvals * self.n
        self.corrected_pvals[self.corrected_pvals>1]=1

    def holm_bonferroni(self):
        notreject = self.pvals >= (self.alpha / np.arange(self.n, 0, -1))
        pre_pvals_corrected = self.pvals * np.arange(self.n, 0, -1)
        self.corrected_pvals  = np.maximum.accumulate(pre_pvals_corrected)
        self.corrected_pvals[self.corrected_pvals>1]=1
        del pre_pvals_corrected

    def fdr_bh(self):
        conf = np.arange(self.n, 0,-1)/self.n
        pre_pvals_corrected = self.pvals * conf
        self.corrected_pvals  = np.maximum.accumulate(pre_pvals_corrected)
        self.corrected_pvals[self.corrected_pvals>1]=1
