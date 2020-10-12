import numpy as np
import matplotlib.pyplot as plt

def Fermi_func(A, Kd=False, DG = False):
	''' This funtion returns the characteristic sigmoidal function describing binding affinity in Ligand+Receptor framework '''
	if(Kd):
		return 1/(1+(Kd/A))
	if(DG):
		return 1/(1+(np.exp(DG)/A))
def Fermi_func2(Kd, A):
	''' This funtion returns the characteristic sigmoidal function describing binding affinity in Ligand+Receptor framework '''
	return 1/(1+(Kd/A))



Kd = 1e-9 #Molar
Kd2 = 1e-5
min_A = 1e-15
max_A = 1e-3
N_A = 6e23
A = np.logspace(np.log10(min_A), np.log10(max_A), 100)

scales = np.array(['linear', 'log'])
for yscale in scales:
	fig, ax = plt.subplots(figsize = (12, 8), gridspec_kw = {'top':0.9, 'bottom':0.15})
	#ax.set_title(r'$K_d = %.1e M$'%Kd, fontsize = 30)
	ax.plot(A, Fermi_func(A, Kd = Kd), linewidth = 4, color = 'darkblue', label =  r'$K_d = %.1e M$'%Kd)
	if (yscale == 'log'):
		ax.plot(A, Fermi_func(A, Kd = Kd2), linewidth = 4, color = 'darkred', label =  r'$K_d = %.1e M$'%Kd2)
		ax.vlines(1e8*1e3/N_A, 0, Fermi_func(1e8*1e3/N_A, Kd),color = 'grey', alpha = 0.8, linestyle = 'dashed')
		ax.hlines(Fermi_func(1e8*1e3/N_A, Kd), min_A, 1e8*1e3/N_A ,color = 'blue', alpha = 0.8, linestyle = 'dashed')
		ax.hlines(Fermi_func(1e8*1e3/N_A, Kd2), min_A, 1e8*1e3/N_A ,color = 'red', alpha = 0.8, linestyle = 'dashed')
		ax.set_title(r'Naive vs. Mature', fontsize = 30)
	if (yscale == 'linear'):
		ax.set_title(r'Mature', fontsize = 30)
	ax.hlines(1, min_A, max_A, linestyle = 'dashed')
	ax.hlines(0.5, min_A, max_A, linestyle = 'dashed')
	ax.vlines(Kd, 0, 1, linestyle = 'dashed')
	ax.set_xscale('log')
	ax.set_yscale(yscale)
	ax.set_xticklabels(['%.1e'%i for i in ax.get_xticks()*(N_A/1e3)])
	ax.set_xlabel('Antigen density [Antigens/mL]', fontsize = 30)
	ax.set_ylabel('Binding Probability', fontsize = 30)
	ax.tick_params(labelsize = 20)
	ax.legend(fontsize = 24, loc = 4)
	fig.savefig('../Figures/binding_probability_'+yscale+'.png')
	plt.show()

A_fix = 1e8*1e3/N_A
A_fix2 = 1e5*1e3/N_A
Kd_interval = np.logspace(-12, -4)
fig, ax = plt.subplots(figsize = (12, 8), gridspec_kw = {'top':0.9, 'bottom':0.15})
#ax.set_title(r'$K_d = %.1e M$'%Kd, fontsize = 30)
ax.plot(Kd_interval, 1/Fermi_func2(Kd_interval, A_fix), linewidth = 4, color = 'black', label =  r'$[A] = %.1e$ copies/mL'%(A_fix*N_A/1e3))
#ax.plot(Kd_interval, 1/Fermi_func2(Kd_interval, A_fix2), linewidth = 4, color = 'grey', label =  r'$[A] = %.1e$ copies/mL'%(A_fix2*N_A/1e3))
#ax.plot(Kd_interval, 1/Fermi_func(A_fix, Kd = Kd2), linewidth = 4, color = 'darkred', label =  r'$K_d = %.1e M$'%Kd2)
ax.vlines(Kd, 0, 1/Fermi_func2(Kd, A_fix), color = 'blue', alpha = 0.8, linestyle = 'dashed')
ax.vlines(Kd2, 0, 1/Fermi_func2(Kd2, A_fix), color = 'red', alpha = 0.8, linestyle = 'dashed')
ax.hlines(1/Fermi_func2(Kd, A_fix), 1e-12, Kd ,color = 'blue', alpha = 0.8, linestyle = 'dashed')
ax.hlines(1/Fermi_func2(Kd2, A_fix), 1e-12, Kd2 ,color = 'red', alpha = 0.8, linestyle = 'dashed')
#ax.hlines(1/Fermi_func2(Kd, A_fix2), 1e-12, Kd ,color = 'blue', alpha = 0.8, linestyle = 'dashed')
#ax.hlines(1/Fermi_func2(Kd2, A_fix2), 1e-12, Kd2 ,color = 'red', alpha = 0.8, linestyle = 'dashed')
ax.set_title(r'Naive vs. Mature', fontsize = 30)
#ax.hlines(1, min_A, max_A, linestyle = 'dashed')
#ax.hlines(0.5, min_A, max_A, linestyle = 'dashed')
#ax.vlines(Kd, 0, 1, linestyle = 'dashed')
ax.set_xscale('log')
ax.set_yscale(yscale)
#ax.set_xticklabels(['%.1e'%i for i in ax.get_xticks()*(N_A/1e3)])
ax.set_xlabel(r'$K_d [M]$', fontsize = 30)
ax.set_ylabel('Clone size', fontsize = 30)
ax.tick_params(labelsize = 20)
ax.legend(fontsize = 24, loc = 4)
fig.savefig('../Figures/clone_size_'+yscale+'.png')
plt.show()

