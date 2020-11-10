import numpy as np
import matplotlib.pyplot as plt
from Bio import Phylo
from io import StringIO
from matplotlib.lines import Line2D
from datetime import datetime, timedelta
import scipy.special as sc
import seaborn as sns
import pickle

class Sequence():
	"""docstring for Sequence"""
	def __init__(self, seq_id, master_sequence, energy_parent, complementary_sequence, Energy_Matrix, parent,  parent_id = 0, Master_Seq = False):
		super(Sequence, self).__init__()

		#A given Sequence object has the following attributes
		self.id = seq_id
		self.Alphabet = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't']
		self.parent = parent
		self.parent_id = parent_id
		self.master_sequence = master_sequence
		self.complementary_sequence = complementary_sequence
		self.active = False
		self.clone_size = 1
		self.energy = 0
		self.energy_parent = energy_parent
		self.delta_energy = 0
		self.sequence = parent
		self.pos_mut = 0
		self.tree_position = 1 # 0 = internal ; 1 = external
		self.hamming_distance = 0
		if not(Master_Seq):
			self.mutate_sequence(Energy_Matrix = Energy_Matrix)
		if (Master_Seq):
			self.energy = self.energy_parent

	def calculate_delta_energy(self, Energy_Matrix, old_letter, new_letter):

		M = Energy_Matrix

		list_comp_seq = list(self.complementary_sequence)

		pos_new_letter = np.where(np.isin(self.Alphabet,new_letter))[0][0]
		pos_old_letter = np.where(np.isin(self.Alphabet,old_letter))[0][0]

		pos_comp_seq = np.where(np.isin(self.Alphabet,list_comp_seq[self.pos_mut]))[0][0]

		self.delta_energy = M[pos_comp_seq][pos_new_letter]-M[pos_comp_seq][pos_old_letter]
		self.energy = self.energy_parent + self.delta_energy


	def mutate_sequence(self, Energy_Matrix):
		""" This function will create a new mutations and give an energy value to the new sequence according to the Energy_Matrix. """
		self.pos_mut = np.random.randint(9)
		list_seq = list(self.sequence)
		old_letter = self.sequence[self.pos_mut]
		self.Alphabet.remove(old_letter)
		new_letter = np.random.choice(self.Alphabet)
		list_seq[self.pos_mut] = new_letter
		self.sequence = "".join(list_seq)
		self.hamming_distance = hamming_distance(self.master_sequence, self.sequence)
		self.Alphabet = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't']
		self.calculate_delta_energy(Energy_Matrix = Energy_Matrix, old_letter = old_letter, new_letter = new_letter)
	
		#Ask Michael about the best way to produce the energy

class Stochastic_simulation():
	"""docstring for Stochastic_simulation"""
	def __init__(self, Sequences, n_linages, T, U, gamma, nu, R, beta, master_Sequence_energy):
		super(Stochastic_simulation, self).__init__()
		self.n_linages = n_linages
		self.Sequences = Sequences
		self.T = T
		self.U = U
		self.gamma = gamma
		self.nu = nu
		self.R = R
		self.beta = beta
		self.master_Sequence_energy = master_Sequence_energy

		self.linages_time_series = np.ones(shape =(n_linages, 1))
		self.antigen_time_series = np.array([1])
		self.time_series = np.array([0])
		self.probabilities = np.zeros((2*n_linages)+1)

	def calculate_probabilities(self):

		# Initialize with the event of antigen growth.
		self.probabilities[0] = self.beta*self.antigen_time_series[-1]

		# fill with Bcell activation and proliferation
		for i in range(1, self.n_linages+1):
			#Activation
			self.probabilities[(2*i)-1] = (self.antigen_time_series[-1]/(self.antigen_time_series[-1]+np.exp(self.master_Sequence_energy + self.Sequences[i-1].energy)))*(1-self.Sequences[i-1].active)
			#Proliferation
			self.probabilities[(2*i)] = self.nu*(self.linages_time_series[i-1,-1] - 1)

	def gillespie_step(self):

		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# 1. Generate 2 random numbers uniformly distributed in (0,1)
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		r1 = np.random.rand()
		r2 = np.random.rand()

		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# 2. Calculate probabilities
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

		self.calculate_probabilities()

		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# 3. Calculate alpha
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

		probabilities_cumsum = np.cumsum(self.probabilities)
		alpha = sum(self.probabilities)

		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# 4. Compute the time until the next event takes place
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		tau = (1/alpha)*np.log(float(1/r1))
		self.time_series = np.append(self.time_series, self.time_series[-1]+tau)

		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# 5. Compute which event takes place
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		transitionIdx   = np.searchsorted(probabilities_cumsum,r2*alpha)
		transition_Type  = transitionIdx % 2
		transition_Agent  = int((transitionIdx+1)/2)

		#antigen profileration event
		if(transition_Agent == 0):
			self.antigen_time_series = np.append(self.antigen_time_series, self.antigen_time_series[-1]+1)
			self.linages_time_series = np.hstack((self.linages_time_series, self.linages_time_series[:,-1].reshape(self.n_linages, 1)))
		else:
			#Activation event
			if(transition_Type==1):
				self.Sequences[transition_Agent-1].active = True
				self.antigen_time_series = np.append(self.antigen_time_series, self.antigen_time_series[-1])
				temp_array = self.linages_time_series[:,-1].reshape(self.n_linages, 1)
				temp_array[transition_Agent-1] = temp_array[transition_Agent-1] + 1
				self.linages_time_series = np.hstack((self.linages_time_series, self.linages_time_series[:,-1].reshape(self.n_linages, 1)))
			#B cell proliferation event
			else:
				self.antigen_time_series = np.append(self.antigen_time_series, self.antigen_time_series[-1])
				temp_array = self.linages_time_series[:,-1].reshape(self.n_linages, 1)
				temp_array[transition_Agent-1] = temp_array[transition_Agent-1] + 1
				self.linages_time_series = np.hstack((self.linages_time_series, temp_array))

	def Gillespie(self):

		while(self.time_series[-1] < self.T):
			self.gillespie_step()

	def plot_energy_matrix(self, ax):

		M = [[-1.06,0.19,-0.23,0.16,-0.08,0.06,0.08,0.04,0.00,-0.08,0.19,-0.02,0.05,0.13,0.69,0.03,-0.19,0.24,0.71,0.00],
	 	[0.19,0.04,-0.42,-0.28,-0.20,-0.14,-0.67,-0.13,0.25,0.19,0.19,0.14,0.46,0.08,0.44,0.65,0.99,0.31,0.00,-0.34],
		[-0.23,-0.42,-0.44,-0.19,-0.30,-0.2,-0.16,0.00,0.03,0.38,0.31,0.29,0.49,0.18,0.27,0.39,-0.16,0.41,0.44,0.20],
		[0.16,-0.28,-0.19,-0.22,-0.41,-0.25,0.02,0.11,-0.22,0.25,0.14,0.21,0.36,0.53,0.35,0.59,0.49,0.42,0.36,0.25],
		[-0.08,-0.20,-0.30,-0.41,-0.27,-0.29,-0.09,0.24,-0.01,0.23,0.20,0.25,0.26,0.30,0.43,0.67,0.16,0.35,0.19,0.42],
		[0.06,-0.14,-0.22,-0.25,-0.29,-0.29,-0.07,0.02,-0.10,0.16,0.25,0.18,0.24,0.50,0.34,0.58,0.19,0.30,0.44,0.09],
		[0.08,-0.67,-0.16,0.02,-0.09,-0.07,-0.12,-0.04,-0.09,-0.18,0.22,0.34,0.08,0.06,0.29,0.24,-0.12,-0.16,0.22,-0.28],
		[0.04,-0.13,0.00,0.11,0.24,0.02,-0.04,-0.06,0.09,0.14,0.13,0.09,-0.20,-0.20,-0.10,0.00,-0.34,-0.25,-0.21,-0.33],
		[0.00,0.25,0.03,-0.22,-0.01,-0.10,-0.09,0.09,-0.13,-0.07,-0.09,-0.06,0.08,0.28,0.26,0.12,0.34,0.43,0.14,0.10],
		[-0.08,0.19,0.38,0.25,0.23,0.16,0.18,0.14,-0.07,-0.38,-0.26,-0.16,-0.06,-0.14,0.25,-0.22,0.20,-0.04,0.11,-0.11],
		[0.19,0.19,0.31,0.14,0.20,0.25,0.22,0.13,-0.09,-0.26,0.03,-0.08,-0.14,-0.11,0.00,-0.29,-0.19,-0.35,-0.09,-0.07],
		[-0.02,0.14,0.29,0.21,0.25,0.18,0.34,0.09,-0.06,-0.16,-0.08,-0.20,-0.14,-0.14,0.26,-0.31,-0.05,0.17,-0.13,0.01],
		[0.05,0.46,0.49,0.36,0.26,0.24,0.08,-0.20,0.08,-0.06,-0.14,-0.14,0.29,-0.25,-0.17,-0.17,-0.02,-0.52,-0.38,-0.42],
		[0.13,0.08,0.18,0.53,0.30,0.50,0.06,-0.20,0.28,-0.14,-0.11,-0.14,-0.25,-0.53,-0.32,-0.30,-0.24,-0.14,-0.33,-0.18],
		[0.69,0.44,0.27,0.35,0.43,0.34,0.29,-0.10,0.26,0.25,0.00,-0.26,-0.17,-0.32,-0.03,-0.15,-0.45,-0.74,-0.97,-0.10],
		[0.03,0.65,0.39,0.59,0.67,0.58,0.24,0.00,0.12,-0.22,-0.29,-0.31,-0.17,-0.30,-0.15,0.04,-0.39,-0.72,-0.76,0.04],
		[-0.19,0.99,-0.16,0.49,0.16,0.19,-0.12,-0.34,0.34,0.20,-0.19,-0.05,-0.02,-0.24,-0.45,-0.39,-0.29,-0.12,0.22,-0.21],
		[0.24,0.31,0.41,0.42,0.35,0.30,-0.16,-0.25,0.43,-0.04,-0.35,0.17,-0.52,-0.14,-0.74,-0.72,-0.12,0.11,0.75,-0.38],
		[0.71,0.00,0.44,0.36,0.19,0.44,0.22,-0.21,0.14,0.11,-0.09,-0.13,-0.38,-0.33,-0.97,-0.76,0.22,0.75,0.25,0.11],
		[0.00,-0.34,0.20,0.25,0.42,0.09,-0.28,-0.33,0.10,-0.11,-0.07,0.01,-0.42,-0.18,-0.10,0.04,-0.21,-0.38,0.11,0.26]]
		sns.heatmap(np.flip(M, axis = 0), ax = ax, cmap=plt.cm.RdBu_r)
		ax.set_title('MJ-Matrix', fontsize = 22)
		ax.tick_params(labelsize = 20)
		ax.set_xticklabels(['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't'])
		ax.set_yticklabels(np.flip(['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't']));
		cbar = ax.collections[0].colorbar
		cbar.ax.tick_params(labelsize=18)

	def plot_antigen_time(self, ax):

		ax.plot(self.time_series, self.antigen_time_series, linewidth  = 4)
		ax.set_yscale('log')
		ax.set_xlabel(r'Time $t$', fontsize = 20)
		ax.set_ylabel(r'Antigen $\rho$', fontsize = 20)
		ax.tick_params(labelsize = 20)
		#handles, labels = ax.get_legend_handles_labels()
		#ax.legend(np.concatenate(([Line2D([0], [0], color='tab:red', linewidth=4, linestyle='solid', ms = 8)],handles)),np.concatenate(([r'$n_b(r, \rho)$'],labels)), loc = 0, fontsize = 20)

	def plot_prob_binding(self, ax):
		rho_array = np.logspace(0, np.log10(max(self.antigen_time_series)), 5)
		colors = plt.cm.Reds(np.linspace(0,1,len(rho_array)))
		for i, rho in enumerate(rho_array): 
			ax.plot(np.linspace(-2,8,10), (1/(1+np.exp(self.master_Sequence_energy + np.linspace(-2,8,10) - np.log(rho)))), linewidth  = 4, color = colors[i], label = r'$\rho = %.0e$'%(rho))
		ax.set_yscale('log')
		ax.set_xlabel(r'Energy $\epsilon$', fontsize = 20)
		ax.set_ylabel(r'Probability of binding $p_b$', fontsize = 20)
		ax.tick_params(labelsize = 20)
		ax.legend(loc = 0, fontsize = 20)

	def stackplot_linages_time(self, ax):
		colors = []
		for i in self.Sequences:
			if(i.active==True):
				colors.append('tab:red')
			else:
				colors.append('indigo')
		ax.stackplot(self.time_series, self.linages_time_series/np.sum(self.linages_time_series, axis = 0), colors=colors, alpha = 0.9);
		#ax.set_xscale('log')
		#ax.set_yscale('log')
		ax.set_xlabel(r'Time $t$', fontsize = 20)
		ax.set_ylabel(r'B cell Linages', fontsize = 20)
		ax.tick_params(labelsize = 20)
		#ax.legend(loc = 0, fontsize = 20)

	def hist_sequences_hamming_distance(self, Sequences, ax):

		rho_array = np.logspace(0, np.log10(max(self.antigen_time_series)), 5)
		colors = plt.cm.Reds(np.linspace(0,1,len(rho_array)))
		data_distances = ax.hist([Sequences[i].hamming_distance for i in range(int(len(Sequences)))], bins = range(10), align = 'left', label = r'$S(d)$', color = 'lightsteelblue', alpha = 0.5)
		ax.plot(data_distances[1][0:-1], sc.comb(9, data_distances[1][0:-1])*((20-1)**data_distances[1][0:-1]), linewidth = 4 , color = 'lightsteelblue', alpha = 0.6)

		ax.hist([self.Sequences[i].hamming_distance for i in range(int(len(self.Sequences)))], bins = range(10), align = 'left', label = r'$US(d)$', color = 'indigo', alpha = 0.6)
		ax.plot(data_distances[1][0:-1], self.U*sc.comb(9, data_distances[1][0:-1])*((20-1)**data_distances[1][0:-1]), linewidth = 4 , color = 'indigo', alpha = 0.6)

		ax.hist([self.Sequences[i].hamming_distance for i in range(int(len(self.Sequences))) if self.Sequences[i].active], bins = range(10), align = 'left', label = r'Activated Linages', color = 'tab:red', alpha = 0.8)
		for i, rho in enumerate(rho_array):
			ax.plot(data_distances[1][0:-1], self.U*sc.comb(9, data_distances[1][0:-1].astype(int))*((20-1)**data_distances[1][0:-1])*(1/(1+np.exp(self.master_Sequence_energy + data_distances[1][0:-1] - np.log(rho)))) , color = colors[i], linestyle = 'dashed', linewidth = 3)

		ax.set_ylim(0.1, 2e5)    
		ax.set_yscale('log')
		ax.set_xlabel(r'Hamming Distance $d$', fontsize = 20)
		ax.set_ylabel(r'Number of linages', fontsize = 20)
		ax.tick_params(labelsize = 20)
		ax.legend(loc = 0, fontsize = 20)

	def hist_sequences_energy(self, Sequences, ax):

		rho_array = np.logspace(0, np.log10(max(self.antigen_time_series)), 5)
		colors = plt.cm.Reds(np.linspace(0,1,len(rho_array)))
		data_energies = ax.hist([Sequences[i].energy for i in range(int(len(Sequences)))], bins = 20, align = 'left', label = r'$S(\epsilon)$', color = 'lightsteelblue', alpha = 0.5)
		#ax.plot(data_energies[1][0:-1], sc.comb(9, data_energies[1][0:-1])*((20-1)**data_energies[1][0:-1]), linewidth = 4 , color = 'lightsteelblue', alpha = 0.6)

		ax.hist([self.Sequences[i].energy for i in range(int(len(self.Sequences)))], bins = 20, align = 'left', label = r'$US(\epsilon)$', color = 'indigo', alpha = 0.6)
		#ax.plot(data_energies[1][0:-1], self.U*sc.comb(9, data_energies[1][0:-1])*((20-1)**data_energies[1][0:-1]), linewidth = 4 , color = 'indigo', alpha = 0.6)

		ax.hist([self.Sequences[i].energy for i in range(int(len(self.Sequences))) if self.Sequences[i].active], bins = 20, align = 'left', label = r'Activated Linages', color = 'tab:red', alpha = 0.8)
		#for i, rho in enumerate(rho_array):
		#	ax.plot(data_energies[1][0:-1], self.U*sc.comb(9, data_energies[1][0:-1].astype(int))*((20-1)**data_energies[1][0:-1])*(1/(1+np.exp(self.master_Sequence_energy + data_energies[1][0:-1] - np.log(rho)))) , color = colors[i], linestyle = 'dashed', linewidth = 3)

		#ax.set_ylim(0.1, 2e5)    
		ax.set_yscale('log')
		ax.set_xlabel(r'Energy $\epsilon$', fontsize = 20)
		ax.set_ylabel(r'Number of linages', fontsize = 20)
		ax.tick_params(labelsize = 20)
		ax.legend(loc = 0, fontsize = 20)

	def plot_k_largest_linages(self, k, ax):

		#Calculate array of the frequencies of the largest k linages
		Seq_states = [i.active for i in self.Sequences]
		Seq_sizes = self.linages_time_series[:,-1]
		k = k
		biggest_k_linages_sizes = np.sort(Seq_sizes)[-k:]
		Pos = np.array([i for i, j in enumerate(Seq_sizes) if np.isin(j,biggest_k_linages_sizes)])
		biggest_k_linages = self.linages_time_series[Pos,:]
		#for i in range(1,int(len(self.linages_time_series[0,:]))):
		#    biggest_k_linages = np.vstack((biggest_k_linages, self.linages_time_series[Pos,i]))
		#biggest_k_linages_freq = np.transpose(biggest_k_linages)/np.sum(np.transpose(biggest_k_linages), axis = 0)
		biggest_k_linages_freq = biggest_k_linages/np.sum(biggest_k_linages, axis = 0)

		ax.stackplot(self.time_series, biggest_k_linages_freq, alpha = 0.9);
		#ax[0].set_yscale('log')
		ax.set_xlabel(r'Time $t$', fontsize = 20)
		ax.set_ylabel(r'k largest linages', fontsize = 20)
		ax.tick_params(labelsize = 20)

		return biggest_k_linages_freq

	def plot_entropy_k_largest_linages(self, k, biggest_k_linages_freq, ax):

		#Calculate array of the frequencies of the largest k linages
		#Seq_states = [i.active for i in self.Sequences]
		#Seq_sizes = self.linages_time_series[:,-1]
		#k = k
		#biggest_k_linages_sizes = np.sort(Seq_sizes)[-k:]
		#Pos = np.array([i for i, j in enumerate(Seq_sizes) if np.isin(j,biggest_k_linages_sizes)])
		#biggest_k_linages = self.linages_time_series[Pos,:]
		#for i in range(1,int(len(self.linages_time_series[0,:]))):
		#    biggest_k_linages = np.vstack((biggest_k_linages, self.linages_time_series[Pos,i]))
		#biggest_k_linages_freq = np.transpose(biggest_k_linages)/np.sum(np.transpose(biggest_k_linages), axis = 0)
		#biggest_k_linages_freq = biggest_k_linages/np.sum(biggest_k_linages, axis = 0)

		#Calculate entropy
		entropy = [np.sum(-1*biggest_k_linages_freq[:,t]*np.log(biggest_k_linages_freq[:,t])) for t in range(int(len(self.time_series)))]
		ax.plot(self.time_series, entropy, linewidth = '4', color = 'indigo')
		#ax[1].set_yscale('log')
		ax.set_xlabel(r'Time $t$', fontsize = 20)
		ax.set_ylabel(r'Entropy', fontsize = 20)
		ax.tick_params(labelsize = 20)
		
def print_raw_file(Sequences, filename):

	file = open(filename, 'w+')
	for i in range(len(Sequences)):
		np.savetxt(file, np.array([Sequences[i].parent_id, Sequences[i].id]), fmt = '%d', delimiter = ' ', newline = ' ')
		file.write("\n")
	file.close()

def generate_newick_format(filename):

	file  = np.loadtxt(filename, dtype = int)
	n_f = '0()'
	print(file)

def hamming_distance(chaine1, chaine2):

    return sum(c1 != c2 for c1, c2 in zip(chaine1, chaine2))

def find_complementary_seq(sequence, Energy_Matrix):

	M = Energy_Matrix
	Alphabet = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't']
	list_seq = list(sequence)
	list_comp_seq = []
	for i in list_seq:
		pos_i = np.where(np.isin(Alphabet,i))[0][0]
		list_comp_seq.append(Alphabet[np.where(np.isin(M[pos_i],min(M[pos_i])))[0][0]])
	comp_seq = "".join(list_comp_seq)
	return comp_seq

def calculate_energy(Energy_Matrix, seq1, seq2):

	M = Energy_Matrix
	Alphabet = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't']
	list_seq1 = list(seq1)
	list_seq2 = list(seq2)
	Energy = 0
	for i in range(9):
		pos_i = np.where(np.isin(Alphabet,list_seq1[i]))[0][0]
		pos_j = np.where(np.isin(Alphabet,list_seq2[i]))[0][0]
		Energy += M[pos_i][pos_j]
	return Energy

def generate_Sequences(n_seq):

	M = [[-1.06,0.19,-0.23,0.16,-0.08,0.06,0.08,0.04,0.00,-0.08,0.19,-0.02,0.05,0.13,0.69,0.03,-0.19,0.24,0.71,0.00],
 	[0.19,0.04,-0.42,-0.28,-0.20,-0.14,-0.67,-0.13,0.25,0.19,0.19,0.14,0.46,0.08,0.44,0.65,0.99,0.31,0.00,-0.34],
	[-0.23,-0.42,-0.44,-0.19,-0.30,-0.2,-0.16,0.00,0.03,0.38,0.31,0.29,0.49,0.18,0.27,0.39,-0.16,0.41,0.44,0.20],
	[0.16,-0.28,-0.19,-0.22,-0.41,-0.25,0.02,0.11,-0.22,0.25,0.14,0.21,0.36,0.53,0.35,0.59,0.49,0.42,0.36,0.25],
	[-0.08,-0.20,-0.30,-0.41,-0.27,-0.29,-0.09,0.24,-0.01,0.23,0.20,0.25,0.26,0.30,0.43,0.67,0.16,0.35,0.19,0.42],
	[0.06,-0.14,-0.22,-0.25,-0.29,-0.29,-0.07,0.02,-0.10,0.16,0.25,0.18,0.24,0.50,0.34,0.58,0.19,0.30,0.44,0.09],
	[0.08,-0.67,-0.16,0.02,-0.09,-0.07,-0.12,-0.04,-0.09,-0.18,0.22,0.34,0.08,0.06,0.29,0.24,-0.12,-0.16,0.22,-0.28],
	[0.04,-0.13,0.00,0.11,0.24,0.02,-0.04,-0.06,0.09,0.14,0.13,0.09,-0.20,-0.20,-0.10,0.00,-0.34,-0.25,-0.21,-0.33],
	[0.00,0.25,0.03,-0.22,-0.01,-0.10,-0.09,0.09,-0.13,-0.07,-0.09,-0.06,0.08,0.28,0.26,0.12,0.34,0.43,0.14,0.10],
	[-0.08,0.19,0.38,0.25,0.23,0.16,0.18,0.14,-0.07,-0.38,-0.26,-0.16,-0.06,-0.14,0.25,-0.22,0.20,-0.04,0.11,-0.11],
	[0.19,0.19,0.31,0.14,0.20,0.25,0.22,0.13,-0.09,-0.26,0.03,-0.08,-0.14,-0.11,0.00,-0.29,-0.19,-0.35,-0.09,-0.07],
	[-0.02,0.14,0.29,0.21,0.25,0.18,0.34,0.09,-0.06,-0.16,-0.08,-0.20,-0.14,-0.14,0.26,-0.31,-0.05,0.17,-0.13,0.01],
	[0.05,0.46,0.49,0.36,0.26,0.24,0.08,-0.20,0.08,-0.06,-0.14,-0.14,0.29,-0.25,-0.17,-0.17,-0.02,-0.52,-0.38,-0.42],
	[0.13,0.08,0.18,0.53,0.30,0.50,0.06,-0.20,0.28,-0.14,-0.11,-0.14,-0.25,-0.53,-0.32,-0.30,-0.24,-0.14,-0.33,-0.18],
	[0.69,0.44,0.27,0.35,0.43,0.34,0.29,-0.10,0.26,0.25,0.00,-0.26,-0.17,-0.32,-0.03,-0.15,-0.45,-0.74,-0.97,-0.10],
	[0.03,0.65,0.39,0.59,0.67,0.58,0.24,0.00,0.12,-0.22,-0.29,-0.31,-0.17,-0.30,-0.15,0.04,-0.39,-0.72,-0.76,0.04],
	[-0.19,0.99,-0.16,0.49,0.16,0.19,-0.12,-0.34,0.34,0.20,-0.19,-0.05,-0.02,-0.24,-0.45,-0.39,-0.29,-0.12,0.22,-0.21],
	[0.24,0.31,0.41,0.42,0.35,0.30,-0.16,-0.25,0.43,-0.04,-0.35,0.17,-0.52,-0.14,-0.74,-0.72,-0.12,0.11,0.75,-0.38],
	[0.71,0.00,0.44,0.36,0.19,0.44,0.22,-0.21,0.14,0.11,-0.09,-0.13,-0.38,-0.33,-0.97,-0.76,0.22,0.75,0.25,0.11],
	[0.00,-0.34,0.20,0.25,0.42,0.09,-0.28,-0.33,0.10,-0.11,-0.07,0.01,-0.42,-0.18,-0.10,0.04,-0.21,-0.38,0.11,0.26]]

	Alphabet = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't']
	antigen_sequence = "".join(np.random.choice(Alphabet, 9))
	print('Antigen Seq: ' + antigen_sequence + '\n')

	master_sequence = find_complementary_seq(sequence = antigen_sequence , Energy_Matrix =  M)
	print('Master Seq: ' + master_sequence + '\n')

	master_sequence_energy = calculate_energy(Energy_Matrix = M, seq1 = master_sequence, seq2 = antigen_sequence)

	print('Master Seq energy: ', master_sequence_energy, '\n')

	Master_Sequence = Sequence(seq_id = 0, parent = master_sequence, energy_parent = master_sequence_energy, Master_Seq=True, master_sequence = master_sequence, complementary_sequence = antigen_sequence, Energy_Matrix = M)
	Sequences = np.array([Master_Sequence])
	sequences = np.array([master_sequence])
	zero_date = datetime(2000, 1, 1)

	file = open('../Text_files/file_parent_daughter2.txt', 'w+')
	file_1 = open('../Text_files/timedNodeFile2.txt', 'w+')
	file_2 = open('../Text_files/strainFile2.txt', 'w+')
	#np.savetxt(file_1, np.array(['Node', 'Parent', 'Time', 'SigClade']), fmt='%d')
	#file.write("\n")
	np.savetxt(file, np.array([str(Master_Sequence.id)+'\t', str(Master_Sequence.parent_id)]), fmt='%s', delimiter='', newline='', header = 'node\t parent\n', comments='')
	file.write("\n")
	np.savetxt(file_1, np.array([str(Master_Sequence.id)+'\t', str(Master_Sequence.parent_id)+'\t', str(0)+'\t' , str(0)]), fmt='%s', delimiter='', newline='', header = 'Node\t Parent\t Time\t SigClade\n', comments='')
	file_1.write("\n")
	n_seq = n_seq
	for i in range(1, n_seq):
	    succ = False 
	    while(succ == False):
	        parent = np.random.choice(Sequences)
	        new_seq = Sequence(seq_id = i, parent = parent.sequence, energy_parent = parent.energy,  parent_id = parent.id, master_sequence = master_sequence, complementary_sequence = antigen_sequence, Energy_Matrix = M)
	        #check if the new sequence is already in the tree. Here we can check for other condicions like that the energy is higher than the parent.
	        if not(np.isin(new_seq.sequence, sequences)):
	            parent.tree_position = 0    
	            Sequences = np.append(Sequences, new_seq)
	            sequences = np.append(sequences, new_seq.sequence)
	            succ = True
	    np.savetxt(file, np.array([str(new_seq.id)+'\t', str(parent.parent_id)]), fmt='%s', delimiter='', newline='')
	    file.write("\n")
	    np.savetxt(file_1, np.array([str(new_seq.id)+'\t', str(parent.parent_id)+'\t', str(new_seq.hamming_distance)+'\t', str(0)]), fmt='%s', delimiter='', newline='')
	    file_1.write("\n")
	    
	for i in range(1, n_seq):
	    if(Sequences[i].tree_position==1):
	        new_date = zero_date + timedelta(Sequences[i].hamming_distance)
	        np.savetxt(file_2, np.array([str(Sequences[i].id)+'\t', 'A/Germany/'+str(i)+'/2020'+'\t', 'Germany'+'\t', 'EPI_ISL_'+str(i)+'\t', str(new_date.year)+'-'+str(new_date.month)+'-'+str(new_date.day)+'\t', '0']),fmt='%s', delimiter='', newline='')
	        file_2.write("\n")
        
	file.close()
	file_1.close()
	file_2.close()

	return Sequences

def plot_histogram_hamming_distance(Sequences, ax):

	distances = np.array([i.hamming_distance for i in Sequences])
	data_distances = np.histogram(distances, bins=range(int(max(distances))))

	#ax.plot(data_distances[1][0:-1], scipy.special.comb(9, data_distances[1][0:-1]), linewidth = 4 , label = 'Binary')
	ax.plot(data_distances[1][0:-1], sc.comb(9, data_distances[1][0:-1])*((20-1)**data_distances[1][0:-1]), color = 'steelblue', linewidth = 4 , label = '20-Alphabet')
	#ax.plot(data_distances[1][0:-1], scipy.special.comb(9, data_distances[1][0:-1])*((4-1)**data_distances[1][0:-1]), linewidth = 4 , label = '4-Alphabet')

	#ax.plot(data_distances[1][0:-1], np.exp(4*data_distances[1][0:-1]), linewidth = 4, label = r'$e^{\lambda r}$')
	ax.plot(data_distances[1][0:-1], data_distances[0], linewidth = 4, label = 'Data', linestyle = 'dashed')

	ax.set_yscale('log')
	#ax.set_ylim(1,1e8)
	ax.set_xlabel(r'Hamming Distance $r$', fontsize = 20)
	ax.set_ylabel(r'', fontsize = 20)
	ax.tick_params(labelsize = 20)
	handles, labels = ax.get_legend_handles_labels()
	ax.legend(np.concatenate(([],handles)),np.concatenate(([],labels)), loc = 2, fontsize = 20)

	return distances

def plot_histogram_energy(Sequences, ax):

	energies = np.array([i.energy for i in Sequences])
	data_energies = np.histogram(energies, bins=50)

	ax.plot(data_energies[1][0:-1], data_energies[0], linewidth = 4, color = 'indianred', label = 'Data', linestyle = 'dashed')
	ax.set_yscale('log')
	#ax.set_ylim(1,1e10)
	ax.set_xlabel(r'Energy $r$', fontsize = 20)
	ax.set_ylabel(r'', fontsize = 20)
	ax.tick_params(labelsize = 20)
	handles, labels = ax.get_legend_handles_labels()
	ax.legend(np.concatenate(([],handles)),np.concatenate(([],labels)), loc = 2, fontsize = 20)

	return energies

def plot_scatter_hamming_distance_energy(distances, energies, ax):

	ax.scatter(distances, energies, color = 'indigo', s = 25, marker = '^')
	ax.hlines(energies[0],0,9, linestyle = 'dashed', label = r'Master Seq. energy $\epsilon_0$')
	#ax.set_yscale('log')
	#ax.set_ylim(1,1e10)
	ax.set_xlabel(r'Hamming distance $d$', fontsize = 20)
	ax.set_ylabel(r'Energy $\epsilon$', fontsize = 20)
	ax.tick_params(labelsize = 20)
	ax.legend(loc = 0, fontsize = 20)



