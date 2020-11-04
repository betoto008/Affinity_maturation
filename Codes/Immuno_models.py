import numpy as np
import matplotlib.pyplot as plt
from Bio import Phylo
from io import StringIO
from matplotlib.lines import Line2D
from datetime import datetime, timedelta
import scipy.special as sc

class Sequence():
	"""docstring for Sequence"""
	def __init__(self, seq_id, parent, parent_id = 0, Master_Seq = False):
		super(Sequence, self).__init__()

		#A given Sequence object has the following attributes
		self.id = seq_id
		self.Alphabet = ['a', 'b', 'c', 'd']
		self.parent = parent
		self.parent_id = parent_id
		self.active = False
		self.clone_size = 1
		self.energy = 1
		self.sequence = parent
		self.pos_mut = 0
		self.tree_position = 1 # 0 = internal ; 1 = external
		self.hamming_distance = 0

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
			self.probabilities[(2*i)-1] = (self.antigen_time_series[-1]/(self.antigen_time_series[-1]+np.exp(self.master_Sequence_energy + self.Sequences[i-1].hamming_distance)))*(1-self.Sequences[i-1].active)
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

		for i in range(self.T):
			self.gillespie_step()

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
			ax.plot(np.array(range(9)), (1/(1+np.exp(self.master_Sequence_energy + np.array(range(9)) - np.log(rho)))), linewidth  = 4, color = colors[i], label = r'$\rho = %.0e$'%(rho))
		ax.set_yscale('log')
		ax.set_xlabel(r'Hamming Distance $d$', fontsize = 20)
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

	def hist_sequences(self, Sequences, ax):

		rho_array = np.logspace(0, np.log10(max(self.antigen_time_series)), 5)
		colors = plt.cm.Reds(np.linspace(0,1,len(rho_array)))
		data_distances = ax.hist([Sequences[i].hamming_distance for i in range(int(len(Sequences)))], bins = range(10), align = 'left', label = r'$S(d)$', color = 'lightsteelblue', alpha = 0.5)
		ax.plot(data_distances[1][0:-1], sc.comb(9, data_distances[1][0:-1])*((4-1)**data_distances[1][0:-1]), linewidth = 4 , color = 'lightsteelblue', alpha = 0.6)

		ax.hist([self.Sequences[i].hamming_distance for i in range(int(len(self.Sequences)))], bins = range(10), align = 'left', label = r'$US(d)$', color = 'indigo', alpha = 0.6)
		ax.plot(data_distances[1][0:-1], self.U*sc.comb(9, data_distances[1][0:-1])*((4-1)**data_distances[1][0:-1]), linewidth = 4 , color = 'indigo', alpha = 0.6)

		ax.hist([self.Sequences[i].hamming_distance for i in range(int(len(self.Sequences))) if self.Sequences[i].active], bins = range(10), align = 'left', label = r'Activated Linages', color = 'tab:red', alpha = 0.8)
		for i, rho in enumerate(rho_array):
			ax.plot(data_distances[1][0:-1], self.U*sc.comb(9, data_distances[1][0:-1].astype(int))*((4-1)**data_distances[1][0:-1])*(1/(1+np.exp(self.master_Sequence_energy + data_distances[1][0:-1] - np.log(rho)))) , color = colors[i], linestyle = 'dashed', linewidth = 3)

		ax.set_ylim(0.1, 2e5)    
		ax.set_yscale('log')
		ax.set_xlabel(r'Hamming Distance $d$', fontsize = 20)
		ax.set_ylabel(r'Number of linages', fontsize = 20)
		ax.tick_params(labelsize = 20)
		ax.legend(loc = 0, fontsize = 20)
		
def mutate_sequence(Sequence):
		""" This function will create a new mutations and give an energy value to the new sequence according to the Energy_Matrix. """
		Alphabet = ['a', 'b', 'c', 'd']
		Sequence.pos_mut = np.random.randint(9)
		list_seq = list(Sequence.sequence)
		old_letter = Sequence.sequence[Sequence.pos_mut]
		Alphabet.remove(old_letter)
		new_letter = np.random.choice(Alphabet)
		list_seq[Sequence.pos_mut] = new_letter
		Sequence.sequence = "".join(list_seq)
		Sequence.hamming_distance = hamming_distance('aaaaaaaaa', Sequence.sequence)
		Alphabet = ['a', 'b', 'c', 'd']

		#Ask Michael about the best way to produce the energy

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
    
def generate_Sequences(n_seq):

	Master_Sequence = Sequence(seq_id = 0, parent = 'aaaaaaaaa', Master_Seq=True)
	Sequences = np.array([Master_Sequence])
	sequences = np.array(['aaaaaaaaa'])
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
	        new_seq = Sequence(seq_id = i, parent = parent.sequence, parent_id = parent.id)
	        mutate_sequence(new_seq)
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


