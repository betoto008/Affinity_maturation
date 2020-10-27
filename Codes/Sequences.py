import numpy as np
import matplotlib.pyplot as plt

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
		self.energy = 0
		self.sequence = parent
		self.pos_mut = np.random.randint(9)
		self.tree_position = 1 # 0 = internal ; 1 = external
		self.hamming_distance = hamming_distance('aaaaaaaaa', self.sequence)
		if not(Master_Seq):
			self.create_sequence()


	def create_sequence(self):
		""" This function will create a new mutations and give an energy value to the new sequence according to the Energy_Matrix. """
		
		list_seq = list(self.sequence)
		old_letter = self.sequence[self.pos_mut]
		self.Alphabet.remove(old_letter)
		new_letter = np.random.choice(self.Alphabet)
		list_seq[self.pos_mut] = new_letter
		self.sequence = "".join(list_seq)
		self.Alphabet = ['a', 'b', 'c', 'd']

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
    