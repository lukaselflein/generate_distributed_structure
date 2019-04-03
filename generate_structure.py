""" Generate atomic structure following a given distribution. 
Copyright 2019 Simulation Lab
University of Freiburg
Author: Lukas Elflein <elfleinl@cs.uni-freiburg.de>
"""

import numpy as np
import matplotlib.pyplot as plt

def exponential(x, rate=1):
	return rate * np.exp(-1 * rate * x)


def get_histogram(struc):
	x = struc[:,0]
	y = struc[:,1]
	z = struc[:,2]

	histograms = []	
	for dim_pos in (x, y, z):
		bins = 10
		hist = np.histogram(dim_pos, bins=bins)
		histograms += [hist]

	return histograms


def plot_dist(histogram, name):
	hist, bins = histogram
	width = 0.7 * (bins[1] - bins[0])
	center = (bins[:-1] + bins[1:]) / 2

	# normalize
	hist = hist.astype(float)
	print(hist)
	hist /= hist.sum()
	print(hist)
	
	exp = exponential(center)
	
	fi, ax = plt.subplots()
	ax.bar(center, hist, align='center', width=width)
	ax.plot(center, exp)
	plt.title(name)
	plt.savefig(name + '.png')
	#plt.show()


def invert_histogram(histogram):
	pass


def distance_distributions(dist_1, dist_2):
	pass


def generate_structure(distribution, box=np.array([10, 10, 10])):
	atom_distance = 1
	atom_count = 100
	atom_positions = []
	for i in range(atom_count + 1):
		x = np.random.random_integers(0, box[0])
		y = np.random.random_integers(0, box[1])
		z = np.random.random_integers(0, box[2])

		atom_positions += [[x, y, z]]

	atom_positions = np.array(atom_positions)

	return atom_positions


def main():
	struc = generate_structure(exponential)
	histx, histy, histz = get_histogram(struc)
	plot_dist(histz, 'z')


if __name__ == '__main__':
	main()
