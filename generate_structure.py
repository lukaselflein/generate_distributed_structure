""" Generate atomic structure following a given distribution. 
Copyright 2019 Simulation Lab
University of Freiburg
Author: Lukas Elflein <elfleinl@cs.uni-freiburg.de>
"""

import numpy as np
import matplotlib.pyplot as plt

def exponential(x, rate=0.2, cutoff=10):
	return rate * np.exp(-1 * rate * x)


def test():
    import doctest
    doctest.testmod()


def get_histogram(struc, n_bins=20, is_normalized=True):
	x = struc[:,0]
	y = struc[:,1]
	z = struc[:,2]

	histograms = []	
	for dim_pos in (x, y, z):
		n_bins = n_bins
		hist, bins = np.histogram(dim_pos, bins=n_bins, density=is_normalized)

		histograms += [(hist, bins)]
	return histograms


def get_centers(bins):
	"""Return the center of the provided bins.

	Example:
	>>> get_centers(bins=np.array([0.0, 1.0, 2.0]))
	array([0.5, 1.5])
	"""	
	return (bins[:-1] + bins[1:]) / 2


def get_nearest_pos(array, value):
	"""Find the value of an array clostest to the second argument.

	Example:
	>>> get_nearest_pos(array=[0, 0.5, 1.0], value=0.6)
	1
	"""	
	array = np.asarray(array)
	pos = np.abs(array - value).argmin()
	return pos


def inverted_distribution(distribution, p, support=None):
	"""Inverts a distribution x->p, and returns the x-value belongig to the provided p.

	Arguments:
	distribution: a function x -> p; x should be approximatable by a compact support
	p: an output of the distribution function, probablities in (0,1) are preferrable
	
	>>> inverted_distribution(lambda x: x**2, 16)
	4.0
	"""
	if support is None:
		# Define the x-values to evaluate the function on
		support = np.arange(0,10,0.1)

	# Calculate the histogram of the distribution
	hist = distribution(support)

	# If the p is not in the image of the support, get the nearest hit instead
	nearest_pos = get_nearest_pos(hist, p)

	# Get the x-value belonging to the probablity value provided in the input
	x = support[nearest_pos]
	return x


def plot_dist(histogram, name):
	hist, bins = histogram
	width = 0.7 * (bins[1] - bins[0])
	centers = get_centers(bins)

	# normalize
	hist = hist.astype(float)
	
	exp = exponential(centers)
	
	fi, ax = plt.subplots()
	ax.bar(centers, hist, align='center', width=width)
	ax.plot(centers, exp)
	plt.title(name)
	plt.xlabel('Distance ' + name)
	plt.savefig(name + '.png')


def generate_structure(distribution, box=np.array([10, 10, 10])):
	atom_count = 100
	atom_positions = []
	for i in range(atom_count + 1):
		# z is distributed according to the given distribution
		# To approximate this, we insert an atom with probablity dis(z) at place z.
		# This we do by inverting the distribution, and sampling uniformely from distri^-1:
		p = np.random.uniform()
		z = inverted_distribution(distribution, p)

		# x and y are uniformely distributed
		x = np.random.random_integers(0, box[0])
		y = np.random.random_integers(0, box[1])
		x = 1
		y = 2

		atom_positions += [[x, y, z]]
	print(atom_positions)
	atom_positions = np.array(atom_positions)

	return atom_positions
	

def main():
	struc = generate_structure(exponential)
	histx, histy, histz = get_histogram(struc)
	print(histx)
	
	plot_dist(histz, 'z')
	# plot_dist(histx, 'x')
	# plot_dist(histy, 'y')

	

	if False:
		centers = np.arange(0,10,0.1)
		fig, ax = plt.subplots()
		ax.plot(centers, exponential(centers))
		plt.show()
	 
		inverted_histz = invert_histogram(histz)

		centers = get_centers(histz[1])
		exp = exponential(centers)
		centers, exponential(centers)

		fig, ax = plt.subplots()
		ax.plot(exponential(centers), centers)
		plt.show()


if __name__ == '__main__':
	# Run doctests
	test()
	# Execute everything else
	main()
