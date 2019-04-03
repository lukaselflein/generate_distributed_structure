""" Generate atomic structure following a given distribution. 
Copyright 2019 Simulation Lab
University of Freiburg
Author: Lukas Elflein <elfleinl@cs.uni-freiburg.de>
"""

import numpy as np
import matplotlib.pyplot as plt

def exponential(x, rate=0.3, cutoff=10):
	"""Exponential distribution."""
	return rate * np.exp(-1 * rate * x)


def uniform(x, *args, **kwargs):
	"""Uniform distribution."""
	return np.ones(x.shape)


def test():
	"""Run docstring unittests"""
	import doctest
	doctest.testmod()


def pdf_to_cdf(pdf):
	"""Transform partial distribution to cumulative distribution function

	>>> pdf_to_cdf(pdf=np.array([0.5, 0.3, 0.2]))
	array([0.5, 0.8, 1. ])
	"""
	cdf = np.cumsum(pdf)
	cdf /= cdf[-1]

	return cdf


def get_centers(bins):
	"""Return the center of the provided bins.

	Example:
	>>> get_centers(bins=np.array([0.0, 1.0, 2.0]))
	array([0.5, 1.5])
	"""	
	bins = bins.astype(float)
	return (bins[:-1] + bins[1:]) / 2


def get_nearest_pos(array, value):
	"""Find the value of an array clostest to the second argument.

	Example:
	>>> get_nearest_pos(array=[0, 0.25, 0.5, 0.75, 1.0], value=0.55)
	2
	"""	
	array = np.asarray(array)
	pos = np.abs(array - value).argmin()
	return pos


def get_histogram(struc, bins=np.arange(0, 11)):
	"""Slice the list of atomic positions, aggregate positions into histogram."""
	# Extract x/y/z positions only
	x, y, z = struc[:, 0], struc[:, 1], struc[:, 2]

	histograms = []	
	for dimension_position in (x, y, z):
		hist, bins = np.histogram(dimension_position, bins=bins, density=True)
		# Normalize the histogram for all values to sum to 1
		hist /= sum(hist)

		histograms += [(hist, bins)]
	return histograms


def plot_dist(histogram, name, reference_distribution=None):
	"""Plot histogram with an optional reference distribution."""
	hist, bins = histogram
	width = 0.95 * (bins[1] - bins[0])
	centers = get_centers(bins)

	fi, ax = plt.subplots()
	ax.bar(centers, hist, align='center', width=width, label='Empirical distribution')

	if reference_distribution is not None:
		ref = reference_distribution(centers)
		ref /= sum(ref)
		ax.plot(centers, ref, color='red', marker='o', label='Reference distribution')

	plt.title(name)
	plt.legend()
	plt.xlabel('Distance ' + name)
	plt.savefig(name + '.png')


def quartile_function(distribution, p, support=None):
	"""Inverts a distribution x->p, and returns the x-value belonging to the provided p.

	Assumption: The distribution to be inverted must have a strictly increasing CDF!
	Also see 'https://en.wikipedia.org/wiki/Quantile_function'.

	Arguments:
	distribution: a function x -> p; x should be approximatable by a compact support
	p: an output of the distribution function, probablities in (0,1) are preferrable
	"""
	if support is None:
		# Define the x-values to evaluate the function on
		support = np.arange(0,1,0.01)

	# Calculate the histogram of the distribution
	hist = distribution(support)

	# Sum the distribution to get the cumulatative distribution
	cdf =  pdf_to_cdf(hist)

	# If the p is not in the image of the support, get the nearest hit instead
	nearest_pos = get_nearest_pos(cdf, p)

	# Get the x-value belonging to the probablity value provided in the input
	x = support[nearest_pos]
	return x


def quartile_sampler(distribution, support):
	"""Wrapper for quartile_function."""
	# z is distributed according to the given distribution
	# To approximate this, we insert an atom with probablity dis(z) at place z.
	# This we do by inverting the distribution, and sampling uniformely from distri^-1:
	p = np.random.uniform()
	sample = quartile_function(distribution, p, support=support)
	
	return sample


def rejection_sampler(distribution, support, max_tries=1000):
	"""Sample distribution by drawing from support and keeping according to distribution.
	
	Draw a random sample from our support, and keep it if another random number is
	smaller than our target distribution at the support location.

	Arguments
	distribution: The target distribution, as a histogram over the support
	support: locations in space where our distribution is defined
	max_tries: how often the sampler should attempt to draw before giving up. 
		   If the distribution is very sparse, increase this parameter to still get results.

	Returns
	sample: a location which is conistent (in expectation) with being drawn from the distribution.
	"""

	for i in range(max_tries):
		sample = np.random.choice(support)
		
		# Keep sample with probablity of distribution
		if np.random.random() < distribution(sample):
			return sample

	raise RuntimeError('Maximum of attempts max_tries {} exceeded!'.format(max_tries))
	

def generate_structure(distribution, box=np.array([10, 10, 10])):
	"""Generate an atomic structure.
	
	Z coordinates are distributed according to the given distribution.
	To construct the positions, we insert an atom with probality distribution(z) at place z.
	This sampling is done by inverting the distribution, and sampling uniformely from distri^-1.
	
	X and Y coordinates are drawn uniformely.

	Arguments:

	"""
	atom_count = 50000
	atom_positions = []

	# We define which positions in space the atoms can be placed
	support = {}
	# Using the box parameter, we construct a grid inside the box
	# With gridpoints every 0.1 units:
	grid_density = 1
	# This results in a 100x100x100 grid:
	support['x'] = np.arange(0, box[0], grid_density)
	support['y'] = np.arange(0, box[1], grid_density)
	support['z'] = np.arange(0, box[2], grid_density)
	
	# For every atom, draw random x, y and z coordinates
	for i in range(atom_count + 1):
		# Z coordinate is distributed non-uniformely
		z = rejection_sampler(distribution, support['z'])

		# x and y are uniformely distributed
		x = rejection_sampler(uniform, support['x'])
		y = rejection_sampler(uniform, support['y'])

		atom_positions += [[x, y, z]]
	atom_positions = np.array(atom_positions)

	return atom_positions
	

def main():
	"""Generate an example structure, plot the distributions in this example."""
	struc = generate_structure(exponential)
	histx, histy, histz = get_histogram(struc)
	h = get_histogram(struc)
	
	plot_dist(histz, 'z', reference_distribution=exponential)
	plot_dist(histx, 'x', reference_distribution=uniform)
	plot_dist(histy, 'y', reference_distribution=uniform)


if __name__ == '__main__':
	# Run doctests
	test()
	# Execute everything else
	main()
