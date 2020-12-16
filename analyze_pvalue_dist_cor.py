import csv
from scipy.stats import gaussian_kde
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def analyze_p_dist(pfam_family_id, vector_cor_w_dist_loc, group_similar_aas):
	"""
	"""
	if group_similar_aas == False:
		out_loc = vector_cor_w_dist_loc + pfam_family_id + "-vector-dist_aa.csv"
		out_image = vector_cor_w_dist_loc + pfam_family_id + "-kde_aa.png"
	elif group_similar_aas == True:
		out_loc = vector_cor_w_dist_loc + pfam_family_id + "-vector-dist_groups.csv"
		out_image = vector_cor_w_dist_loc + pfam_family_id + "-kde_groups.png"

	pvalues, dists, vector_dists = [], [], []
	with open(out_loc, "r") as f:
		reader = csv.reader(f)
		next(reader, None)

		for row in reader:
			pvalue = row[4]
			dist = row[6]

			if len(dist) > 0:
				v_s, v_e = int(row[0]), int(row[2])
				pvalue = float(row[4])
				dist = float(row[6])

				pvalues.append(pvalue)
				dists.append(dist)
				vector_dists.append(v_e - v_s)
	f.close()

	xy = np.vstack([pvalues, dists, vector_dists])
	kde = gaussian_kde(xy)

	num_points = 100
	pvalue_plot = np.linspace(0, max(pvalues), num_points)
	dist_plot = np.linspace(0, max(dists), num_points)
	vector_dist_plot = np.linspace(0, max(vector_dists), num_points)

	# find max
	kde_max = 0.0
	for pvalue in pvalue_plot:
		for dist in dist_plot:
			for vector_dist in vector_dist_plot:
				xy = np.vstack([pvalue, dist, vector_dist])
				kde_value = kde(xy)[0]
				if kde_value > kde_max:
					kde_max = kde_value

	# cut off values below kde threshold
	kde_fract_min = 0.25
	plot_x, plot_y, plot_z, plot_c = [], [], [], []
	for pvalue in pvalue_plot:
		for dist in dist_plot:
			for vector_dist in vector_dist_plot:
				xy = np.vstack([pvalue, dist, vector_dist])
				kde_value = kde(xy)[0]				

				if kde_value >= kde_max * kde_fract_min:
					plot_x.append(pvalue)
					plot_y.append(dist)
					plot_z.append(vector_dist)
					plot_c.append(kde_value)

	fig = plt.figure()
	ax = Axes3D(fig)

	#for i in range(len(plot_x)):
	#	x, y, z, k = plot_x[i], plot_y[i], plot_z[i], plot_c[i] 
	#	ax.scatter(x, y, z, marker = 's', c = k, alpha = k / max(plot_c), cmap = 'Greys', s = 10)
	ax.scatter(plot_x, plot_y, plot_z, marker = 's', alpha = 0.5, c = plot_c, cmap = 'Greys', s = 10)
	#plt.colorbar()
	ax.set_xlabel("P-Value")
	ax.set_ylabel("Distance (Angstroms)")
	ax.set_zlabel("Vector Distance (aas)")
	ax.set_title("Gaussian KDE Fit to " + pfam_family_id + " P-Value/Distance Data" + \
							"\n" + "KDE Values Above " + str(kde_fract_min) + " of KDE Max" + "\n" + \
							"Amino Acids Grouped? " + str(group_similar_aas))

	plt.savefig(out_image, bbox_inches='tight')
	plt.close()

	return