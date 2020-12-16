from config import pkl_plot_folder
import matplotlib.pyplot as plt
import os
import pickle

print(pickle.format_version)

for filename in os.listdir(pkl_plot_folder):
	if "pickle" in filename:
		figx = pickle.load(open(pkl_plot_folder + filename, "rb"))

		#figx.show()
		plt.show()
