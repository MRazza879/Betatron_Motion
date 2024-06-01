#!/usr/bin/env python3
#Author: Marta Razza
#This script reads the files containing the values of the dynamical indicator: Lyapunov error(LE) and Reversibility error method (REM) and plots them.

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys

# Function to read the file
def reading_file(filename):
    try:
        data = pd.read_csv(filename, header=None, sep='\t')
        return data
    except FileNotFoundError:
        print(f"Error: File '{filename}' does not exist.")
        sys.exit(1)
    except pd.errors.EmptyDataError:
        print(f"Error: File '{filename}' is empty.")
        sys.exit(1)
    except Exception as e:
        print(f"Error: An error occurred while reading '{filename}': {e}")
        sys.exit(1)

def main():
    # Verifying the number of arguments
    if len(sys.argv) < 3:
        print("Please provide the name of 2 different files to be read. The order is: lyapunov and reversibility.")
        return

    lyapunov_file = sys.argv[1]
    reversibility_file = sys.argv[2]

    process_indicators(lyapunov_file, "Lyapunov_error_4D_Hénon_map")
    process_indicators(reversibility_file, "Reversibility_error_method_4D_Hénon_map")

# Class to plot the dynamical indicators
class Indicators:
    def __init__(self, filename, title):
        self.filename = filename
        self.data = reading_file(filename)
        self.data.columns = ['x0', 'p0', 'indicator']
        self.title = title

    def plot(self):
        x = np.array(self.data['x0'])
        p = np.array(self.data['p0'])
        indicator = np.array(self.data['indicator'])

        # Creation of the figure with higher resolution
        fig, ax = plt.subplots(figsize=(8, 6), dpi=300)

        # Definition of the color map
        cmap = plt.colormaps['viridis']

        # Creation of the color map with smaller markers
        sc = ax.scatter(x, p, c=indicator, cmap=cmap, marker='s', s=1)

        # Adding the color bar
        cbar = plt.colorbar(sc)

        # Graph settings
        ax.set_title(self.title)
        ax.set_xlabel("x")
        ax.set_ylabel("p")
        #if you want to change the dimension of lattice remember to change these values below
        ax.set_xlim([0, 0.45])
        ax.set_ylim([0, 0.45])

    def save_plot(self, filename):
        plt.savefig(filename, dpi=300)
        
def process_indicators(filename, title):
    indicators = Indicators(filename, title)
    indicators.plot()
    indicators.save_plot(f'{title}.png')
    

if __name__ == "__main__":
    main()


