#!/usr/bin/env python3
# Script to extract 2D plots for tracking, Lyapunov error and reversibility error method with data obtained with 2DHenonMapc.c
# Author: Marta Razza

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
from matplotlib.colors import LogNorm

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
        print("Please provide the name of 2 different files to be read. The order is: tracking and lyapunov.")
        return

    tracking_file = sys.argv[1]
    lyapunov_file = sys.argv[2]

    process_tracking(tracking_file)
    process_indicators(lyapunov_file, "Lyapunov_error")

class Tracking:
    def __init__(self, filename):
        self.filename = filename
        self.data = reading_file(filename)
        self.data.columns = ['x', 'y']

    def plot(self):
        plt.figure(figsize=(10, 10), dpi=300)
        plt.scatter(self.data['x'], self.data['y'], c='blue', marker='.', s=1)
        plt.xlabel('x')
        plt.ylabel('p')
        plt.xlim([-1, 1])  # Set x-axis limits
        plt.ylim([-1, 1])  # Set y-axis limits
        plt.title('Phase-Space Portrait of the 2D Henon Map')

    def save_plot(self, filename):
        plt.savefig(filename)

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
        cmap = plt.colormaps['jet']
        norm = LogNorm(vmin=1, vmax=1E3)

        # Creation of the color map with smaller markers
        sc = ax.scatter(x, p, c=indicator, cmap=cmap, norm=norm, marker='s', s=1)

        # Adding the color bar
        cbar = plt.colorbar(sc)
        cbar.set_label('log10(indicator)')

        # Graph settings
        ax.set_title(self.title)
        ax.set_xlabel("x")
        ax.set_ylabel("p")
        ax.set_xlim([-1, 1])
        ax.set_ylim([-1, 1])

    def save_plot(self, filename):
        plt.savefig(filename, dpi=300)

def process_tracking(filename):
    tracking = Tracking(filename)
    tracking.plot()
    tracking.save_plot('tracking.png')

def process_indicators(filename, title):
    indicators = Indicators(filename, title)
    indicators.plot()
    indicators.save_plot(f'{title}.png')

if __name__ == "__main__":
    main()
