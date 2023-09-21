# Function to save figure in a new folder as svg file
# INPUT: figure and name of new folder
# OUTPUT: figure saved in new folder in current directory

import os
import matplotlib.pyplot as plt

def svg(figure, foldername, figname):
    currentdir = os.getcwd()
    path = os.path.join(currentdir, str(foldername))
    if not os.path.isdir(path):
        os.makedirs(path)
        
    figure.savefig(os.path.join(path, figname + '.svg'))
        
def png(figure, foldername, figname):
    currentdir = os.getcwd()
    path = os.path.join(currentdir, str(foldername))
    if not os.path.isdir(path):
        os.makedirs(path)
        
    figure.savefig(os.path.join(path, figname + '.png'))

def jpeg(figure, foldername, figname):
    currentdir = os.getcwd()
    path = os.path.join(currentdir, str(foldername))
    if not os.path.isdir(path):
        os.makedirs(path)
        
    figure.savefig(os.path.join(path, figname + '.jpeg'))
    
    