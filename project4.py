from pylab import *
import os, shutil, sys

os.chdir("build-project4-Desktop-Debug")

def plotting(subdir, title, anititle):
    infile = open(title,'r')

    if os.path.isdir(subdir):
        shutil.rmtree(subdir)
    os.mkdir(subdir)
    os.chdir(subdir)

    counter = 0

    for line in infile:
        splitline = line.split()
        x = zeros(len(splitline))
        for i in range(len(splitline)):
            x[i] = float(splitline[i])
        plot(x)
        savefig("tmp%04d.png"%counter)
        counter += 1
    
        #os.system("mencoder 'mf://_tmp*.png' -mf type=png:fps=10 -ovc lavc -lavcopts vcodec=wmv2 -oac copy -o animation.mpg")

    os.chdir(os.pardir) # Move to build-dir

plotting("tempFE", "ForwardEuler.txt", "FEuler.mpg")



