from pylab import *
import os, shutil, sys

os.chdir("build-project4-Desktop-Debug")

def plotting(subdir, Title, anititle):

    infile = open(Title,'r')

    if os.path.isdir(subdir):
        shutil.rmtree(subdir)
    os.mkdir(subdir)
    os.chdir(subdir)

    counter = 0
    
    xt1 = []; xt2 = []
    for line in infile:
        splitline = line.split()
        u = zeros(len(splitline))
        for i in range(len(splitline)):
            u[i] = float(splitline[i])
            if counter == 5:
                xt1.append(u[i])
            elif counter == 30:
                xt2.append(u[i])

        x = linspace(0,1,len(u))

        plot(x,u,'r',hold=False)
        title(Title)
        xlim(0,1); ylim(0,1)
        xlabel("x")
        ylabel("u")
        savefig("tmp%04d.png"%counter)
        counter += 1
        
    os.system("mencoder 'mf://tmp*.png' -mf type=png:fps=10 -ovc lavc -lavcopts vcodec=wmv2 -oac copy -o %s"%anititle)
    infile.close()
    os.chdir(os.pardir) # Move to build-dir
    return xt1, xt2, x

xt11, xt21, x = plotting("tempFE", "ForwardEuler.txt" , "FEuler.mpg")
xt12, xt22, x = plotting("tempBE", "BackwardEuler.txt", "BEuler.mpg")
xt13, xt23, x = plotting("tempCN", "CrankNicolson.txt", "CrankN.mpg")
xt14, xt24, x = plotting("tempAn", "Analytical.txt"   , "Analytic.mpg")

plot(x, xt11, label="Forward Euler",  hold=False)
plot(x, xt12, label="Backward Euler", hold=True)
plot(x, xt13, label="Crank Nicolson", hold=True)
plot(x, xt14, label="Analytical",     hold=True)
title("Comparing solvers at t1. dx = 0.1 and dt = 0.01")
legend(loc="upper right")
xlim(0,1); ylim(0,1)
xlabel("x")
ylabel("u")
savefig("solverst1.png")
show()

plot(x, xt21, label="Forward Euler",  hold=False)
plot(x, xt22, label="Backward Euler", hold=True)
plot(x, xt23, label="Crank Nicolson", hold=True)
plot(x, xt24, label="Analytical",     hold=True)
title("Comparing solvers at t1")
xlim(0,1); ylim(0,1)
legend(loc="upper right")
xlabel("x")
ylabel("u")
savefig("solverst2.png")
show()
