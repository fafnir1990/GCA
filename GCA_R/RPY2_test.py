import rpy2
import rpy2.robjects as robjects
import rpy2.interactive.packages

from rpy2.robjects.packages import importr
robjects.r["set.seed"](1)
copula = importr("copula")
sp3D = importr("scatterplot3d")

myCop = copula.archmCopula(family="clayton", dim =7, param=2)

#mydata = read.table("~/projekt_indywidualny/GCA2/input")

# myMvd = copula.mvdc(copula=myCop, margins=robjects.StrVector(["norm", "norm", "norm","norm", "norm", "norm","norm"]),
#               paramMargins = robjects.vectors.ListVector([
#                   robjects.vectors.ListVector([mean=0.66,sd = 0.1]), robjects.vectors.ListVector(mean=0.66,sd = 0.1),
#                   robjects.vectors.ListVector(mean=0.66,sd = 0.1), robjects.vectors.ListVector(mean=0.66,sd = 0.1),
#                   robjects.vectors.ListVector(mean=0.66,sd = 0.1), robjects.vectors.ListVector(mean=0.66,sd = 0.1),
#                   robjects.vectors.ListVector(mean=0.66,sd = 0.1) ]))
myMvd = copula.mvdc(copula=myCop, margins=robjects.StrVector(["norm", "norm", "norm","norm", "norm", "norm","norm"]),
              paramMargins = robjects.vectors.ListVector([

                  robjects.vectors.FloatVector([0.66, 0.1]), robjects.vectors.FloatVector([0.66, 0.1]),
                  robjects.vectors.FloatVector([0.66, 0.1]), robjects.vectors.FloatVector([0.66, 0.1]),
                  robjects.vectors.FloatVector([0.66, 0.1]), robjects.vectors.FloatVector([0.66, 0.1]),
                  robjects.vectors.FloatVector([0.66, 0.1]) ]))

x = copula.rmvdc(myMvd, 1000)

#contour(myMvd, dmvdc, xlim = c(0, 1), ylim = c(1, 4))

robjects.r['jpeg']('rplot.jpg')
copula.contour(x)
robjects.r['dev.off']()

'''
robjects.r['jpeg']('rplot2.jpg')
sp3D.scatterplot3d(x)
robjects.r['dev.off']()
'''
#robjects.r['save'](args)
