__author__ = 'rav'
import numpy as np
from numpy import genfromtxt
from rpy2.robjects.vectors import ListVector, FloatVector
from rpy2 import robjects
from rpy2.robjects import r
from rpy2.robjects.packages import importr
import rpy2.robjects.numpy2ri



rpy2.robjects.numpy2ri.activate()

# Tutaj mozna wczytac plik csv
my_data = genfromtxt('data.csv', delimiter=',')
#print type(my_data)
my_data = np.asarray(map(list, zip(*my_data)))
#print type(my_data)




data = my_data
# Zmiana nazwy mapy ciepla wejsciowych danych
robjects.r['jpeg']('rplot3.jpg')
r.heatmap(data)
robjects.r['dev.off']()
print "DATA"
print my_data
mean_sd = []

col1 = data[:,0]
col2 = data[:,1]
print 'BBB'
print my_data[1]

#obliczenie sredniej i odchylenia stanadardowego
for i in range(len(my_data)):

    mean_sd.append([np.mean(my_data[i]), np.std(my_data[i])])
n = len(mean_sd)

my_data1 = my_data

#data = robjects.r['c'](data)
my_data = robjects.r['t'](data)
print my_data[0], my_data[1], my_data[2], my_data[3]
print my_data


'''
for column in map(list, zip(my_data)):
    mean_sd.append([np.mean(column), np.std(column)])
'''


print "meansd"
print mean_sd

list_dict = [{'mean' : i[0], 'sd' : i[1]} for i in mean_sd]
list_Margins = [("V" + str(i + 1), ListVector(list_dict[i])) for i in range(n)]
#print list_Margins



robjects.r["set.seed"](1)

#import pakietow bezposrednio z R (biblioteki doinstalowane)
copula = importr("copula")
sp3D = importr("scatterplot3d")

#family pozwala na zmiane kopuli ["clayton", "gumbel", "frank"]
myCop = copula.archmCopula(family="gumbel", dim=n, param=2)


#Dane pochodzace z dwoch pierwszych kolumn tabeli
families = ["clayton", "gumbel", "frank"]
for f in families:
    myCop2 = copula.archmCopula(family=f,  param=2, dim=3)
    # paramMargins2 = ListVector([   ("V1", ListVector([{'mean': mean_sd[0][0], 'sd': mean_sd[0][1]}])),
    #                                ("V2", ListVector([{'mean': mean_sd[1][0], 'sd': mean_sd[1][1]}]))])
    paramMargins2 = ListVector([("V1", ListVector({'mean' : mean_sd[0][0], 'sd' : mean_sd[0][1]})),
                                ("V2", ListVector({'mean' : mean_sd[1][0], 'sd' : mean_sd[1][1]})),
                                ("V3", ListVector({'mean' : mean_sd[2][0], 'sd' : mean_sd[2][1]}))])
    margins2 = robjects.StrVector(["norm", "norm", "norm"])
    myMvd2 = copula.mvdc(copula=myCop2, margins=margins2, paramMargins=paramMargins2)
    X2 = copula.rmvdc(myMvd2, 4)
    robjects.r['jpeg'](f + '.jpg')
    xlim = r['c'](-10,30)
    ylim = r['c'](-10, 30)
    #zlim = r['c'](-10, 30)
    #copula.contour(myMvd2, copula.dmvdc,  xlim=xlim, ylim=ylim, zlim=zlim)
    sp3D.scatterplot3d(X2)
    robjects.r['dev.off']()



#paramMargins = ListVector(map(FloatVector, mean_sd))
paramMargins = ListVector(list_Margins)



#ListVector([("V1", ListVector([{"mean": 0.66, "sd": 0.1}])),
#            ("V2", ListVector([{"mean": 0.66, "sd": 0.1}]))])





margins = robjects.StrVector(n * ["norm"])

myMvd = copula.mvdc(copula=myCop, margins=margins, paramMargins=paramMargins)
print myMvd
# myMvd = copula.mvdc(copula=myCop, margins=robjects.StrVector(["norm", "norm", "norm","norm", "norm", "norm","norm"]),
#               paramMargins = robjects.vectors.ListVector([
#                   robjects.vectors.ListVector([mean=0.66,sd = 0.1]), robjects.vectors.ListVector(mean=0.66,sd = 0.1),
#                   robjects.vectors.ListVector(mean=0.66,sd = 0.1), robjects.vectors.ListVector(mean=0.66,sd = 0.1),
#                   robjects.vectors.ListVector(mean=0.66,sd = 0.1), robjects.vectors.ListVector(mean=0.66,sd = 0.1),
#                   robjects.vectors.ListVector(mean=0.66,sd = 0.1) ]))



# myMvd = copula.mvdc(copula=myCop, margins=robjects.StrVector(["norm", "norm", "norm","norm", "norm", "norm","norm"]),
#               paramMargins = robjects.vectors.ListVector([
#
#                   robjects.vectors.FloatVector([0.66, 0.1]), robjects.vectors.FloatVector([0.66, 0.1]),
#                   robjects.vectors.FloatVector([0.66, 0.1]), robjects.vectors.FloatVector([0.66, 0.1]),
#                   robjects.vectors.FloatVector([0.66, 0.1]), robjects.vectors.FloatVector([0.66, 0.1]),
#                   robjects.vectors.FloatVector([0.66, 0.1]) ]))

x = copula.rmvdc(myMvd, 4)
print x
#contour(myMvd, dmvdc, xlim = c(0, 1), ylim = c(1, 4))
xlim = r['c'](0, 1)
ylim = r['c'](0, 1)
robjects.r['jpeg']('rplot.jpg')
copula.contour(x,  xlim = xlim, ylim = ylim)
#copula.contour(myMvd, copula.dmvdc,  xlim = xlim, ylim = ylim)
robjects.r['dev.off']()

params = []
map(params.extend, mean_sd)
params.append(2)
#print robjects.vectors.FloatVector(params)
qual = copula.loglikMvdc(robjects.vectors.FloatVector(params), x, myMvd)

abc = copula.fitCopula(myCop, my_data, method="itau")


mm = r['apply'](my_data, 2, r['mean'])
vv  = r['apply'](my_data, 2, r['var'])
b1 = r['c'](mm[1]**2/vv[1], vv[1]/mm[1])
b2 = r['c'](mm[2]**2/vv[2], vv[2]/mm[2])

# Just First Column with All rows: R my_data[, 1], but I use numpy before transopose
print 'KORELACJA'
#print r['cor'](my_data1[0], my_data1[1], method = "kendall")[0] * 2
a = r['sin'](r['cor'](my_data1[0], my_data1[1], method = "kendall")[0] * np.pi/2)
start = r['c'](params)
#start = r['c'](b1, b2, a)

def_ = copula.fitMvdc(my_data, myMvd, start = start)
#print qual
#print abc
print def_

print 'MY MVD'
print myMvd
print 'MY COP'
print myCop
