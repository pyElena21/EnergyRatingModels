'''
A collection of examples to test the modules inside the package along with import/read/write 
to CSV file functions

'''


import pvmodels as pvm
import pvchar as pvc


# Some coefficients to test the modules.

# testing pvsat model
def testpvsat():     
    ermodule = {'a1':0.12136, 'a2':-1.187*(10**(-5)), 'a3':0.00240, 'a':-0.0010254}
    power = pvm.pvsat(ermodule, [506.2,598.0,588.7],[15.2,18.7,19.5],1.677)
    print(power)     

#testing the characterisation matrix of pvsat model
def testpvsatMat(geff,tmod,power,area):
    pvsatcoeff = pvc.pvsatMat(geff,tmod,power,area)
    return pvsatcoeff

# testing the characterisation matrix for esti model
def testestierMat(Pstc,geff,tmod,power):
    esticoeff = pvc.estierMat(Pstc,geff,tmod,power)
    return esticoeff

#testing esti pv model
def testestier():
    
    ermodule = {'k1': 0.262947978  ,'k2':  0.166107167  ,'k3':  -0.010676762   ,'k4': -0.004493653  ,'k5':  -1.148000642  ,'k6':  1.03E-05
}
    esticoeff = pvm.estier(245.0, ermodule, [12.7574542,708.2537941,61.79614066],[39.06155432,35.9935333,17.5441417])
    return esticoeff










