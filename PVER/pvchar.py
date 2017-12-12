'''
Building the characterisation matrices (dictionaries) for the pv modules
This module develops the fitting procedures for the pvmodels.
The coefficients are automatically entered into a dictionary.

----------------------------------------------------------
 Note: The way to declare a dictionary is the following: 
----------------------------------------------------------

dict = {"name":value, "name2": value2,...."nameN": valueN}

References:
-----------

[1] Friesen, R. Gottschalg, H.G.Beyer, S. Williams, A. Guerin de Montgareuil, N. van der Borg, W.G.J.H.M. van Sark, 
T. Huld B. Muller A.C. de Keizer Y. Niu Intercomparison Of Different Energy Prediction Methods -
Within The European Project Performance - Results Of The 1st Round Robin, 22nd European Photovoltaic Solar Energy Conference,2007

[2] Friesen, G. Dittmann, S. Williams S.Gottschalg R.Beyer H.G.
Guerin de Montgareuil, A.Van Der Borg, N.J.C.M.Burgers, A.R.Kenny, R.P.Huld, T.Muller, B
Reise, C. Kurnik, J. Topic, M., Intercomparison of Different Energy Prediction Methods within the European Project "Performance" - 
Results on the 2nd Round Robin,24th European Photovoltaic Solar Energy Conference, 2009 

[3] Dittmann, S.,Friesen, G., Williams, S., Betts, T.R., Gottschalg, R., Beyer, H.G., Montgareuil, A. Guerin de,
Borg, N.v.d., Burgers, A.R., Huld, T., Muller, B., Reise, C., Kurnik, J., Topic, M., Zdanowicz, T. F. Fabero, 
Results Of The 3rd Modelling Round Robin Within The European Project Performance - Comparison 
Of Module Energy Rating Methods, 25th European Photovoltaic Solar Energy Conference and Exhibition / 5th World Conference on Photovoltaic Energy Conversion, Valencia, Spain, September 2010

'''

try:
    import numpy as np
except ImportError():
    raise ImportError("You need to install numpy")

try:
    from scipy.optimize.minpack import curve_fit
except ImportError():
    raise ImportError("You need to install scipy") 
  



def pvsatMat(Geff,Tmod,power,area, **kwargs):
    
    ''' 
    Origin: H2M
    Coefficients: four 
    Result: the four coefficients
    
    Parameters
    ------------
      
               
    Geff:     1d array or list
              indoor irradiance measurements (Watt/m2)
               
    Tmod:     1d array or list
              indoor module temperature measurements.
              
              
    power:    1d array or list
              power (Watt)
    
              
    Area:     numpy float
              the area of the module
              
              
    Returns:  a dictionary with the four coefficients of the model
    --------
    
    Note: Coefficient symbols are consistent to the source
    
    
    '''
    
    
    Geff = np.array((Geff),dtype = float)
    Tmod = np.array((Tmod),dtype = float)
    x = np.array((Geff,Tmod),dtype = float)
    power = np.array(power)
    
    # quality control
    Geff[Geff<0.0] = 0.0001   
    power[power<0.0] = 0.0
    
    
    def _pvs(x,a,a1,a2,a3):   
        p = (a1 + a2*x[0] +a3*np.log(x[0])) * (1+a*(x[1]-25))* area * x[0]
        return p           
        
    
    params, pcov = curve_fit(_pvs,x,power) # currently works with Lavenberg- Marquardt optimisation algorithm
    ermodule = {"a":params[0], "a1":params[1], "a2":params[2], "a3":params[3]}
    
    print(ermodule)     
    
    return ermodule





def estierMat(Pstc,Geff,Tmod,power, **kwargs):

    
    ''' 
    Origin: JRC (Joint Research Centre)
    Coefficients: six 
    Result: coefficients (fitting data)
    
    Parameters
    ------------
        
               
    Geff:     numpy float or 1darray or list
              indoor irradiance measurements (Watt/m2)
               
    Tmod:     numpy float or 1darray or list
              indoor module temperature.
              
              
    Pstc:     numpy float. 
              The maximum power at standard testing conditions STC.
    
    power:    numpy float or 1d array
              power returned for the given irradiance, module temperature (Watt)
              
              
    Returns:  
    
    ermodule: dict 
             A dictionary containing the coefficients of the module
             It should contain the six fields required for this model 
    --------
    
    Note: Coefficient symbols are consistent to the source
    
    Reference:
    ----------
    Huld, Thomas, Friesen, Gabi, Skoczek, Artur, Kenny, Robert P., Sample, Tony, Field, Michael
    Dunlop, Ewan D.A, power-rating model for crystalline silicon PV modules, Solar Energy Materials 
    and Solar Cells, 2011
    
    
    '''          
     
    
    
    Geff = np.array((Geff),dtype = float)
    Tmod = np.array((Tmod),dtype = float)
    
    x = np.array((Geff, Tmod), dtype = float)
    power = np.array((power), dtype = float)
    
    Geff[Geff<0] = 0.0001
    power[power<0] = 0.0
    
    def _est(x,k1,k2,k3,k4,k5,k6):
        p = Pstc*0.001*x[0]*(1+k1*np.log(0.001*x[0])+k2*((np.log(0.001*x[0]))**2)+
                +(x[1]-25)*(k3+(k4)*(np.log(0.001*x[0])+k5*((np.log(0.001*x[0])**2)))+ k6*((x[1]-25)**2)))       
        return p
      
    
    params, pcov = curve_fit(_est,x,power) # currently works with Lavenberg- Marquardt optimisation algorithm
    ermodule = {"k1":params[0], "k2":params[1], "k3":params[2], "k4":params[3], "k5":params[4], "k6":params[5] }
    
    #power = power.values[0]
    print(ermodule)
    return ermodule




def zenitMat(Geff,Tmod,power, **kwargs):
    
    
    ''' 
    Origin: ISE
    Coefficients: five 
    Result: power (power surface)
    
    Parameters
    ------------
    
               
    Geff:     numpy float or 1darray or list
              indoor irradiance measurements (Watt/m2)
               
    Tmod:     numpy float or 1darray or list
              module temperature indoor measurements in celsius.
              
              
    power:    numpy float or 1d array or list
              power returned for the given irradiance, module temperature (Watt)
              
              
    Returns:  
    
    ermodule: dict 
    A dictionary containing the coefficients of the module
    It should contain the six fields required for this model 
    --------
    
    Note: Coefficient symbols are consistent to the source
    
    '''
    
    Geff = np.array((Geff),dtype = float)
    Tmod = np.array((Tmod),dtype = float)
    x = np.array((Geff,Tmod),dtype = float)
    power = np.array(power)
    
    Geff[Geff<0.0] = 0.0001   
    power[power<0.0] = 0.0

    
    
    def _zen(x,a,b,c,d,e):
    
        psi = ((np.log(x[0]+e))**2)/(x[0]+1)
    
        p = a*(x[0]**2) + b*np.log(x[0]+1) * Geff+c*(psi-1)*x[0] + d*(x[1]-25)
        
        return p
    
    
    params, pcov = curve_fit(_zen,x,power)
    ermodule = {"a":params[0], "b":params[1], "c":params[2], "d":params[3], "e":params[4]}
    print(ermodule)
    
    return ermodule


def matrixmetMat(Geff,Tmod,Tcell,power, **kwargs):
    
    ''' 
    Origin: SUPSI
    Coefficients: six 
    Result: coefficients 
    
    Parameters
    ------------
       
               
    Geff:     numpy float or 1darray or list
              indoor irradiance measurements (Watt/m2)
               
    Tmod:     numpy float or 1darray or list
              module temperature indoor measurements.
              
              
    Tcell:    numpy float or 1darray or list
              cell temperature indoor measurements
    
    power:    numpy float or 1d array or list
              power returned for the given irradiance, module temperature (Watt)
    
    Returns:
    --------
              ermodule: dict 
              A dictionary containing the coefficients of the module
              It should contain the six fields required for this model 
    
    
    '''
    

    Tcell = np.array(Tcell,dtype = float)
    Tmod = np.array(Tmod, dtype = float)
    Geff = np.array((Geff),dtype = float)
    power = np.array((power), dtype = float)
    
    Geff[Geff<0.0] = 0.0001   
    power[power<0.0] = 0.0
    deltat = Tmod - Tcell 
    
    x = np.array((Geff, Tmod, deltat), dtype = float)
    
    
    
    def _mat(x,im_stc,a_im, vm_stc,b_vm,c0,c1):
        
        imax = (im_stc*x[0]/1000) * (1+a_im*(x[2]+x[1]-25))
        vmax = vm_stc+c0*np.log(x[0]/1000) + c1*((np.log(x[0]/1000))**2) + b_vm*(x[2]+x[1]-25)
        p = imax*vmax
        return p
    
    
    params, pcov = curve_fit(_mat,x,power)
    ermodule = {"im_stc":params[0], "a_im":params[1], "vm_stc":params[2], "b_vm":params[3], "c0":params[4], "c1":params[5]}
    print(ermodule)
    
    return ermodule
    


def mat61853_1(Geff,Tmod,power, **kwargs):
    
    '''
    Parameters
    ------------
       
               
    Geff:     numpy float or list
              indoor irradiance measurements (Watt/m2)
               
    Tmod:     numpy float or list
              module temperature indoor measurements.
              
              
     power:    numpy float or 1d array
              power returned for the given irradiance, module temperature (Watt)
    
    
    Returns:
    --------
              ermodule: dict 
              A dictionary containing the coefficients of the module
              It should contain the six fields required for this model 
    
    
    ''' 
    
    Tmod = np.array(Tmod, dtype = float)
    Geff = np.array((Geff),dtype = float)
    power = np.array((power), dtype = float)
    
    Geff[Geff<0.0] = 0.0001   
    power[power<0.0] = 0.0
    
    x = np.array((Geff, Tmod), dtype = float)
    
    
    
    
    
    
    def _matiec(x,a0,a1,a2,a3,a4):
        
        #power = (a0 + a1*Geff + a2*(Geff**2))*(a4 + a5*Tmod)
        
        p = (a0 + a1*x[0] + a2*(x[0]**2))*(a3 + a4*x[1]) 
        return p
    
    params, pcov = curve_fit(_matiec,x,power)
    ermodule = {"a0":params[0], "a1":params[1], "a2":params[2], "a3":params[3], "a4":params[4]}
    
    print(ermodule)
    
    return ermodule





































        

