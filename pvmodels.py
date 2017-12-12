''' 
Last updated 12-10-2017 by Elena Koumpli

This file contains module models that have been proposed in Energy Rating round robins (RR) and in the IEC
Standard 61853-1 BUT now has been updated to more parametric models. 

The energy rating equations found in here are for the *modelling procedure* 

Also see pvchar to use these models for the *characterisation procedure* i.e. the extraction of the model 
coefficients. 

-----------
Comments
-----------

The methods can be divided into two groups:

1) Group one (CREST, ECN, UU and INES) determines separately, 
the module efficiency (for 25°C and different irradiance levels) denoted as η(G,25°C) and 
the temperature coefficient as a constant TC or a curve TC(G) and treat these effects as 
independent of each other. 
2) The second group (SUPSI, JRC, H2M and ISE) describes the whole power surface P(G,T) or η(G,T) 
as a single equation.

*Only the methods of the second group are described here due to ambiguousness of models in group 1*

-----------
Main references
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

[4] IEC 61853-1, Photovoltaic (PV) module performance testing and energy rating, 2011

Dependecies
------------

You will find pandas and numpy in pre-compiled packages at http://www.lfd.uci.edu/~gohlke/pythonlibs/


'''

try:
    import numpy as np
except ImportError():
    raise ImportError("You need to install numpy")

try: 
    import pandas as pd
except ImportError():
    raise ImportError("You need to install pandas --- Try: pip install pandas")



def pvsat(ermodule,Geff,Tmod,area, **kwargs):
    
    ''' 
    Origin: H2M
    Coefficients: four 
    Result: efficiency and/or power
    
    Parameters
    ------------
    
    ermodule: dict 
              A dictionary containing the coefficients of the module
              It should contain the four fields required for this model    
    
    
    data :  A dataframe containing the following columns
    
        
        Datetime: Datetime series
        
        
        Geff:     Effective irradiance i.e. irradiance reaching the surface of the module
                  Either calculated or measured. The values need to be *positive*
                   
                   
        Tmod:     module temperature in Celcius.
                  Either calculated(see pvtherm) or measured
                  
    Area:     numpy float
              the area of the module
              
              
    Returns:  A dataframe of datetime, irradiance, temperature and power float or series
    --------
    
    Note: Coefficient symbols are consistent to the source
    
    Also check: pvchar 
    
    '''        
    
    #this will change in the future - depends on other modules in the package 
    
    try: 
       
        # Naming the coefficients
        a1 = ermodule['a1'] 
        a2 = ermodule['a2']
        a3 = ermodule['a3']
        a = ermodule['a']
        
    except KeyError: 
    
        raise Exception('pvsat method needs 4 coefficients to run: a,a1,a2,a3') 
    
        
    Geff = np.array([Geff])
    Geff[Geff<0.0] = 0.0
    Tmod = np.array([Tmod])
    
      
    # power is calculated at each point and a new dataframe is returned
                 
    power = (a1 + a2*Geff +a3*np.log(Geff)) * (1+a*(Tmod-25))*area*Geff 
    power[power<0.0] = 0.0
   
         
    return power[0]





def estier(Pstc,ermodule,Geff, Tmod, **kwargs):

    
    ''' 
    Origin: JRC (Joint Research Centre)
    Coefficients: six 
    Result: power
    
    Parameters
    ------------
    
    ermodule: dict 
              A dictionary containing the coefficients of the module
              It should contain the six fields required for this model    
        
        
        Geff:     Effective irradiance i.e. irradiance reaching the surface of the module
                  Either calculated or measured. The values need to be *positive*
                   
                   
        Tmod:     module temperature in Celcius.
                  Either calculated(see pvtherm) or measured
              
    Pstc:     numpy float. 
              The maximum power at standard testing conditions STC.
              
              
    Returns:  A dataframe of datetime, irradiance, temperature and power float or series
    --------
    
    Note: Coefficient symbols are consistent to the source
    
    Reference:
    ----------
    Huld, Thomas, Friesen, Gabi, Skoczek, Artur, Kenny, Robert P., Sample, Tony, Field, Michael
    Dunlop, Ewan D.A, power-rating model for crystalline silicon PV modules, Solar Energy Materials 
    and Solar Cells, 2011
    
    
    '''        

    try: 
        
        # Naming the coefficients
        
        k1 = ermodule['k1'] 
        k2 = ermodule['k2']
        k3 = ermodule['k3']
        k4 = ermodule['k4']
        k5 = ermodule['k5']
        k6 = ermodule['k6']
        
        
    except KeyError: 
    
        raise Exception('ESTI-ER method needs the following 6 coefficients to run: k1,k2,k3,k4,k5,k6' ) 
    
    
    Geff = np.array([Geff])
    Geff[Geff<0.0] = 0.0
    Tmod = np.array([Tmod])
        
    
    
    power = Pstc*0.001*Geff*(1+k1*np.log(0.001*Geff)+k2*((np.log(0.001*Geff))**2)+
            (Tmod-25)*(k3+(k4)*(np.log(0.001*Geff)+k5*((np.log(0.001*Geff)**2)))+ k6*((Tmod-25)**2)))       
                 
    power[power<0.0] = 0.0
    
    
    return power[0]




def zenit(ermodule,Geff, Tmod, **kwargs):
    
    
    ''' 
    Origin: ISE
    Coefficients: five 
    Result: power (power surface)
    
    Parameters
    ------------
    
    ermodule: dict 
              A dictionary containing the coefficients of the module
              It should contain the six fields required for this model    
        
        
        Geff:     Effective irradiance i.e. irradiance reaching the surface of the module
                  Either calculated or measured. The values need to be *positive* (Watt/m2)
                   
                   
        Tmod:     module temperature in Celcius.
                  Either calculated(see pvtherm) or measured
              
                 
              
    Returns:  power float or series
    --------
    
    '''
    
    try: 
       
        # Naming the coefficients
        a = ermodule['a'] 
        b = ermodule['b']
        c = ermodule['c']
        d = ermodule['d']
        e = ermodule['e']
            
    except KeyError: 
    
        raise Exception('Zenit method needs 5 coefficients to run: a,b,c,d,e') 
    
    
    Geff = np.array([Geff])
    Geff[Geff<0.0] = 0.0
    Tmod = np.array([Tmod])
    
    
    psi = ((np.log(Geff+e))**2)/(Geff+1)
    power = a*(Geff**2) + b*np.log(Geff+1) * Geff+c*(psi-1)*Geff + d*(Tmod-25)
    power[power<0.0] = 0.0
    
    
    return power[0]


def matrixmet(ermodule,Geff,Tmod,Tcell,**kwargs):
    
    ''' 
    Origin: SUPSI
    Coefficients: six 
    Result: I, V, power 
    
    Parameters
    ------------
    
    ermodule: dict 
              A dictionary containing the coefficients of the module
              It should contain the six fields required for this model    
               
  
        
        
        Geff:     Effective irradiance i.e. irradiance reaching the surface of the module
                  Either calculated or measured. The values need to be *positive* (Watt/m2)
                   
                   
        Tmod:     module temperature in Celcius.
                  Either calculated(see pvtherm) or measured
              
        Tcell:    numpy float or list
                  cell temperature
               
    
              
              
    Returns:  Power float or series
    --------
    
    Note: Coefficient symbols are consistent to the source
    

    '''
    
    try: 
    
        # Naming the coefficients
        im_stc = ermodule['im_stc'] 
        a_im = ermodule['a_im']
        vm_stc = ermodule['vm_stc']
        b_vm = ermodule['b_vm']
        c0 = ermodule['c0']
        c1 = ermodule['c1']
        
        
    except KeyError: 
    
        raise Exception('matrix method needs 6 coefficients to run:im_stc,a_im,vm_stc,b_vm,c0,c1') 
    
    
    Geff = np.array([Geff])
    Geff[Geff<0.0] = 0.0
    Tmod = np.array([Tmod])
    Tcell = np.array([Tcell])
    deltat = Tcell - Tmod 
    
    imax = (im_stc*Geff/1000) * (1+a_im*(deltat+Tmod-25))
    vmax = vm_stc+c0*np.log(Geff/1000) + c1*((np.log(Geff/1000))**2) + b_vm*(deltat+Tmod-25)
    
    power = imax*vmax
    power[power<0] = 0.0
    
    return power[0]
    
    


def iec61853_1(ermodule,Geff,Tmod, **kwargs):
    
    ''' 
    Origin: 61853-1
    Coefficients: varies with the degree of the polynomial regression --> here used 5
    2nd degree polynomial irradiance and (1D) linear temperature dependance for pmax.
     
    Result: power 
    
    Parameters
    ------------
    
    ermodule: dict 
              A dictionary containing the coefficients of the module
              It should contain the six fields required for this model    

        
        Geff:     Effective irradiance i.e. irradiance reaching the surface of the module
                  Either calculated or measured. The values need to be *positive* (Watt/m2)
                   
                   
        Tmod:     module temperature.
                  Either calculated(see pvtherm) or measured
    
    Returns
    --------
    
    power float or series
 
    
    '''

    try:
        
        a0 = ermodule['a0'] 
        a1 = ermodule['a1']
        a2 = ermodule['a2']
        a3 = ermodule['a3']
        a4 = ermodule['a4']
        #a5 = ermodule['a5']
    
    except KeyError:
        
        raise Exception('Current interpolation method needs 5 coefficients to run: a, a1, a2, a3, a4')
    
    Geff = np.array([Geff])
    Geff[Geff <0.0] = 0.0
    Tmod = np.array([Tmod])
    
    # 2D polynomial function for irradiance and 1D linear for temperature dependance @  Blago Mihailov
    # Caveats: Waht happens at non-linear devices - the dataset needs to be sliced to different levels of irradiance.
    
    power = (a0 + a1*Geff + a2*(Geff**2))*(a3 + a4*Tmod)
    power[power<0.0] = 0.0
    
    
    return power[0]
     


def rplaton(ermodule,Geff,Tmod, **kwargs):
    
    '''
    Model described in [5] and used in fault detection
     
    Result: current, voltage, power at maximum power point 
    
    Parameters
    ------------
    
    ermodule: dict 
              A dictionary containing the coefficients of the module
              It should contain the six fields required for this model    

        
        Geff:     Effective irradiance i.e. irradiance reaching the surface of the module
                  Either calculated or measured. The values need to be *positive* (Watt/m2)
                   
                   
        Tmod:     module temperature.
                  Either calculated(see pvtherm) or measured
                  
        
        
        
    
    Returns
    --------
    
    power series 
 
    
    
    [5] R. Platon, J. Martel, N. Woodruff, et al., Online Fault Detection in PV Systems,
    in IEEE TRANSACTIONS ON SUSTAINABLE ENERGY
    
    
    '''
    
    try:
        
        a0 = ermodule['a0'] 
        a1 = ermodule['a1']
        a2 = ermodule['a2']
        a3 = ermodule['a3']
        
    
    except KeyError:
        
        raise Exception('Current interpolation method needs 4 coefficients to run: a, a1, a2, a3')
    
    
    Geff = np.array([Geff])
    Geff[Geff <0.0] = 0.0
    Tmod = np.array([Tmod])
    

    
    power = Geff*(a0 + a1*Geff +a2*np.log(Geff))*(1+a3*(Tmod-25))
    power[power<0.0] = 0.0
    

    return  power[0]


    


