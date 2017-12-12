'''

Author: Elena Koubli


This module comprises some of the most used PV thermal models in literature for the translation 
of ambient to module temperature. 

Check pvtester to have an idea how to run this module.

Main reference:
----------------

Segado, Patricia, Mora J. Carretero, Sidrach-de-Cardona, Mariano, Models to predict the operating 
temperature of different photovoltaic modules in outdoor conditions,Progress In Photovoltaics_ 
Research And Applications,2014



'''

try:
    import numpy as np
except ImportError():
    raise ImportError("You need to install numpy")


def servantThermal(thrmcoeff,Tamb,Geff,Ws):
    
    ''' 
    Parameters
    ------------
    
    thrmcoeff: dict 
              A dictionary containing the coefficients of the module
              It should contain the three fields required for this model    
               
    Geff:     numpy float or 1d array or list (Watt/m2)
              Effective irradiance i.e. irradiance reaching the surface of the module
              Either calculated or measured. The values need to be *positive*
               
    Tamb:     numpy float or 1d array or list
              ambient temperature.
              
              
    Ws:       numpy float or 1d array or list
              Windspeed
    
    
              
              
    Returns: Module temperature 
    
    Tmod:     or 1d array or series
              module temperature.
              
              
    
    '''   
    
    try:
        
        D = thrmcoeff['D']
        E = thrmcoeff['E']
        F = thrmcoeff['F']
    
    except KeyError:
        
        raise Exception('the required coefficients for this model are 3: D, E, F')
        
    Tmod = Tamb + D*Geff*(1+E*Tamb)*(1-F*Ws)
    
    return Tmod 
    


def rossThermal(thrmcoeff,Tamb,Geff):
    
    ''' 
    Parameters
    ------------
    
    thrmcoeff: dict 
              A dictionary containing the coefficients of the module
              It should contain the three fields required for this model    
               
    Geff:     numpy float or 1d array or list (Watt/m2)
              Effective irradiance i.e. irradiance reaching the surface of the module
              Either calculated or measured. The values need to be *positive*
               
    Tamb:     numpy float or 1d array or list
              ambient temperature (Celsius).
              
              
              
    Returns: Module temperature 
    
    Tmod:     or 1d array or list
              module temperature (in celcius).
            
              
    Note: The coefficients depend on the module mounting as shown below
    
    

    1- well cooled (k= 0.02)
    2- free standing (k=0.0208)
    3- flat on roof  (0.026)
    4- not so well cooled (k=0.0342)
    5- transparent PV (k=0.0455)
    6- facade integrated (k=0.0538)
    7- on sloped roof (k=0.0563)
                  

    
    References:
    -----------
    
    Ross, Jr, R.G., Interface design considerations for terrestrial solar cell modules,12th Photovoltaic 
    Specialists Conference, Baton Rouge, La. USA, 1976
    
              
    '''
    
    ross = thrmcoeff['ross']  
    Tmod = Tamb + ross*Geff
    
    
    return Tmod 




def noctThermal(thrmcoeff, Tamb, Geff, units = "celcius"):
    
    ''' 
    Parameters
    ------------
    
    thrmcoeff: dict 
              A dictionary containing the coefficients of the module
              It should contain the three fields required for this model    
               
    Geff:     numpy float or 1d array or list (Watt/m2)
              Effective irradiance i.e. irradiance reaching the surface of the module
              Either calculated or measured. The values need to be *positive*
               
    Tamb:     numpy float or 1d array or list
              ambient temperature.
              
              
              
    Returns: Module temperature 
    ---------
    
    Tmod:     or 1d array or list
              module temperature (in celcius).
              
    '''          
    
    noct = thrmcoeff['noct']
    
    if units.lower() == "celcius":
        
        Tmod = Tamb + (Geff/800)*(noct-20)
         
    elif units.lower() == "kelvin":
        
        Tmod = Tamb + (Geff/800)*(noct-293.0)
     
    
    return Tmod # NOCT model in Celsius
    



def kingThermal(thrmcoeff, Tamb, Geff, Ws):
    
    ''' 
    Parameters
    ------------
    
    thrmcoeff: dict 
              A dictionary containing the coefficients of the module
              It should contain the three fields required for this model    
               
    Geff:     numpy float or 1d array or list (Watt/m2)
              Effective irradiance i.e. irradiance reaching the surface of the module
              Either calculated or measured. The values need to be *positive*
               
    Tamb:     numpy float or 1d array or list
              ambient temperature (in Celsius).
              
    Ws:       numpy float or 1d array or list 
              Wind speed
              
    Returns: Module temperature 
    
    Tmod:     or 1d array or list
              module temperature (in celcius).
              
    
    Taken from pvlib- python
    -------------------------
    Typical parameters:
    -------------------
    
    Module Type                           Mount              a            b         deltat

    Glass/cell/glass                   Open rack            -3.47     -0.0594          3
    Glass/cell/glass                   Close roof mount     -2.98     -0.0471          1
    Glass/cell/glass                   Insulated back       TBD         TBD
    Glass/cell/polymer sheet           Open rack            -3.56     -0.0750          3
    Glass/cell/polymer sheet`          Close roof mount      TBD         TBD
    Glass/cell/polymer sheet           Insulated back        -2.81     -0.0455         0
    22X Linear Concentrator            Tracker               -3.23      -0.130         13
    
    
    
    
             
    Reference:
    ----------
    Kratochvil, Jay A Boyson, William Earl King, David L, Photovoltaic array performance model,
    SAND2004-3535, Sandia National Laboratories (SNL), 2004
    
    '''
    
    try: 
                
        m = thrmcoeff['m']
        n = thrmcoeff['n']
    
    except KeyError:
        
        raise Exception('Two coefficients are needed for this model: m,n')
    
    
    Tmod = Tamb+Geff*np.exp(m+n*Ws)
    
    return Tmod

  

def kingModToCell(deltat, Geff, Tmod):
    
    ''' 
    Parameters
    ------------
    This equation returns cell temperature from the back of module temperature
    
    thrmcoeff: dict 
              A dictionary containing the coefficients of the module
              It should contain the three fields required for this model    
               
    Geff:     numpy float or 1d array or list (Watt/m2)
              Effective irradiance i.e. irradiance reaching the surface of the module
              Either calculated or measured. The values need to be *positive*
                         
     
    deltaT:  numpy float or 1d array or list
            Is taken for different configurations of module types (see documentation string of KingThermal)
              
    Returns: cell temperature 
    
    
    Tcell:     or 1d array or list
              module temperature (in celcius).
              
    
    
    *Check kingThermal for input coefficients*
    
              
    Reference:
    
    Kratochvil, Jay A Boyson, William Earl King, David L, Photovoltaic array performance model,
    SAND2004-3535, Sandia National Laboratories (SNL), 2004
    
    '''
    
    
    Tcell = Tmod + (Geff/1000)*deltat
    
    
    
    return Tcell



# Need to add the draft 61853 model here


def therm61853_1(thrmcoeff, Tamb, Geff, Ws):
    
    ''' 
    Parameters
    ------------
    This equation returns module temperature proposed in the IEC standard 61853 -2 (draft) 
    
    thrmcoeff: dict 
              A dictionary containing the coefficients of the module
              It should contain the three fields required for this model    
               
    Geff:     numpy float or 1d array or list
              Effective irradiance i.e. irradiance reaching the surface of the module
              Either calculated or measured. The values need to be *positive*
    
    Tamb:     numpy float or 1d array or list
              ambient temperature (in Celsius)
                         
    Ws:       numpy float or 1d array or list 
              Wind speed       
    
    Returns: module temperature 
    ---------
    
    Tmod:     numpy float or 1d array or list
              module temperature (in celcius).            
    
            
    Reference:
    
    D.Faiman "Assessing the outdoor operating temperature of photovoltaic modules" in prog. photovoltaics Res. Appl. pp 307-315, 2008
    
    '''
    
    try:
        
        u0 = thrmcoeff['u0']
        u1 = thrmcoeff['u1']
        
    except KeyError:
        
        raise Exception('Two coefficients are needed for this model: u0,u1')
    
    Tmod = Tamb + Geff*(u0 + u1*Ws)
    
    return Tmod



