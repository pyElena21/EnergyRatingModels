'''
Gives the annual yield of a module

'''

import pvmodels as pvm
import pvtherm as pvt

try: 
    import pandas as pd
except ImportError():
    raise ImportError("You need to install pandas --- Try: pip install pandas")



try: 
    import numpy as np
except ImportError():
    raise ImportError("You need to install numpy")



def pvOutput(datafile, ermodule, Pstc, col_names = ['datetime', 'Geff','Tmod'],  area = 1, thrmcoeff = None, output_file = None,
             mod_temp = True, modelpv = 'iec61853_1', modeltrm = 'therm61853_1', thermcoeff = None, aggr_result = None,):
    
    '''
    Parameters
    ------------
    
    datafile: At the moment it's a file with 3 or four 4 columns of datetime, irradiance, temperature 
                and/or wind speed. If wind speed is not provided but required for modelling an error is 
                raised.
    
    ermodule: A dictionary with the characterisation coefficients
                Also check pvchar
    Pstc : Nominal capacity of the module
    
    col_names: A list of strings
                This defines the column names as appeared in the input CSV file for easy parsing
    
    
    area : The area of the module 
    
    thrmcoeff : The coefficients of the thermal model
                Also check pvtherm
    output_file: csv file
                The output path of the file if it is not None
    
    mod_temp: Boolean
            if True it means that module temperature is given otherwise it is ambient temperature given and 
            module temperature will be calculated with one of the chosen models
            
    modelpv: String. it must be one of the following:
              
              pvsat 
              esti 
              zenit 
              matrixmet   
              iec61853
              
              Also see pvmodels
            
    modeltrm : String. it must be one of the following:
                servantThermal 
                rossThermal 
                noctThermal 
                kingThermal 
                therm61853_1
               
               Also see pvtherm
                
    thermcoeff : The coefficients of the thermal model
    
    
    aggr_result: If None it means there is no aggregation applied to the final result. 
                At the moment this functionality is not enabled.
    
    
    Returns
    -------
    
    Power output in W
    
    NOTE: Currently the function does not do any aggregation.
    The results are given in the same resolution as irradiance and temperature
    
    
    
    '''
    
    
    # Read from file 
    
    data = pd.read_csv(datafile, index_col = None)
    
    # do not consider NAs
    
    print('Dropping NA values...')
    data.dropna(inplace = True)
    
    # check datetime data
    
    datetime = data[col_names[0]]
    
    # converting datetime to a pandas datetime object
    if not isinstance(datetime, pd.DatetimeIndex):
        datetime = pd.to_datetime(datetime, box = True) 

       
    # assign each column
    Geff = data[col_names[1]]
    
    if not mod_temp:
        
        Tamb = data[col_names[2]]
        
        if  modeltrm == 'therm61853_1':
            
            try:
                Ws = data[col_names[4]]
                Tmod = pvt.therm61853_1(thrmcoeff, Tamb, Geff, Ws)
                
            except (IndexError,KeyError):
                
                raise Exception ('You must provide wind speed for this model: Wind column name or series was not found')
               
        
        elif modeltrm == 'servantThermal':
            
            try:
                Ws = data[col_names[4]]
            
                Tmod = pvt.servantThermal(thrmcoeff, Tamb, Geff, Ws)
            
            except (IndexError,KeyError):
                
                raise Exception ('You must provide wind speed for this model: Wind column name or series was not found')
               
            
        elif modeltrm == 'rossThermal':
        
            Tmod = pvt.rossThermal(thrmcoeff, Tamb, Geff)
            
        elif modeltrm == 'noctThermal':
            
            Tmod = pvt.noctThermal(thrmcoeff, Tamb, Geff)
            
        elif modeltrm == 'kingThermal': 
            
            try:
            
                Ws = data[col_names[4]]
                Tmod = pvt.kingThermal(thrmcoeff, Tamb, Geff, Ws)
            
            except (IndexError,KeyError):
                
                raise Exception ('You must provide wind speed for this model: Wind column name or series was not found')
               
        else:
            raise Exception('You need to add module temperature or choose a thermal model from pvtherm module')
            
    else:
        Tmod = data[col_names[2]]
         
    if modelpv == 'pvsat': 
    
        pmax = pvm.pvsat(ermodule, Geff, Tmod, area)
        
    elif modelpv == 'esti':
        
        pmax = pvm.estier(Pstc, ermodule, Geff, Tmod)
   
    elif modelpv == 'zenit':
        
        pmax = pvm.zenit(ermodule, Geff, Tmod)
    
    elif modelpv == 'matrixmet':
        
        Tcell = data['Tcell']
        pmax = pvm.matrixmet(ermodule, Geff, Tmod, Tcell)
    
    elif modelpv == 'iec61853':
    
        pmax = pvm.iec61853_1(ermodule, Geff, Tmod)
    else:
        
        raise Exception('you must choose a valid pv model from pvmodel module')
    
    
  
    data['power'] = pmax
    #result = pd.DataFrame({'datetime':datetime, 'G(W/m2)':Geff, 'Tmod (oC)':Tmod, 'power':pmax})
    
    
    
    if output_file is None:
        
        print(data)
        
    else:
        
        data.to_csv(output_file)

    return data



# for adding aggregated results: To be added in later
def _findRes(DatetimeIndex):
    
    iter = 1
    res = None
    while res is None:
   
        i = np.random.randint(0,DatetimeIndex.shape[0]-1) 
        dates = pd.DatetimeIndex([DatetimeIndex[i-1],DatetimeIndex[i], DatetimeIndex[i+1]])  
        res = pd.infer_freq(dates)
        iter+= 1
        if iter > 50:
            
            print('No frequency found')
            break   
    return res


def _spltString(resolution):           
    
    import string 
    
    for let in resolution:
         
        if let in list(string.ascii_letters):
             
            num = resolution.strip(let)
            
            if num == '':
                
                num = 1.0          

    return  int(num), let













