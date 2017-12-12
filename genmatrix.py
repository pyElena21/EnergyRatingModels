'''

This module is used to generate matrices per the IEC 61853-1 to test the ER models

'''
try:
    from pvlib import pvsystem
except ImportError:
    raise ImportError("Latest version of pvlib-python is needed for this module")
    

try:
    import pandas as pd
except ImportError:
    raise ImportError("You need to install Pandas for this module")


def chooseCECModuleCoeffs(file = 'sam-library-cec-modules-2015-6-30.csv', module_model = '1Soltech_1STH_245_WH'):
    
    '''
    This function helps choose a module from CEC module database and its output can be used to generate matrices for a 
    large collection of modules.
    
    Parameters
    ----------
    
    file:         String  
                  A path to file that contains the module database. Could be the original or a (shorter) copy of that 
    
    module_model: String 
                  A module name from the database  
    
    
    Returns
    -------
    
    module_params_stc:  A dictionary of the STC parameters of the modules needed for the one diode model
    
    Also see generateMatrix
    
    ''' 
    
    module_coeffs = pvsystem.retrieve_sam(name = "CECMod",samfile = file)
    
    try:
        test_module = module_coeffs[module_model]
        module_params_stc = test_module.to_dict()
        
    except KeyError:
        raise KeyError('The chosen model does not exist or is typed incorrectly')
    
    return module_params_stc
    
    

def generateMatrix(module_params_stc={}, file_input_data = None, col_names = ['G','T'], EgRef = 1.121, dEgdT = -0.000126):
    
    
    '''
    Generates the 61853-1 matrix based on IV diode modelling and the module reference parameters.  
    
    Parameters
    ----------
    
    file input_data: CSV file 
    
    It is the data given in a csv file format having two columns 1: Irradiance and 2 temperature.
    
    Col_names :  A list of strings
    
    This defines the column names as appeared in the CSV file
    
    module_params_stc : default None 
                        If None then example parameters are provided for this
                        The module parameters at STC (standard testing conditions)
                        a_ref  : Diode modified ideality factor
                        
                        Calculated by n*Ns*k*T_ref/q: n diode ideality factor, Ns number of cells in series, k boltzmann constant,
                                                    T_ref temperature at STC, q electron charge. 
                        
                        IL_ref  : Light (photo-)current 
                        I0_ref  : Diode saturation current
                        Rsh_ref  : Shunt resistance (parallel)
                        Rs_ref  : Series resistance 
                        alpha_sc : Temperature coefficient for short circuit current (needed for translation to other conditions)
                        V_oc_ref : open circuit voltage at STC
                        
    In this module the following conditions are taken into account                    
    conditions : The G-T matrix of 61853 -1. 
                 G (irradiance),T(temperature) data taken from a csv file. 
                 Be careful of the data, they should have titled columns as 'G','T' 
                 
                 
                 
                 if none is given then the following matrix is used for 22 conditions of G,T: 
                 
                   G\T| 15  25  50  75
                   -------------------------
             1100    |  -    x   x   x
             1000    |  x    x   x   x
              800    |  x    x   x   x
              600    |  x    x   x   x
              400    |  x    x   x   -
              200    |  x    x   -   -
              100    |  x    x   -   -
                 
    EgRef:  float in eV
            The bandgap of the material at reference conditions 
    
    dEgdT:  float
            gradient of Eg with temperature
            
    
    
    Returns
    ------
    
    matrixdf: A Pandas dataframe with 3 columns of irradiance, module temperature and maximum power
    
    '''
    
    AMref = 1

    k =1.38E-23 
    q = 1.602E-19
 
    if not module_params_stc:
        
        # for a 60 cell module P_ref = 230 W
        module_params_stc['I_L_ref']   = 8.44  
        module_params_stc['I_o_ref']   = 1.04E-09
        module_params_stc['R_sh_ref']  = 319
        module_params_stc['R_s']   = 0.367
        module_params_stc['a_ref'] = 1.05*60*298.0* k/q      #n*Ns*k*T_ref/q
        module_params_stc['V_oc_ref'] = 36.9
        module_params_stc['alpha_sc'] = 0.005

        
    alpha_isc = module_params_stc['alpha_sc']
    
    
        
    
    if file_input_data is None:
        
        print('Generating default matrix ...')
        
        G = [100,200,400,600,800,1000,1100]
        T = [15,25,50,75]
       
        # 61853 -1 matrix  22 conditions
        
        conditions = [(G[0],T[0]),(G[0],T[1]),
                      (G[1],T[0]),(G[1],T[1]),
                      (G[2],T[0]),(G[2],T[1]), (G[2],T[2]),
                      (G[3],T[0]), (G[3],T[1]),(G[3],T[2]),(G[3],T[3]),
                      (G[4],T[0]), (G[4],T[1]),(G[4],T[2]),(G[4],T[3]),
                      (G[5],T[0]), (G[5],T[1]),(G[5],T[2]),(G[5],T[3]),
                      (G[6],T[1]), (G[6],T[2]),(G[6],T[3])]
       
        
        df =  pd.DataFrame(conditions, columns = ['G','T'])
        
        G = df['G']
        T = df['T'] 
       
    else: # write the given data frame as list of tuples               
        
        print('Make sure the format is as requested...')
        
        df = pd.read_csv(file_input_data, index_col = None)
          
        G = df[col_names[0]]
        T = df[col_names[1]]
        
        
        
    module_params_list = pvsystem.calcparams_desoto(G, T, alpha_isc, 
                                        module_params_stc, EgRef, dEgdT, AMref)
        
    module_params = {
                         'I_L': module_params_list[0], 
                         'I_o':module_params_list[1], 
                         'R_s': module_params_list[2], 
                         'R_sh': module_params_list[3], 
                         'nNsVth': module_params_list[4]
                         }
        
    dfparams = pvsystem.singlediode(module_params_stc, 
                                  module_params['I_L'], 
                                  module_params['I_o'], 
                                  module_params['R_s'], 
                                  module_params['R_sh'], 
                                  module_params['nNsVth']
                                  )
    pmax = dfparams['p_mp'] 
    imp = dfparams['i_mp']
    vmp = dfparams['v_mp']
    
    #print(dfparams['v_mp'])
        
    matrixdf = pd.DataFrame({'G (W/m2)':G,'T(oC)': T, 'P(W)': pmax})
    
    print(matrixdf)
    return matrixdf
    
    



