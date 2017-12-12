'''

Steps followed to run the code internally.

Any input data can be changed from here

'''

from genmatrix import generateMatrix, chooseCECModuleCoeffs
import pvyield
import pvchar as pvc
import os

def main(datafile='C:\\Users\\elek2_backup\\Helena\\DATA_WORKSPACES\\Photoclass\\outdoor_data_example.csv'):
    
    ''' 
    
    That is an example of how to use the code.
    For more details refer to the corresponding modules

    datafile : csv file to save the ouput data
    
    Also See: generateMatrix, pvchar and pv yield
    
    
    '''
    # Generate a matrix
    
    # choose a module from the CEC module database

    os.chdir('..')
    samfile = os.path.join(os.getcwd(),'Test Files//sam-library-cec-modules-2015-6-30.csv')
    
    module_coeffs = chooseCECModuleCoeffs(file = samfile, module_model = '1Soltech_1STH_245_WH')
    
    # Generate matrix (default values are for Silicon)
    
    matrix = generateMatrix(module_coeffs, file_input_data = None, col_names = ['G','T'], EgRef = 1.121, dEgdT = -0.000126)
    
    # Interpolate matrix using a method - stc might be acquired from matrix
    
    P_stc = matrix.loc[matrix['G (W/m2)']==1000].loc[matrix['T(oC)']==25]['P(W)'].item() #matrix['T(oC)'] == 25.0]#['T(oC)'==25])
    
    # Choose the 61853-1 method for the coefficients
    
    Geff = matrix['G (W/m2)']
    Tmod = matrix['T(oC)']
    power = matrix['P(W)']
    
    ermodule = pvc.pvsatMat(Geff, Tmod, power, area = 1)
    
    # calculate pvyield
    
    # define your datafile (this input might change in the future datafile is added in function for convenience 
    # but is not very good practice)
    
    
    
    
    output = pvyield.pvOutput(datafile, ermodule, P_stc, area = 1, thrmcoeff = None, output_file = None,
             mod_temp = True, modelpv = 'pvsat', modeltrm = 'therm61853_1', thermcoeff = None, aggr_result = None)
    
    
    
    return output


if __name__ == '__main__':
    #pass
    main()
    
