

class Module(object):
    
    '''
    This class defines a PV module object with all its essential characteristics
    
    '''

    def __init__(self, IL_ref, I0_ref, Rsh_ref, Rs_ref, alpha_isc, p_stc, N_series = 1, n_factor= None, a_mod_factor = None, type = None, **kwargs):
        
        
        '''
        Requires essential modelling parameters and STC data and returns a dictionary which 
        can be updated with new characteristics
        
        '''
        
        self.IL_ref = IL_ref
        self.I0_ref = I0_ref
        self.Rsh_ref = Rsh_ref
        self.Rs_ref = Rs_ref
        self.alpha_isc = alpha_isc
        self.p_stc = p_stc
        self.n_factor = n_factor # diode ideality factor
        self.a_mod_factor = a_mod_factor # modified diode ideality factor
        self.N_series = N_series
        self.type = type 
                
    
    def getParams(self, **kwargs):
            
        module_params = {'IL_ref': self.IL_ref,'I0_ref': self.I0_ref,
                         'Rs_ref':self.Rs_ref,'Rsh_ref': self.Rsh_ref, 
                         'n':self.n_factor, 'a':self.a_mod_factor}
          
        module_params.update(**kwargs)
        
        return module_params
            
            
        
        
        
       
    

 







    