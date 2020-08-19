# Fill a few variables 
import numpy as np

def fill_values(Zci, deltai, gti):

    global Zc
    global delta
    global gas_type
    global gamma 

    Zc = Zci
    delta = deltai
    gas_type = gti
    gamma = np.float64(1.) + delta

    return 
