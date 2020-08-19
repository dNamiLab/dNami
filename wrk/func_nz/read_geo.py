import numpy as np

def read_geo(rescale=True):

    #Read in the information 
    #print(' Reading geometry information ...')
    geo_info = np.loadtxt('geo.txt',dtype = np.float64)

    #Geometry information
    A_L = geo_info[0]
    dA_x = geo_info[1]
    lx = geo_info[2]
    
    if rescale:
        #Redim dAdx so that Lx = 1
        dA_x = dA_x * lx

    return A_L, dA_x 


def read_lx():

    #Read in the information 
    geo_info = np.loadtxt('geo.txt',dtype = np.float64)

    #Geometry information
    A_L = geo_info[0]
    dA_x = geo_info[1]
    lx = geo_info[2]

    return lx 
