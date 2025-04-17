# -*- coding: utf-8 -*-
"""
Created on Mon Apr 14 14:49:34 2025

@author: ethan
"""

import numpy as np

def ReadCombined(filename):
    """This function is designed to read and extract
    data from the combined MW-M31 data files.

    PARAMETERS
    ----------
        filename: 'string'
            The name of the combined data text file.
    OUTPUTS
    -------
        data: 'np.ndarray'
            The array of data extracted from the file.
            Includes mass, ptype, position, and 
            velocity data.
        
    """
    #generate the data array from the file
    data = np.genfromtxt(filename,dtype=None,names=True)

    return data