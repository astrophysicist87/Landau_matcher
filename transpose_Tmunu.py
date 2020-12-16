#! /usr/bin/env python

import numpy as np
import sys
#import os.path

if len(sys.argv) != 5:
    print "Usage: python transpose_Tmunu.py [filename] [x grid size] [y grid size] [new filename]"
else:
    filename = sys.argv[1]
    xGridSize = int(sys.argv[2])
    yGridSize = int(sys.argv[3])
    
    #newfilename = os.path.dirname(filename) + '/transposed_' + os.path.basename(filename) 
    newfilename = sys.argv[4]
    
    data = np.loadtxt(filename)
    dataShape = data.shape
    data = data.reshape([xGridSize, yGridSize, dataShape[-1]])
    
    # transpose x and y axes
    data=np.transpose(data,axes=(1,0,2))
    
    # and save the result
    np.savetxt(newfilename, data.reshape(dataShape), fmt='%12.8g')
    print 'Saved transposed Tmunu to', newfilename