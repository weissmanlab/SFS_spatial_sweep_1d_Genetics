import numpy as np
import math
import random
import matplotlib.pyplot as plt
import cv2
import os
import matplotlib.cm as cm
import itertools

data = np.loadtxt('L=50_N=100_s=0.500000_m=0.100000_tfinal=3000_l0=4_0.txt')
time_pts = data.shape[0]   #Number of frames
dimensions = data.shape[1] #Dimension of each frame
data = data.reshape(time_pts/dimensions, dimensions, dimensions)  ## ((tfinal, L, L))

time = data.shape[0]
print(time)
img_array = []


""""Plotting"""""
for i in range(0,time,10): 
        print(i)
        plt.figure(i)
        plt.imshow(data[i])
        plt.colorbar()
        #plt.xlabel('Deme Number')
        #plt.ylabel('Frequency of Beneficial Mutation in deme')
        plt.savefig(str(i))
        plt.close(i)
        img = cv2.imread(str(i)+".png")
        height, width, layers = img.shape
        size = (width,height)
        img_array.append(img)
        os.remove(str(i)+".png")
    
fourcc = cv2.cv.CV_FOURCC(*'XVID')  #Note: Check version compatiblity

video = cv2.VideoWriter('2-D_wave_L=50_N=100_s=0.50_m=0.1_tfinal=198_l0=4.avi', fourcc, float(1), (width, height)) 
for i in range(len(img_array)):
    video.write(img_array[i])
video.release()
cv2.destroyAllWindows()
