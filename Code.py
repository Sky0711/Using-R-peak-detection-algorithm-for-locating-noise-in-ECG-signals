

import sys
import numpy as np
import scipy.signal
import scipy.ndimage

# Function for locating the R-peaks
def detect_beats(
		ecg,	# The raw ECG signal
		rate,	# Sampling rate in HZ
		# Window size in seconds to use for 
		ransac_window_size=5.0,
		# Low frequency of the band pass filter
		lowfreq=5.0,
		# High frequency of the band pass filter
		highfreq=15.0,
		):
	

	ransac_window_size = int(ransac_window_size*rate)

	lowpass = scipy.signal.butter(1, highfreq/(rate/2.0), 'low')
	highpass = scipy.signal.butter(1, lowfreq/(rate/2.0), 'high')
	# TODO: Could use an actual bandpass filter
	ecg_low = scipy.signal.filtfilt(*lowpass, x=ecg)
	ecg_band = scipy.signal.filtfilt(*highpass, x=ecg_low)
	
	# Square (=signal power) of the first difference of the signal
	decg = np.diff(ecg_band)
	decg_power = decg**2
	
	# Robust threshold and normalizator estimation
	thresholds = []
	max_powers = []
	for i in range(int(len(decg_power)/ransac_window_size)):
		sample = slice(i*ransac_window_size, (i+1)*ransac_window_size)
		d = decg_power[sample]
		thresholds.append(0.5*np.std(d))
		max_powers.append(np.max(d))

	threshold = np.median(thresholds)
	max_power = np.median(max_powers)
	decg_power[decg_power < threshold] = 0

	decg_power /= max_power
	decg_power[decg_power > 1.0] = 1.0
	square_decg_power = decg_power**2

	shannon_energy = -square_decg_power*np.log(square_decg_power)
	shannon_energy[~np.isfinite(shannon_energy)] = 0.0

	mean_window_len = int(rate*0.125+1)
	lp_energy = np.convolve(shannon_energy, [1.0/mean_window_len]*mean_window_len, mode='same')
	#lp_energy = scipy.signal.filtfilt(*lowpass2, x=shannon_energy)
	
	lp_energy = scipy.ndimage.gaussian_filter1d(lp_energy, rate/8.0)
	lp_energy_diff = np.diff(lp_energy)

	zero_crossings = (lp_energy_diff[:-1] > 0) & (lp_energy_diff[1:] < 0)
	zero_crossings = np.flatnonzero(zero_crossings)
	zero_crossings -= 1
	return zero_crossings

# Function which detects anomalies in given signal	
def detect_anomaly(my_array1#the sample with higer frequency (like 360Hz or 256Hz)#
                  ,my_array2#the sample with 128HZ frequency#
                  ):
         NumberOf_Peaks = len(my_array1)
         i = 0; k = 0; NumberOfAnomalies = 0;
         peak_result=[]
                  
         while my_array1[i]!= my_array1[-1]:
          i=i+1
          k=k+1
          if my_array1[i]>(my_array2[k]+33): #go through both arrays and compare#
           print (' R Peak Number=%d  Sample Value=%d'%(i, my_array2[k]))
	   peak_result.insert(i, my_array2[k])
	   i=i-1
	   NumberOfAnomalies = NumberOfAnomalies+1
	    
	 
	 print ('\n Number of found anomalies=%d'%(NumberOfAnomalies))
	 print (' Number of detected R peaks = %d'%(NumberOf_Peaks))
	 return peak_result
	 
def all_anomalies(peak_256, peak_360):

         ListOfAllAnomalies = np.unique(peak_360 + peak_256)
         
         return ListOfAllAnomalies	  
	 
#Returns only R-peaks that contain anomaly	
def detect_R_peaks_with_anomaly(my_array1#the sample with higer frequency (like 360Hz or 256Hz)#
                  ,my_array2#the sample with 128HZ frequency#
                  ):
         NumberOf_Peaks = len(my_array1)              
         i = 0; k = 0;
         peak_result = []
         #unique_results = []
                  
         while my_array1[i]!= my_array1[-1]:
	  i=i+1
          k=k+1
          if my_array1[i]>(my_array2[k]+33): #search both arrays and compare#
	     peak_result.insert(i,[i, my_array2[k]])
	     i=i-1

         #unique_results = np.unique(peak_result)
	 return peak_result
	 
	 
#Returns the list of all R-peaks	
def sort_all_anomalous_R_peaks(AbnormalRpeaks_256 #Anomalous R-peks derived from 256HZ#
                              ,AbnormalRpeaks_360 #Anomalous R-peks derived from 360HZ#
                              ,ListOfAllAnomalies
                               ):
  
         res = AbnormalRpeaks_360 + AbnormalRpeaks_256
         raz = np.setdiff1d(res,ListOfAllAnomalies)

         print(raz)
         
 
                     
#Main Function	
if __name__ == '__main__':
	
	
	rate = float(360) #Frequency of ORIGINAL signal (360 in this case) 
	
	ecg = np.loadtxt(sys.stdin)

	my_array128 = detect_beats(ecg, 128) #Process ECG data with 128Hz and save them into array#
	my_array256 = detect_beats(ecg, 256) #Process ECG data with 256Hz and save them into array#
	my_array360 = detect_beats(ecg, 360) #Process ECG data with 360Hz and save them into array#
	peak_256=[]
	peak_360=[]

        Number = len(my_array360) # Number of found R-peaks in processd signal

        print '\n --------- \n 256HZ values \n ---------'
	peak_256 = detect_anomaly(my_array256,my_array128) #Searching for anomalies by using the 256HZ frequency#
	
	print '\n --------- \n 360HZ values \n ---------'
	peak_360 = detect_anomaly(my_array360,my_array128) #Searching for anomalies by using the 360HZ frequency#
	
	
	
	
	print '\n ----------------------------'
	TotalNumberOfAnomalies = len(all_anomalies(peak_256, peak_360))
	print ' \n TotalNumberOfAnomalies = %d \n'%(TotalNumberOfAnomalies)
        print ' Real number of R peaks in the signal = %d'%(Number) 
        
        
        AbnormalRpeaks_256 = detect_R_peaks_with_anomaly(my_array256,my_array128)
        AbnormalRpeaks_360 = detect_R_peaks_with_anomaly(my_array360,my_array128)
        print(AbnormalRpeaks_256)
        print(AbnormalRpeaks_360)
	
	AllAnomalies = all_anomalies(peak_256, peak_360)
	print(AllAnomalies)
	
	
	sort_all_anomalous_R_peaks(AbnormalRpeaks_256, AbnormalRpeaks_360, AllAnomalies)
        

	
	import matplotlib.pyplot as plt  #Prints R-peaks and marks the anomalies on the graph#
	dt = 1.0/rate
	t = np.linspace(0, len(ecg)*dt, len(ecg)) #Gives us durration of the signal, time(t).
	plt.plot(t, ecg)
	plt.xlabel("Time measured in seconds (s)",fontsize=13,fontweight='bold')
	plt.ylabel(" 10 mV range ",fontsize=13,fontweight='bold')
        
	plt.title('Sampling frequency = 256Hz',fontsize=18,ha='right',fontweight='bold',color='red')
	plt.suptitle('Sampling frequency = 360Hz',fontsize=18,ha='right', fontweight='bold',color='green')	
	plt.scatter(t[peak_256], ecg[peak_256],s=120,edgecolors='red',linewidth='3',color='white') # Prints 256Hz results
	plt.scatter(t[peak_360], ecg[peak_360],s=60,edgecolors='green',linewidth='3',color='white') #Prints 360Hz results
	plt.show()

	

