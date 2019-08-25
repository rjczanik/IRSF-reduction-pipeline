# -*- coding: utf-8 -*-
"""
Created on Mon Aug 12 09:08:29 2013

@author: robert
"""
import time
#import astropy as astp
#from astropysics import phot
#import astropysics as astpy
import pyfits
import numpy as np
import scipy as sp
import scipy.ndimage as ndimg
import pyraf as pyr
from pyraf import iraf
#from phot import aperphot
#from matplotlib.pylab import subplots,close
#from matplotlib import cm
import matplotlib.pyplot as plt
import pysex as pys
import os
import pp
import pywcs
#from ds9 import *
#print ds9_targets()
#d=ds9()

canstars = 150

standard_name = ['Pers9129', 'Pers9132', 'Pers9136', 'Pers9137','Pers9140']
standard_ra = [120.31291,126.39804,133.58916,138.96002,146.42983]
standard_dec = [-50.327051, -39.098971, -54.802914, -36.542614, -45.827126]
standard_J = [10.914, 11.949, 12.486, 11.153, 11.409]
standard_dJ = [0.006, 0.006, 0.008, 0.007, 0.011]
standard_H = [10.585, 11.669, 12.214, 10.891, 11.085]
standard_dH = [0.006, 0.005, 0.008, 0.007, 0.008]
standard_K = [10.496, 11.609, 12.142, 10.836, 11.022]
standard_dK = [0.009, 0.004, 0.011, 0.010, 0.012]

dates = ['050418','050422','050425','120104','120105','120107','120108','120109','120112','120113','120114','120115','120116']
fltrs = ['j','h','k']
#fltrs = ['k']

#fig,ax=subplots(1,1)
#ax.set_aspect('equal')
#ax.set_xlim(0,1240)
#ax.set_ylim(0,1240)
#ax.hold(True)

#FWHM = det_area = 5#int(raw_input("Please enter the detection area's radius:"))
#if det_area%2==0:
#    det_area=det_area+1

##############################FIND STARS ON THE IMAGE AND ELIMINATE BAD SOURCES#####

def stardetection(image_loc):

	global FWHM
	FWHM = det_area = 5
	image = pyfits.getdata(image_loc)
#Plotting the raw image
	plt.xlabel('X')
	plt.ylabel('Y')
	plot_image = np.log10(image)
	ax1.set_title('The raw scanned image')
	ax1.imshow(plot_image, cmap = plt.cm.binary)
	plt.draw()
	plt.show()
	image[np.isnan(image)] = -9999
	agtergrond=np.median(image)
	raw_stars = pys.run(image,params=['X_IMAGE','Y_IMAGE'],conf_args={'ANALYSIS_THRESH':2,'FILTER':'Y','FILTER_NAME':'default.conv','DETECT_TYPE':'CCD','DEBLEND_NTHRESH':32,'DEBLEND_MINCONT':0.005,'CLEAN':'Y','CLEAN_PARAM':1.0,'MASK_TYPE':'CORRECT','SEEING_FWHM':'1.2','STARNNW_NAME':'default.nnw','BACK_SIZE':32,'BACK_FILTERSIZE':3,'GAIN':5.0,'DETECT_THRESH':10,'PIXEL_SCALE':1.0,'THRESH_TYPE':'ABSOLUTE','DETECT_MINAREA':det_area})
	starsx=raw_stars[1].tonumpy()
	starsy=raw_stars[0].tonumpy()
#	for n in range(0, len(starsx)): #Filter the stars that are too close to the edge
#		if ( ( starsx[n] >= ( image[0,:] -10 ) ) or ( starsy[n] >= ( image[:,] - 10 ) ) ):
#			starsx = np.delete( starsx, n )
#			starsy = np.delete( starsy, n )
	starsx = np.delete( starsx, (starsx > np.size(image[ 0, : ]) - 13) )
	starsy = np.delete( starsy, (starsy > np.size(image[ :, 0 ]) - 13) )
	candidates = np.zeros((len(starsx),3),dtype=int)
	for n in range(0,len(starsx)-1):
		if image[starsy[n],starsx[n]] > agtergrond:
			candidates[n,0]=int(starsx[n])
			candidates[n,1]=int(starsy[n])
			candidates[n,2]=image[starsy[n]-1,starsx[n]-1]
	print candidates
	print np.size(candidates[:,0])
	print agtergrond
#Drawing the stars that were detected with Source Extractor as a whole collection in the candidates array
	ax2.imshow(plot_image, cmap = plt.cm.binary)
	ax2.set_title('Stars detected with SEXTRACTOR')
	ax2.scatter(starsx,starsy)
	plt.draw()
	plt.show()

#Using insert sorting to search for the brightest sources    
	end = len(starsx)
#    print end
	temp = np.empty(3,dtype=float)
	
	for i in range(1,end):
#        print 'i='+str(i)
		temp[0] = candidates[i,0]
		temp[1] = candidates[i,1]
		temp[2] = candidates[i,2]
		j = i - 1
		while ( j >= 0 ) and ( candidates[j,2] > temp[2] ):
			candidates[j+1,0] = candidates[j,0]
			candidates[j+1,1] = candidates[j,1]
			candidates[j+1,2] = candidates[j,2]
                        
			j = j-1
		candidates[j+1,0] = temp[0]
		candidates[j+1,1] = temp[1]
		candidates[j+1,2] = temp[2]

	print candidates 
    
#Slicing the array to the stars that will be used for the FWHM  and the PSF model
	candidates = np.copy(candidates[end-1.5*canstars:end-0.5*canstars:1,:])
	print candidates
	
#Purely experimental on stars chosen by hand:
	
#	if(filter == 'j'):
#	del candidates
#	candidates = np.zeros((19,3), dtype = float)
#	candidates[0,0] = 292
#	candidates[0,1] = 890
#	candidates[0,2] = image[890,292]
#	candidates[1,0] = 406
#	candidates[1,1] = 905
#	candidates[1,2] = image[905,406]
#	candidates[2,0] = 470
#	candidates[2,1] = 781
#	candidates[2,2] = image[781,470]
#	candidates[3,0] = 195
#	candidates[3,1] = 753
#	candidates[3,2] = image[753,195]
#	candidates[4,0] = 361
#	candidates[4,1] = 712
#	candidates[4,2] = image[712,361]
#	candidates[5,0] = 763
#	candidates[5,1] = 671
#	candidates[5,2] = image[671,763]
#	candidates[6,0] = 738
#	candidates[6,1] = 556
#	candidates[6,2] = image[556,738]
#	candidates[7,0] = 421
#	candidates[7,1] = 819
#	candidates[7,2] = image[819,421]
#	candidates[8,0] = 542
#	candidates[8,1] = 628
#	candidates[8,2] = image[628,542]
#	candidates[9,0] = 738
#	candidates[9,1] = 556
#	candidates[9,2] = image[556,738]
#	candidates[10,0] = 632
#	candidates[10,1] = 47
#	candidates[10,2] = image[47,632]
#	candidates[11,0] = 115
#	candidates[11,1] = 302
#	candidates[11,2] = image[302,115]
#	candidates[12,0] = 192
#	candidates[12,1] = 177
#	candidates[12,2] = image[177,192]
#	candidates[13,0] = 292
#	candidates[13,1] = 890
#	candidates[13,2] = image[890,292]
#	candidates[14,0] = 356
#	candidates[14,1] = 73
#	candidates[14,2] = image[73,356]
#	candidates[15,0] = 341
#	candidates[15,1] = 185
#	candidates[15,2] = image[185,341]
#	candidates[16,0] = 115
#	candidates[16,1] = 302
#	candidates[16,2] = image[302,115]
#	candidates[17,0] = 82
#	candidates[17,1] = 173
#	candidates[17,2] = image[173,82]
#	candidates[18,0] = 632
#	candidates[18,1] = 47
#	candidates[18,2] = image[47,632]
#	candidates[19,0] = 213
#	candidates[19,1] = 572
#	candidates[19,2] = image[572,213]
#	candidates[20,0] = 693
#	candidates[20,1] = 862
#	candidates[20,2] = image[862,693]
#	candidates[21,0] = 96
#	candidates[21,1] = 394
#	candidates[21,2] = image[394,96]
#	candidates[22,0] = 707
#	candidates[22,1] = 523
#	candidates[22,2] = image[523,707]
#	candidates[23,0] = 685
#	candidates[23,1] = 431
#	candidates[23,2] = image[431,685]
	
#	[[356,73], [341, 185], [115,302], [82, 173], [632,47]]

    #Searching for stars that are not too bright or too faint in the list:
#    for e in candidates:
#        print e
#        if ( ( e[2] < np.median(image) - 500 ) or ( e[2] > 28000 ) ):
#            np.delete(candidates, e)
#            print 'x '+str(e)
#    print candidates
#    can_dim = np.size(candidates[0,:])    
    
    #Determination of the extent of each of the 30 brightest stars' profiles
#    fwhm_arr = np.empty((canstars-1,2),dtype = int)
#    for j in range(canstars-1):
#        print j
#        dx = 0
#        limit_indi = False
#        noise_factor = 50
#        while (limit_indi == False):
#            print candidates[j,0]-dx, candidates[j,0]+dx, image[candidates[j,0],candidates[j,1]]
#            if (image[candidates[j,0]-dx,candidates[j,1]] > (agtergrond+noise_factor)) or (image[candidates[j,0]+dx,candidates[j,1]] > (agtergrond+noise_factor)):
#                limit_indi = False
#                dx = dx + 1
#            elif(  ):
                
#            else:
#                fwhm_arr[j,0] = dx
#                limit_indi = True


#Removing the candidates stars that are too close to the edge

	for candi in candidates:
#		print candi, candi[0], candi[1], np.size(image[:,0]), np.size(image[0,:])
		if ( ( candi[0] < 30 ) or ( candi[1] < 30 ) or ( candi[1] > ( np.size(image[:,0]) - 30 ) ) or ( candi[0] > ( np.size(image[0,:]) - 30 ) ) ):
			print 'XXX'
			candidates = np.delete(candidates, candi, 0)
	del candi
	print candidates, np.size(candidates[:,0])

#Determination of the extent for each of the candidate stars' profile. 
	fwhm_arr = np.zeros((2,canstars),dtype = int)
	counter = 0
	for j in candidates:
		print j
		dx = 0
		limit_indi = False
		noise_factor = 20
		print counter, j[1]-dx, j[0]-dx, image[j[1]-1,j[0]-1], image[j[1]-1,j[0]-1], np.size(image[0,:]), np.size(image[:,0])
		while ( ( limit_indi == False ) and ( dx < 12 ) ):
			if ( image[ j[1] -1 - dx, j[0] -1 - dx ] > ( agtergrond + noise_factor ) ):# or ( image[ j[1] -1 + dx, j[0] -1 + dx ] > ( agtergrond + noise_factor ) ):
				limit_indi = False
				dx = dx + 1
			elif( ( j[0] - dx ) < 50 or ( j[0] + dx ) > ( np.size(image[:,0]) - 50 ) or ( j[1] - dx ) < 50 or ( j[0] + dx ) > ( np.size(image[:,1]) - 50 ) or ( image[ j[1] -1 - dx, j[0] -1 ] == -9999 ) or ( image[ j[1] -1, j[0] -1 - dx ] == -9999 )):#or ( image[ j[1] -1 + dx, j[0] -1 ] == -9999 )  or ( image[ j[1] -1, j[0] -1 + dx ] == -9999 ) )
				np.delete(candidates, j)
#				fwhm_arr = np.copy(fwhm_arr[0:np.size(fwhm_arr)-1])
				limit_indi = True
			else:
				fwhm_arr[0,counter-1] = dx+1
				limit_indi = True
		counter = counter + 1
		del dx

#Determining the FWHM of each star's profile
	for i in range(0,counter-1):
		profile_range = fwhm_arr[0,i]
		print i,profile_range, candidates[ i, 2 ], (candidates[ i, 2 ] - agtergrond)/2
		for dy in range(0,profile_range):
#			print dy, image[candidates[i,1]-1,candidates[i,0]-1]-image[candidates[i,1]-1-dy,candidates[i,0]-1-dy], (candidates[ i, 2 ] - image[ candidates[ i, 1 ]-1 + profile_range, candidates[ i, 0 ]-1 + profile_range ])/2, candidates[i,2] - image[candidates[i, 1]-1 - profile_range, candidates[i,0]-1 - profile_range]
			if ( ( image[ candidates[i,1] - dy, candidates[i,0] - dy ] ) >= ( ( candidates[i,2] - image[candidates[i, 1] - profile_range, candidates[i,0] - profile_range] ) /2 ) ):
				fwhm_arr[1,i] = dy
	print fwhm_arr
	
	good_fwhm = np.where( fwhm_arr != 0 )
	fwhm_arr = np.copy(fwhm_arr[good_fwhm])
	for l in fwhm_arr:
		if (l == 0):
			fwhm_arr = np.delete(fwhm_arr, l)

	print fwhm_arr

	FWHM = np.median(fwhm_arr)*2.0#
	print fwhm_arr, FWHM

	ax3.set_title('FWHM candidates') #Plotting the stars that were used for the FWHM selection and the construction of the PSF
	ax3.imshow(plot_image, cmap = plt.cm.binary)
	ax3.scatter(candidates[:,0],candidates[:,1], color='orange')
	plt.draw()
	plt.show()


	if os.path.exists(datedir+runnumber+filter+comment+'_FindPar.par'):
		os.remove(datedir+runnumber+filter+comment+'_FindPar.par')
	if os.path.exists(datedir+runnumber+filter+fitsname+'.coo'):
		os.remove(datedir+runnumber+filter+fitsname+'.coo')

	print fwhm_arr, FWHM

	GAIN      = 5.0           #e-/ADU
	READNOISE = 30.0#0.45          #Readnoise for a single
	NUMREADS  = 250              #Number of reads in the ramp
	NUMFILES  = 10              #Number of files used to make dither 
	GAIN      = NUMFILES*GAIN         #If N images are averaged, GAIN=N*GAIN
	READNOISE = READNOISE * np.sqrt(NUMREADS)
	print 'Using FWHM of ' + str(FWHM)
	SKY = np.median(image)
	SIGMA = 3.0#np.std(image)
	MAX = 30000

#FINDPARS ##############################################################
	DataParsFileName =  datedir+runnumber+filter+comment+'_DataPar.par'
#    if OverWrite == 1 or not os.path.exists(DataParsFileName) :
#	iraf.datapars.setParam('threshold', SIGMA)
	iraf.datapars.setParam('scale',1.0)#0.453)
	iraf.datapars.setParam('fwhmpsf',FWHM/2.0)   #FWHM of biggest star
	iraf.datapars.setParam('sigma',SIGMA)   #Std. Dev of sky pixels
#	iraf.datapars.setParam('datamin',SKY-np.sqrt(SKY))  #Min good data value ~Sky-3Sigma
	iraf.datapars.setParam('datamin',150)  #Min good data value ~Sky-3Sigma
	iraf.datapars.setParam('noise','poisson')
	iraf.datapars.setParam('datamax',MAX-100)  #Max good data value
	iraf.datapars.setParam('epadu',GAIN)     #e- per ADU; calculated above
	iraf.datapars.setParam('readnoise',READNOISE)#readnoise
#KEYWORDS IN FITSHEADER
	iraf.datapars.setParam('exposure','EXPOS')    #Integration time
	iraf.datapars.setParam('airmass', 'AIRMASS')    #Airmass at middle
	iraf.datapars.setParam('filter', 'FILTER')    #Filter used
	iraf.datapars.setParam('obstime','TIME_UTC')    #UT at middle
	iraf.datapars.saveParList(filename=DataParsFileName)
	iraf.datapars.setParList(ParList=DataParsFileName)

#FINDPARS ##############################################################
	FindParsFileName =  datedir+runnumber+filter+comment+'_FindPar.par'
    #if OverWrite == 1 or not os.path.exists(FindParsFileName) :
	iraf.findpars.setParam('threshold',SIGMA)
	iraf.findpars.setParam('nsigma', 1.0)
	iraf.findpars.setParam('ratio', 1.0)
	iraf.findpars.setParam('theta', 0.0)
	iraf.findpars.setParam('sharplo', 0.2)
	iraf.findpars.setParam('sharphi',  1.0)
	iraf.findpars.setParam('roundhi', 1.0)
	iraf.findpars.setParam('roundlo', -1.0)
	iraf.findpars.saveParList(filename=FindParsFileName)
	iraf.findpars.setParList(ParList = FindParsFileName)


    #Using the newly calculated FWHM to detect the stars on the image
#######################This is for SExtractor###############################
#    raw_stars = pys.run(image,params=['X_IMAGE','Y_IMAGE'],conf_args={'ANALYSIS_THRESH':1,'FILTER':'Y','FILTER_NAME':'default.conv','DETECT_TYPE':'CCD','DEBLEND_NTHRESH':32,'DEBLEND_MINCONT':0.005,'CLEAN':'Y','CLEAN_PARAM':1.0,'MASK_TYPE':'CORRECT','SEEING_FWHM':'1.2','STARNNW_NAME':'default.nnw','BACK_SIZE':32,'BACK_FILTERSIZE':3,'GAIN':5.0,'DETECT_THRESH':10,'PIXEL_SCALE':0.453,'THRESH_TYPE':'RELATIVE','DETECT_MINAREA':FWHM,'INTERP_MAXXLAG':16,'INTERP_MAXYLAG':16,'INTERP_TYPE':'ALL'})

#########################This is for DAOFIND################################

	iraf.daophot.daofind( datedir+runnumber+filter+fitsname, output = datedir+runnumber+filter+fitsname+'.coo', starmap = datedir+runnumber+filter+fitsname+'.starmap.', skymap = datedir+runnumber+filter+fitsname+'.skymap.', boundary = 'nearest', constant = 0.0, interactive = 'no', wcsout = ')_.wcsout', cache = 'yes', verify = 'no', update = 'yes', verbose = 'yes' )

	raw_stars = np.loadtxt(fname=datedir+runnumber+filter+fitsname+'.coo', dtype = float, usecols = (0,1))

	print raw_stars

	starsx=raw_stars[:,0]
	starsy=raw_stars[:,1]
    
	badstarsx = []
	badstarsy = []
	i=0
	brdr_limitx=det_area*4.0
	brdr_limity=det_area*4.0
	while ( i < len(starsx) ) or ( i < len(starsy) ):
#        print 'Source number '+str(i)+' at location '+str(starsx[i])+' '+str(starsy[i])
		scan = np.copy(image)
		if ( image[starsy[i]-1][starsx[i]-1] < np.sqrt(agtergrond) ) :
#            print 'Star at '+str(starsx[i])+' '+str(starsy[i])+' is dimmer than the noise '+str(agtergrond)
			badstarsx.append(starsx[i])
			badstarsy.append(starsy[i])
			starsx = np.delete(starsx,i)
			starsy = np.delete(starsy,i)
		else:
			r = 1
			s = 1
                    
		i=i+1
	del i
#    print str(len(starsx))+' good sources'
	stars=[starsx,starsy]
	print str(np.size(stars[0])) + ' stars detected'
    
	ax4.set_title('Stars found with DAOFIND') #Plotting the number of stars that were detected with DAOFIND, using the FWHM from the candidates
	plt.xlabel('X')
	plt.ylabel('Y')
	plt.imshow(plot_image, cmap = plt.cm.binary)
	plt.scatter(starsx, starsy)
	plt.scatter(badstarsx,badstarsy,marker='+',c='r')
	plt.draw()
	plt.show()
	
	return stars, image, agtergrond-20, candidates, counter
    
##############################DO PIPELINE FOR ALL OF THE IMAGES#############

def allimages():
    for date in dates:
        for fltr in fltrs:
            lstloc = '../../nir/reduced/'+date+'/'
            with open(lstloc+fltr+'list', "r") as f:
                lines = f.readlines()
                for l in lines:
                    stars,image = stardetection(lstloc+l)
                    plot(stars, image)
                    plt.cla()

#########################PERFORM AUTOMATED PSF PHOTOMETRY####################
def photometry(stars, image,agtergrond,brt_src):

#Check if previous parameter files exists and if they should be removed
	if (OverWrite == 1):
		os.remove(datedir+runnumber+filter+comment+'_daofind.par')
		os.remove(datedir+runnumber+filter+comment+'_CenterPar.par')
		os.remove(datedir+runnumber+filter+comment+'_FitSkyPar.par')
		os.remove(datedir+runnumber+filter+comment+'_FindPar.par')
#        os.remove(datedir+runnumber+filter+comment+'_DataPar.par')
		os.remove(datedir+runnumber+filter+comment+'_DAOPar.par')

#KEY PARAMETERS Paramters to use
#    FWHM      = 15		  #Full Width at Half Max 
	GAIN      = 5.0 		  #e-/ADU
	READNOISE = 30#0.45		  #Readnoise for a single
	NUMREADS  = 250			  #Number of reads in the ramp
	NUMFILES  = 250 		  #Number of files used to make dither 
	GAIN      = NUMFILES*GAIN         #If N images are averaged, GAIN=N*GAIN
	READNOISE = READNOISE * np.sqrt(NUMREADS)
	print 'Using FWHM of ' + str(FWHM)
	SKY = np.median(image)
	SIGMA = 4.0#np.std(image)
	MAX = 30000

#PHOTPARS #############################################################
#    PhotParsFileName = datedir+runnumber+filter+comment+'_PhotPar.par'
#    if OverWrite == 1 or not os.path.exists(PhotParsFileName) :
#        iraf.photpars.setParam('image',datedir+runnumber+PhotParsFileName)
#        iraf.photpars.setParam('coords',datedir+runnumber+'stars.coo')
#        iraf.photpars.setParam('apertures',FWHM)
#        iraf.photpars.saveParList(filename=PhotParsFileName)
#    else :
#        iraf.photpars.setParList(ParList = PhotParsFileName)
            
#CENTERPARS #############################################################
	CenterParsFileName = datedir+runnumber+filter+comment+'_CenterPar.par'
    #if OverWrite == 1 or not os.path.exists(CenterParsFileName) :
	iraf.centerpars.setParam('calgorithm','centroid') #Recommended for first pass
	iraf.centerpars.setParam('cbox',0.25*FWHM)     #Centroiding box
	iraf.centerpars.setParam('cthreshold','0')
	iraf.centerpars.saveParList(filename=CenterParsFileName)
	iraf.centerpars.setParList(ParList=CenterParsFileName)
    #else :
        #iraf.centerpars.setParList(ParList=CenterParsFileName)

#SKYPARS ################################################################
	FitSkyParsFileName = datedir+runnumber+filter+comment+'_FitSkyPar.par'
#    if OverWrite == 1 or not os.path.exists(FitSkyParsFileName) :
	iraf.fitskypars.setParam('salgorithm','median')
	iraf.fitskypars.setParam('skyvalue',SKY) #Average value of the background sky
	iraf.fitskypars.setParam('annulus',FWHM)#+3.0) #The inner radius of the sky anulus
	iraf.fitskypars.setParam('dannulus',4.0) #Width of the sky fitting annulus
	iraf.fitskypar.setParam('smaxiter',15) #Maximum number of sky fitting iterations
	iraf.fitskypar.setParam('snreject',50) #Maximum number of sky fitting rejection iterations
	iraf.fitskypar.setParam('sloreject',3.0) #Lower K-sigma rejection limit in sky sigma
	iraf.fitskypar.setParam('shireject',3.0) #Higher K-sigma rejection limit in sky sigma
	iraf.fitskypar.setParam('khist',3.0) #Half width of histogram in sky sigma
	iraf.fitskypar.setParam('binsize',0.1) #Binsize of histogram in sky sigma
	iraf.fitskypar.setParam('smooth','Yes') #Use boxcar smoothing on sky histogram
	iraf.fitskypar.setParam('rgrow',0.5) #Region growing radius in scale units
	iraf.fitskypars.saveParList(filename=FitSkyParsFileName)
	iraf.fitskypars.setParList(ParList=FitSkyParsFileName)
#    else:
#        iraf.fitskypars.setParList(ParList=FitSkyParsFileName)

    #else :
        #iraf.findpars.setParList(ParList = FindParsFileName)

#DATAPARS - http://stsdas.stsci.edu/cgi-bin/gethelp.cgi?datapars.hlp
	DataParsFileName =  datedir+runnumber+filter+comment+'_DataPar.par'
#    if OverWrite == 1 or not os.path.exists(DataParsFileName) :
	iraf.datapars.setParam('scale',1.0)#0.453)
	iraf.datapars.setParam('fwhmpsf',FWHM*2.0)   #FWHM of biggest star
	iraf.datapars.setParam('sigma',SIGMA)   #Std. Dev of sky pixels
	iraf.datapars.setParam('datamin',SKY-np.sqrt(SKY))  #Min good data value ~Sky-3Sigma
#	iraf.datapars.setParam('datamin',150)  #Min good data value ~Sky-3Sigma
	iraf.datapars.setParam('noise','poisson')
	iraf.datapars.setParam('datamax',MAX-100)  #Max good data value
	iraf.datapars.setParam('epadu',GAIN)     #e- per ADU; calculated above
	iraf.datapars.setParam('readnoise',READNOISE)#readnoise
#KEYWORDS IN FITSHEADER
	iraf.datapars.setParam('exposure','EXPOS')	#Integration time
	iraf.datapars.setParam('airmass', 'AIRMASS')	#Airmass at middle
	iraf.datapars.setParam('filter', 'FILTER')	#Filter used
	iraf.datapars.setParam('obstime','TIME_UTC')	#UT at middle
	iraf.datapars.saveParList(filename=DataParsFileName)
	iraf.datapars.setParList(ParList=DataParsFileName)
 #   else:
#        iraf.datapars.setParList(ParList=DataParsFileName)

#DAOPARS ###############################################################
	DAOParsFileName =  datedir+runnumber+filter+comment+'_DAOPar.par'
#    if OverWrite == 1 or not os.path.exists(DAOParsFileName) :
	iraf.daopars.setParam('psfrad',FWHM)
	iraf.daopars.setParam('fitrad',0.8*FWHM)
	iraf.daopars.setParam('function', 'lorentz')
	iraf.daopars.setParam('varorder',1)	  #Try using variable PSF
	iraf.daopars.setParam('nclean',50)
	iraf.daopars.setParam('sannulus',FWHM+1.5)
	iraf.daopars.setParam('wsannulus',5.0)
	iraf.daopars.setParam('fitsky','yes')    #Refit the sky
	iraf.daopars.setParam('recenter','yes')
#    iraf.daopars.setParam('fitrad',FWHM)
	iraf.daopars.setParam('groupsky','yes')
	iraf.daopars.setParam('flaterr',1.0)
	iraf.daopars.setParam('proferr',10.0)
	iraf.daopars.setParam('maxiter',35)
	iraf.daopars.setParam('clipexp',6.0)
	iraf.daopars.setParam('cliprange',2.5)
	iraf.daopars.setParam('critsnratio',1.0)
     
	iraf.daopars.saveParList(filename=DAOParsFileName)
	iraf.daopars.setParList(ParList = DAOParsFileName)
#    else: 
#        iraf.daopars.setParList(ParList = DAOParsFileName)

########################################################################
########################################################DAOFIND
#Set up daofind to go without prompting for input
	iraf.daofind.setParam('image',image)		#Set ImageName
	iraf.daofind.setParam('verify','no')			#Don't verify
	iraf.daofind.setParam('interactive','no')		#Interactive
	iraf.daofind.saveParList(datedir+runnumber+filter+comment+'_daofind.par')
	iraf.daofind.setParList(datedir+runnumber+filter+comment+'_daofind.par')

    #PHOT PARAMETERS
#    PHOTParsFileName =  datedir+runnumber+filter+comment+'_PHOTpar.par'
#    if OverWrite == 1 or not os.path.exists(PHOTParsFileName) :
#        iraf.photpars.setParam('image',datedir+runnumber+fitsname)
#        irar.photpars.setParam('coords',datedit+runnumber+'stars.coo')
#        iraf.daopars.saveParList(filename=PHOTParsFileName)
#    else: 
#        iraf.daopars.setParList(ParList = PHOTParsFileName)
#    iraf.daophot.phot.setParam('datapars',DataParsFileName)
#    iraf.daophot.phot.setParam('centerpars',CenterParsFileName)
#    iraf.daophot.phot.setParam('fitskypars',)

    #Do the photometric task on the image
	print 'Performing photometry on the detected sources:'
	iraf.daophot.phot(datedir+runnumber+filter+fitsname,coords=datedir+runnumber+filter+fitsname+'.coo',output=datedir+runnumber+filter+fitsname+'.mag', interactive = 'no',verify='no')
	print 'Photometry was performed on '+str(np.shape(stars)[1])+' sources'
	print 'Performing photometry on the brightest candidates:'
#    if os.path.exists(datedir+runnumber+filter+fitsname+'.mag'):
#        os.remove(datedir+runnumber+filter+fitsname+'.mag')
	iraf.daophot.phot(datedir+runnumber+filter+fitsname,coords=datedir+runnumber+filter+fitsname+'.brt_src.coo',output=datedir+runnumber+filter+fitsname+'.brt_src.mag', interactive = 'no',verify='no')
#    iraf.daophot.phot(datedir+runnumber+filter+fitsname,coords=datedir+runnumber+filter+fitsname+'.coo',output=datedir+runnumber+filter+fitsname+'.mag')
	print 'Selecting the stars for the PSF model with the pstselect module'
	if os.path.exists(datedir+runnumber+filter+fitsname+'.pst'):
		os.remove(datedir+runnumber+filter+fitsname+'.pst')
#    if  filter == 'k':
#        iraf.daophot.pstselect(datedir+runnumber+filter+fitsname,photfile=datedir+runnumber+filter+fitsname+'.brt_src.mag',pstfile=datedir+runnumber+filter+fitsname+'.pst',verbose='yes', maxnpsf=12)
#    else:
	iraf.daophot.pstselect(datedir+runnumber+filter+fitsname,photfile=datedir+runnumber+filter+fitsname+'.brt_src.mag',pstfile=datedir+runnumber+filter+fitsname+'.pst',verbose='yes', maxnpsf=canstars - 5, interactive = 'no', verify='no')
#    iraf.daophot.pstselect(datedir+runnumber+filter+fitsname,photfile=datedir+runnumber+filter+fitsname+'.mag',pstfile=datedir+runnumber+filter+fitsname+'.pst',verbose=py'yes', maxnpsf=60)
	print 'Plotting the stars that were selected for the PSF model'
	pstdata = np.loadtxt(datedir+runnumber+filter+fitsname+'.pst')
#    plt.scatter(pstdata[:,1], pstdata[:,2], c='r', marker='d', linewidths=2.0)
	print pstdata
	print 'Fitting the PSF on the brightest sources detected on the image:'
	if os.path.exists(datedir+runnumber+filter+fitsname+'.psf.fits'):
		os.remove(datedir+runnumber+filter+fitsname+'.psf.fits')
	if os.path.exists(datedir+runnumber+filter+fitsname+'.psg'):
		os.remove(datedir+runnumber+filter+fitsname+'.psg')
	if os.path.exists(datedir+runnumber+filter+fitsname+'.opst'):
		os.remove(datedir+runnumber+filter+fitsname+'.opst')
	if os.path.exists(datedir+runnumber+filter+fitsname+'.nst'):
		os.remove(datedir+runnumber+filter+fitsname+'.nst')
	if os.path.exists(datedir+runnumber+filter+fitsname+'.brt_src.nrj'):
		os.remove(datedir+runnumber+filter+fitsname+'.brt_src.nrj')
	if os.path.exists(datedir+runnumber+filter+fitsname+'.sub.fits'):
		os.remove(datedir+runnumber+filter+fitsname+'.sub.fits')
	if os.path.exists(datedir+runnumber+filter+fitsname+'.nrj'):
		os.remove(datedir+runnumber+filter+fitsname+'.nrj')
	if os.path.exists(datedir+runnumber+filter+fitsname+'.als'):
		os.remove(datedir+runnumber+filter+fitsname+'.als')
	if os.path.exists(datedir+runnumber+filter+fitsname+'.als.sub.fits'):
		os.remove(datedir+runnumber+filter+fitsname+'.als.sub.fits')
	iraf.daophot.psf(datedir+runnumber+filter+fitsname, photfile=datedir+runnumber+filter+fitsname+'.brt_src.mag', pstfile=datedir+runnumber+filter+fitsname+'.pst', psfimage=datedir+runnumber+filter+fitsname+'.psf.fits', opstfile=datedir+runnumber+filter+fitsname+'.opst', groupfile=datedir+runnumber+filter+fitsname+'.psg', fwhmpsf=4.0*FWHM,update='no', matchbyid='yes', interactive='no', showplots = 'no', verify='no')
	#    iraf.daophot.psf(datedir+runnumber+filter+fitsname, photfile=datedir+runnumber+filter+fitsname+'.mag', pstfile=datedir+runnumber+filter+fitsname+'.pst', psfimage=datedir+runnumber+filter+fitsname+'.psf.fits', opstfile=datedir+runnumber+filter+fitsname+'.opst', groupfile=datedir+runnumber+filter+fitsname+'.psg', fwhmpsf=FWHM)
	print 'Recalculating the PSF according to the nstar module'
	iraf.daophot.nstar(datedir+runnumber+filter+fitsname, groupfile=datedir+runnumber+filter+fitsname+'.psg', psfimage=datedir+runnumber+filter+fitsname+'.psf.fits', nstarfile=datedir+runnumber+filter+fitsname+'.nst',rejfile=datedir+runnumber+filter+fitsname+'.brt_src.nrj',update='yes', cache = 'no', verify = 'no')
	print 'Subtracting the PSF stars that were calculated with nstar'
	iraf.daophot.substar(datedir+runnumber+filter+fitsname,photfile=datedir+runnumber+filter+fitsname+'.nst',psfimage=datedir+runnumber+filter+fitsname+'.psf.fits',subimage=datedir+runnumber+filter+fitsname+'.sub.fits',cache='no',update='yes',verbose='yes', verify = 'no')#,exfile=datedir+runnumber+filter+fitsname+'.opst'
	print 'Using the calculated PSF on all of the detected stars on the science frame'
	iraf.daophot.allstar(datedir+runnumber+filter+fitsname,photfile=datedir+runnumber+filter+fitsname+'.mag',psfimage=datedir+runnumber+filter+fitsname+'.psf.fits',allstarfile=datedir+runnumber+filter+fitsname+'.als',rejfile=datedir+runnumber+filter+fitsname+'.nrj',subimage=datedir+runnumber+filter+fitsname+'.als.sub.fits',cache='no',verify='no',update='yes',verbose='yes')

	if os.path.exists(datedir+runnumber+filter+fitsname+'.fnt_src.coo'):
		os.remove(datedir+runnumber+filter+fitsname+'.fnt_src.coo')
	if os.path.exists(datedir+runnumber+filter+fitsname+'.fnt_src.mag'):
		os.remove(datedir+runnumber+filter+fitsname+'.fnt_src.mag')
	if os.path.exists(datedir+runnumber+filter+fitsname+'.allsrcs.mag'):
		os.remove(datedir+runnumber+filter+fitsname+'.allsrcs.mag')
	if os.path.exists(datedir+runnumber+filter+fitsname+'.allsrcs.als'):
		os.remove(datedir+runnumber+filter+fitsname+'.allsrcs.als')
	if os.path.exists(datedir+runnumber+filter+fitsname+'.allsrcs.nrj'):
		os.remove(datedir+runnumber+filter+fitsname+'.allsrcs.nrj')
	if os.path.exists(datedir+runnumber+filter+fitsname+'.allsrcs.als.sub.fits'):
		os.remove(datedir+runnumber+filter+fitsname+'.allsrcs.als.sub.fits')
    
#Detecting the stars on the subtracted image, detecting the sources that were not subtracted
	print 'Determining the number of sources that were missed by the pipeline before subtraction'
	iraf.findpars.setParam('threshold', SIGMA - 3.0)
	if agtergrond < 2000:
		agtergrond = 2000
#	iraf.datapars.setParam('datamin', agtergrond - 1900)
	iraf.datapars.setParam('datamin', 150)
	iraf.datapars.setParam('fwhmpsf',1.2*FWHM)   #FWHM of biggest star
	iraf.daophot.daofind(datedir+runnumber+filter+fitsname+'.als.sub.fits', output = datedir+runnumber+filter+fitsname+'.fnt_src.coo')#, threshold = SIGMA -3.0)
#Photometry on the newly detected faint sources from the original image
	print 'Performing aperture photometry on the newly discovered faint sources'
	iraf.daophot.phot(datedir+runnumber+filter+fitsname, coords=datedir+runnumber+filter+fitsname+'.fnt_src.coo',output=datedir+runnumber+filter+fitsname+'.fnt_src.mag', interactive = 'no',verify='no')	
#Appending the newly detected stars' photometry file to the old one.
	conlist = open(datedir+runnumber+filter+fitsname+'.allsrcs.mag','w') #Declare the name for the list of the sources of teh first DAOFIND run and the faint sources that needs to be concatenated/added
	src_file = open(datedir+runnumber+filter+fitsname+'.mag', 'r') #Open the first round of DAOFIND's run
	fntsrc_file = open(datedir+runnumber+filter+fitsname+'.fnt_src.mag', 'r') #The fainter sources that were detected only after the first lost were subtracted
	for line in src_file:
		conlist.write(line)
	del line
	for line in fntsrc_file:
		if line[0] != '#':
			conlist.write(line)
#	conlist.write()
	conlist.close()
#	os.chdir(datedir+runnumber)
#	iraf.concatenate.unlearn()
#	iraf.concatenate(infiles=datedir+runnumber+filter+fitsname+'.mag, '+datedir+runnumber+filter+fitsname+'.fnt_src.mag', outfile=datedir+runnumber+filter+fitsname+'.allsrcs.mag',append='yes')#Concatenate the magnitude files into one that will be used with an allstar routine
#	os.chdir('../../../code/python/')
#Performing photometry on the combined 'stars' file
	iraf.daophot.allstar(datedir+runnumber+filter+fitsname,photfile=datedir+runnumber+filter+fitsname+'.allsrcs.mag',psfimage=datedir+runnumber+filter+fitsname+'.psf.fits',allstarfile=datedir+runnumber+filter+fitsname+'.allsrcs.als',rejfile=datedir+runnumber+filter+fitsname+'.allsrcs.nrj',subimage=datedir+runnumber+filter+fitsname+'.allsrcs.als.sub.fits',cache='no',verify='no',update='yes',verbose='yes')
	
    
##########################Determining the appreture correction###############
def aper_cor(fname):
	GAIN      = 5.0 		  #e-/ADU
	READNOISE = 30#0.45		  #Readnoise for a single
	NUMREADS  = 250			  #Number of reads in the ramp
	NUMFILES  = 250 		  #Number of files used to make dither 
	GAIN      = NUMFILES*GAIN         #If N images are averaged, GAIN=N*GAIN
	READNOISE = READNOISE * np.sqrt(NUMREADS)

#Determining the aperture correction from the PSG and NST files:
	skip_arr = []
	if os.path.exists(fname+'.fnt_cnd.sub'):
		os.remove(fname+'.fnt_cnd.sub')
	if os.path.exists(fname+'.fnt_cnd.sub.fits'):
		os.remove(fname+'.fnt_cnd.sub.fits')
	if os.path.exists(fname+'.ap'):
		os.remove(fname+'.ap')
	if os.path.exists(fname+'.res'):
		os.remove(fname+'.res')
	iraf.ptools.txdump(fname+'.nst', fields='id,group,mag', expr='yes')
	iraf.ptools.pselect(infiles=fname+'.nst', outfiles=fname+'.fnt_cnd.sub', exp='MAG<19.5')
	iraf.daophot.substar(fname, photfile=fname+'.fnt_cnd.sub', psfimage=fname+'.psf.fits', subimage=fname+'.fnt_cnd.sub.fits', gain = GAIN, readnoise = READNOISE)
	iraf.ptools.txdump(fname+'.nst', fields='xc,yc', Stdout=fname+'.ap', expr='MAG<19.5')
	iraf.daophot.phot(fname+'.fnt_cnd.sub.fits',coords=fname+'.ap', output = fname+'.res', aper=str(FWHM)+',18.', annulus=20., dannu=5., verbose='Yes', gain = GAIN, readnoise = READNOISE)
	#iraf.ptools.txselect(fname+'.nst', outfiles=fname+'.tnst', fields='id,group,mag')
#	lines = 0
#	for line in open(fname+'.nst'):
#		lines +=1
#		if( lines > 42 and (np.mod(lines, 2) == 0 ) ):
#			skip_arr.append(lines)
#	print skip_arr
#	with open(fname+'.nst', 'r'):
#		lines = (line for line in f if predicate(line))
#		brt_arr = np.genfromtxt(skip_arr)
#	brt_grp = np.loadtxt(fname+'.nst',usecols=[0,1,4], skiprows=skip_arr)

###########################Standard star photometry#########################

def standard_star(filename, standard):
	hdulist = pyfits.open(filename) #Load the header from the fits file
	std_image = pyfits.getdata(filename)
	wcs = pywcs.WCS(hdulist[0].header) #Create a wcs object from the header keywords
	print wcs.wcs.name #Prints the name of the wcs as given in the header
	wcs.wcs.print_contents() #Prints all of the settings that were parsed from the header
	f = open(filename+'.coo', 'w')
	print standard, standard_ra[standard], standard_dec[standard], wcs.wcs_sky2pix([[standard_ra[standard],standard_dec[standard]]],2)
	f.write(str(wcs.wcs_sky2pix([[standard_ra[standard],standard_dec[standard]]],1)[0][0])+'	'+str(wcs.wcs_sky2pix([[standard_ra[standard],standard_dec[standard]]],1)[0][1])) #For some nights' observations the x and y coordinates needs to be shifted by 32 pixels because the fix.py script timmed some images and not others.
	f.close()	

	GAIN      = 5.0 		  #e-/ADU
	READNOISE = 30.0#0.45		  #Readnoise for a single
	NUMREADS  = 250			  #Number of reads in the ramp
	NUMFILES  = 250 		  #Number of files used to make dither 
	GAIN      = NUMFILES*GAIN         #If N images are averaged, GAIN=N*GAIN
	READNOISE = READNOISE * np.sqrt(NUMREADS)
	SKY = np.median(std_image)

#CENTERPARS ################################################################
	iraf.centerpars.setParam('calgorithm','none') #Recommended for first pass
	iraf.centerpars.setParam('cbox',4)     #Centroiding box
	iraf.centerpars.setParam('cthreshold','0')

#SKYPARS ###################################################################
	iraf.fitskypars.setParam('salgorithm','median')
	iraf.fitskypars.setParam('skyvalue',SKY) #Average value of the background sky
	iraf.fitskypars.setParam('annulus',15)#+3.0) #The inner radius of the sky anulus
	iraf.fitskypars.setParam('dannulus',5.0) #Width of the sky fitting annulus
	iraf.fitskypar.setParam('smaxiter',15) #Maximum number of sky fitting iterations
	iraf.fitskypar.setParam('snreject',50) #Maximum number of sky fitting rejection iterations
	iraf.fitskypar.setParam('sloreject',3.0) #Lower K-sigma rejection limit in sky sigma
	iraf.fitskypar.setParam('shireject',3.0) #Higher K-sigma rejection limit in sky sigma
	iraf.fitskypar.setParam('khist',3.0) #Half width of histogram in sky sigma
	iraf.fitskypar.setParam('binsize',0.1) #Binsize of histogram in sky sigma
	iraf.fitskypar.setParam('smooth','Yes') #Use boxcar smoothing on sky histogram
	iraf.fitskypar.setParam('rgrow',0.5) #Region growing radius in scale units

    #else :
        #iraf.findpars.setParList(ParList = FindParsFileName)

#DATAPARS - http://stsdas.stsci.edu/cgi-bin/gethelp.cgi?datapars.hlp
	iraf.datapars.setParam('scale',1.0)#0.453)
	iraf.datapars.setParam('datamin',150)  #Min good data value ~Sky-3Sigma
	iraf.datapars.setParam('noise','poisson')
	iraf.datapars.setParam('datamax','INDEF')  #Max good data value
	iraf.datapars.setParam('epadu',GAIN)     #e- per ADU; calculated above
	iraf.datapars.setParam('readnoise',READNOISE)#readnoise
#KEYWORDS IN FITSHEADER
	iraf.datapars.setParam('exposure','EXPOS')	#Integration time
	iraf.datapars.setParam('airmass', 'AIRMASS')	#Airmass at middle
	iraf.datapars.setParam('filter', 'FILTER')	#Filter used
	iraf.datapars.setParam('obstime','TIME_UTC')	#UT at middle

#DAOPARS ###############################################################
	iraf.daopars.setParam('psfrad',12)
	iraf.daopars.setParam('fitrad',12)
	iraf.daopars.setParam('function','lorentz')
	iraf.daopars.setParam('varorder',1)	  #Try using variable PSF
	iraf.daopars.setParam('nclean',50)
	iraf.daopars.setParam('sannulus',15)
	iraf.daopars.setParam('wsannulus',5)
	iraf.daopars.setParam('fitsky','yes')    #Refit the sky
	iraf.daopars.setParam('recenter','yes')
	iraf.daopars.setParam('groupsky','no')
	iraf.daopars.setParam('flaterr',1.0)
	iraf.daopars.setParam('proferr',10.0)
	iraf.daopars.setParam('maxiter',35)
	iraf.daopars.setParam('clipexp',6.0)
	iraf.daopars.setParam('cliprange',2.5)
	iraf.daopars.setParam('critsnratio',1.0)
	iraf.photpars.setParam('apertures',12.0)
     
	if os.path.exists(filename+'.mag'):
		os.remove(filename+'.mag')
	iraf.daophot.phot(filename,coords=filename+'.coo',output=filename+'.mag', interactive = 'no', update ='no', cache= 'no', verify= 'no')
	
	iraf.daopars.unlearn()
	iraf.datapars.unlearn()
	iraf.centerpars.unlearn()
	iraf.fitskypars.unlearn()
	iraf.photpars.unlearn()
	
#	foo = raw_input("fa")

def calibration(location, stnd_strs):
	iraf.photcal() #Initiate the package that is used for calibration
	iraf.photcal(_doprint=0) #Load the calibration package
	iraf.photcal.unlearn()
	if os.path.exists(location+'stnd_out.list'):
		os.remove(location+'stnd_out.list')

	stnd_list = open(location+'stnd_img.list', 'w') #Declaring the list of standard stars for the observation directory
	img_list = open(location+'stnd_list.list','w') #Declaring the file list for the standard stars

	for star in stnd_strs:
#		stnd_list.write(star+' : ' + location + 'jPers'+ star + 'n01.fits' + ' ' + location + 'hPers'+ star + 'n01.fits' + ' ' + location + 'kPers'+ star + 'n01.fits\n') #Writing the files associated with the observations to the standard stars to an observation list
#		img_list.write(location+'jPers'+star+'n01.fits\n')
#		img_list.write(location+'hPers'+star+'n01.fits\n')
#		img_list.write(location+'kPers'+star+'n01.fits\n')
		stnd_list.write(star+' : ' + 'j140124Pers'+ star + 'n01.fits' + ' ' + 'h140124Pers'+ star + 'n01.fits' + ' ' + 'k140124Pers'+ star + 'n01.fits\n') #Writing the files associated with the observations to the standard stars to an observation list
		img_list.write('j140124Pers'+star+'n01.fits.mag\n') #Writing the list of magnitude files to the photometry file
		img_list.write('h140124Pers'+star+'n01.fits.mag\n')
		img_list.write('k140124Pers'+star+'n01.fits.mag\n')
	stnd_list.close()
	img_list.close()
	
	os.chdir(location) #Going to the location of the observations and performing the matching with the created files
	iraf.mknobsfile('@stnd_list.list', idfilters='J,H,Ks', wrap='Yes',allfilters='No',verify='No',verbose='Yes', observations='stnd_out.list', imsets='stnd_img.list')
	os.chdir('../../../code/python/')

	#NOTE: MKCONFIG DOES NOT WORK CORRECTLY THUS IT IS CREATED MANUALLY, AS FOLLOWS:
	if os.path.exists(location+'irsf.cfg'):
		os.remove(location+'irsf.cfg')
	trans_cfg = open(location+'irsf.cfg', 'w')
	trans_cfg.write('# Declare the PERSSON JHK standards catalog variables\n')
	trans_cfg.write('\n')
	trans_cfg.write('catalog\n')
	trans_cfg.write('\n')
	trans_cfg.write('J                2	#The J magnitude\n')
	trans_cfg.write('error(J)         3	#The error on the J magnitude\n')
	trans_cfg.write('H                4	#The H magnitude\n')
	trans_cfg.write('error(H)         5	#The error on the H magnitude\n')
	trans_cfg.write('Ks               6	#The Ks magnitude\n')
	trans_cfg.write('error(Ks)        7	#The error on the Ks magnitude\n')
	trans_cfg.write('\n')
	trans_cfg.write('#Declare the observations file\n')
	trans_cfg.write('\n')
	trans_cfg.write('observations\n')
	trans_cfg.write('TJ               3	#The time of the J filter\n')
	trans_cfg.write('XJ               4	#The airmass in the J filter\n')
	trans_cfg.write('xJ               5	#The x coordinates in the J filter\n')
	trans_cfg.write('yJ               6	#The y coordinates in the J filter\n')
	trans_cfg.write('mJ               7	#The instrumental magnitude in the J filter\n')
	trans_cfg.write('error(mJ)        8	#Error in the intrumental magnitude for J\n')
	trans_cfg.write('\n')
	trans_cfg.write('TH               10	#The airmass in the J filter\n')
	trans_cfg.write('XH               11	#The airmass in the J filter\n')
	trans_cfg.write('xH               12	#The x coordinates in the J filter\n')
	trans_cfg.write('yH               13	#The y coordinates in the J filter\n')
	trans_cfg.write('mH               14	#The instrumental magnitude in the J filter\n')
	trans_cfg.write('error(mH)        15	#Error in the intrumental magnitude for J\n')
	trans_cfg.write('\n')
	trans_cfg.write('TKs              17	#The time of the H filter\n')
	trans_cfg.write('XKs              18	#The airmass in the H filter\n')
	trans_cfg.write('xKs              19	#The x coordinates in the H filter\n')
	trans_cfg.write('yKs              20	#The y coordinates in the H filter\n')
	trans_cfg.write('mKs              21	#The instrumental magnitude in the J filter\n')
	trans_cfg.write('error(mKs)       22	#Error in the intrumental magnitude for J\n')
	trans_cfg.write('\n')
	trans_cfg.write('\n')
	trans_cfg.write('# Sample transformation section for the Persson et al JHK system\n')
	trans_cfg.write('\n')
	trans_cfg.write('transformation\n')
	trans_cfg.write('\n')
	trans_cfg.write('fit   j1=-4.000, j2=5.000, j3=0.500\n')
	trans_cfg.write('JFIT : mJ = J + j1 + j2 * XJ + j3 * (J - Ks)\n')
	trans_cfg.write('\n')
	trans_cfg.write('fit   h1=1.000, h2=0.100, h3=0.500\n')
	trans_cfg.write('HFIT : mH = H + h1 + h2 * XH + h3 * (H - Ks)\n')
	trans_cfg.write('\n')
	trans_cfg.write('fit   k1=1.000, k2=0.050, k3=-0.500\n')
	trans_cfg.write('KsFIT : mKs = Ks + k1 + k2 * XKs + k3 * (H - Ks)\n')
	trans_cfg.close()
#	iraf.mkconfig(config=location+'irsf.cfg', catalog='Persson.dat', observations=location+'stnd_out.list', transform = 'tPersson.dat', catdir='', verify='Yes', edit='Yes', check='Yes',verbose='Yes') #Creating the configure file that states which variables correspnd to which quantities in the transformation equations to the standard system

	#Setting the parameters for the solution of the parameters in the transformation equations
	iraf.fitparams.unlearn()
	iraf.fitparams.setParam('weighting','photometric')
	iraf.fitparams.setParam('addscatter','Yes')
	iraf.fitparams.setParam('maxiter',15)
	iraf.fitparams.setParam('nreject',0)
	iraf.fitparams.setParam('low_reject', 3.0)
	iraf.fitparams.setParam('high_reject', 3.0)
	iraf.fitparams.setParam('grow',0.0)
	iraf.fitparams.setParam('interactive','Yes')
	iraf.fitparams.setParam('logfile', location+'irsf.ans.log')
	iraf.fitparams.setParam('log_unmatched','Yes')
	iraf.fitparams.setParam('log_fit','Yes')
	iraf.fitparams.setParam('log_results','Yes')
	iraf.fitparams.setParam('catdir','')

	iraf.fitparams(observations=location+'stnd_out.list', catalogs='Persson.dat', config='irsf.cfg', parameters=location+'irsf.ans') #Solving the parameters for the transformation equations for the standard star observations of the night.
	
def match_2mass(location, filename):
	twomass = np.loadtxt('../../nir/mosaic/2massRCW34_16arcmin_by_16arcmin.tbl', usecols=[0,1,8,9,12,13,16,17]) #Loads the 2mass calibrated data 
#	filedata = np.loadtxt(location+'jhkmosaic.fits.obs')
	print twomass

	hdulist = pyfits.open(location+filename) #Load the header from the fits file
#	std_image = pyfits.getdata(filename)
	wcs = pywcs.WCS(hdulist[0].header) #Create a wcs object from the header keywords
	print wcs.wcs.name #Prints the name of the wcs as given in the header
	wcs.wcs.print_contents() #Prints all of the settings that were parsed from the header
	
	twomass_pix = np.empty((np.size(twomass[:,0]), 2), dtype=float)
	for star_ind in range(0,np.size(twomass[:,0])):
		twomass_pix[star_ind, 0] = wcs.wcs_sky2pix( twomass[star_ind, 0] , twomass[star_ind, 1] , 0 )[0] #Transforms the RA and DEC coordinates of the 2mass sources onto the pixel coordinates
		twomass_pix[star_ind, 1] = wcs.wcs_sky2pix( twomass[star_ind, 0] , twomass[star_ind, 1] , 0 )[1] 
	
	print twomass_pix
	
	
	
#	f.write(str(wcs.wcs_sky2pix([[standard_ra[standard],standard_dec[standard]]],1)[0][0])+'	'+str(wcs.wcs_sky2pix([[standard_ra[standard],standard_dec[standard]]],1)[0][1])) #For some nights' observations the x and y coordinates needs to be shifted by 32 pixels because the fix.py script timmed some images and not others.

def calibrate_2mass(location, filename):
	iraf.photcal() #Initiate the package that is used for calibration
	iraf.photcal(_doprint=0) #Load the calibration package
	iraf.photcal.unlearn()

	jimage = 'j'+fitsname
	himage = 'h'+fitsname
	kimage = 'k'+fitsname

	jphot = jimage + '.allsrcs.als'
	hphot = himage + '.allsrcs.als'
	kphot = kimage + '.allsrcs.als'

	alsfile = open(location+fitsname+'.als.lst','w')
	alsfile.write(jphot+'\n')
	alsfile.write(hphot+'\n')
	alsfile.write(kphot)
	alsfile.close()

	imgfile = open(location+fitsname+'.imgset.lst','w')
	imgfile.write('RCW34 : '+jimage+' '+himage+' '+kimage)
	imgfile.close()
	
	shftfile = open(location+fitsname+'.shftfile.lst', 'w') #The file showing the pixel shifts between the images
	shftfile.write('0.0 0.0\n')
	shftfile.write('-1.0 5.0\n')
	shftfile.write('-7.0 -1.0')
	shftfile.close()
	
	appcor = open(location+fitsname+'.apcor.lst','w') #The file giving the list for the aperture corrections
	appcor.write(jimage+'.res\n')
	appcor.write(himage+'.res\n')
	appcor.write(kimage+'.res')
	appcor.close()

	if os.path.exists(location+'jhk'+fitsname+'.obs'):
		os.remove(location+'jhk'+fitsname+'.obs')

	os.chdir(location) #To corredtly build the corresponding sources between each filter the program has to go into the directory containing the files
	iraf.photcal.mkobsfile(photfiles='j'+fitsname+'.allsrcs.als,h'+fitsname+'.allsrcs.als,k'+fitsname+'allsrcs.als',idfilters='J,H,Ks', imsets=fitsname+'.imgset.lst',observations='jhk'+fitsname+'.obs', wrap='Yes',allfilters='Yes',verify='no', verbose='Yes', shifts=fitsname+'.shftfile.lst', apercors=fitsname+'.apcor.lst') #Run the routine for the observation files
	os.chdir('../../')

	
def invert_calibration(location, fitsnames):
#	olddir = os.getcwd()
#	os.chdir(location)
	jimage = 'j'+fitsname
	himage = 'h'+fitsname
	kimage = 'k'+fitsname

	jphot = jimage + '.allsrcs.als'
	hphot = himage + '.allsrcs.als'
	kphot = kimage + '.allsrcs.als'

	alsfile = open(location+fitsname+'.als.lst','w')
	alsfile.write(jphot+'\n')
	alsfile.write(hphot+'\n')
	alsfile.write(kphot)
	alsfile.close()

	imgfile = open(location+fitsname+'.imgset.lst','w')
	imgfile.write(fitsname+' : '+jimage+' '+himage+' '+kimage)
	imgfile.close()
	
	shftfile = open(location+fitsname+'.shftfile.lst', 'w') #The file showing the pixel shifts between the images
	shftfile.write('0.0 0.0\n')
	shftfile.write('0.0 0.0\n')
	shftfile.write('0.0 0.0')
	shftfile.close()
	
	appcor = open(location+fitsname+'.apcor.lst','w') #The file giving the list for the aperture corrections
	appcor.write(jimage+'.res\n')
	appcor.write(himage+'.res\n')
	appcor.write(kimage+'.res')
	appcor.close()

	if os.path.exists(location+'jhk'+fitsname+'.obs'):
		os.remove(location+'jhk'+fitsname+'.obs')

	os.chdir(location) #To corredtly build the corresponding sources between each filter the program has to go into the directory containing the files
	iraf.photcal.mkobsfile(photfiles='@'+fitsname+'.als.lst',idfilters='J,H,Ks', imsets=fitsname+'.imgset.lst',observations='jhk'+fitsname+'.obs', wrap='Yes',allfilters='Yes',verify='no', verbose='Yes', shifts=fitsname+'.shftfile.lst', apercors=fitsname+'.apcor.lst') #Run the routine for the observation files
	os.chdir('../../../code/python/')
	
	#Transforming back to apparent magnitudes from the instrumental magnitudes
	if os.path.exists(location+'jhk'+fitsname+'.calib'):
		os.remove(location+'jhk'+fitsname+'.calib')
	iraf.photcal.invertfit(observations=location+'jhk'+fitsname+'.obs', config='irsf.conf', parameters=location+'irsf.ans', calib=location+'jhk'+fitsname+'.calib', catalogs='Persson.dat')#, print='xJ,yJ') #Inverts to apparent magnitudes and also prints the X and Y coordinates to the calibrated file
#	os.chdir(olddir)
    
#####################################MAIN PROGRAM###########################
plt.cla()
#Initiate the IRAF packages that are going to be used
iraf.digiphot()
iraf.daophot()
iraf.apphot()
iraf.ptools()
#Import the good packages from IRAF
iraf.digiphot(_doprint=0)

iraf.daophot(_doprint=0)
iraf.apphot(_doprint=0)
iraf.ptools(_doprint=0)

#PARAMETERS in the various files
datapars   = iraf.datapars.getParList()
daopars    = iraf.daopars.getParList()
findpars   = iraf.findpars.getParList()
centerpars = iraf.centerpars.getParList()
skypars    = iraf.fitskypars.getParList()
photpars   = iraf.photpars.getParList()


#time.sleep(0.5)

#allimages()

#stars,image,agtergrond = stardetection('../../nir/raw/131028/n02/jo0.RCW34n02.fits')
datedir = 'G259_deep_new/images/'
runnumber = ''
fitsname = 'G259_deep.fits'
#filter = 'j'
comment = 'G259_deep.fits'
#FWHM = 0
standard_name1 = '9129'
standard_name2 = '9132'
standard_name3 = '9136'
standard_name4 = '9137'
standard_name5 = '9140'
standard_name6 = '9129'

for filter in fltrs:
	if ( filter == 'j' ):
		plt.figure(0)
		jimage, ( ( ax1, ax2 ), ( ax3 , ax4 ) ) = plt.subplots(2, 2, sharex='col', sharey='row')
		jimage.suptitle('The images for the sources that were detected in the J band')
	if ( filter == 'h' ):
		plt.figure(1)
		himage, ( ( ax1, ax2 ), ( ax3 , ax4 ) ) = plt.subplots(2, 2, sharex='col', sharey='row')
		himage.suptitle('The images for the sources that were detected in the H band')
	if ( filter == 'k' ):
		plt.figure(2)
		kimage, ( ( ax1, ax2 ), ( ax3 , ax4 ) ) = plt.subplots(2, 2, sharex='col', sharey='row')
		kimage.suptitle('The images for the sources that were detected in the K band')
	if os.path.exists(datedir+runnumber+filter+fitsname+'.mag'):
		os.remove(datedir+runnumber+filter+fitsname+'.mag')
	if os.path.exists(datedir+runnumber+filter+fitsname+'.brt_src.mag'):
		os.remove(datedir+runnumber+filter+fitsname+'.brt_src.mag')

	if os.path.exists(datedir+filter+'140110Pers'+standard_name1+'n01.fits.coo'):
		os.remove(datedir+filter+'140110Pers'+standard_name1+'n01.fits.coo')
	if os.path.exists(datedir+filter+'140110Pers'+standard_name1+'n01.fits.mag'):
		os.remove(datedir+filter+'140110Pers'+standard_name1+'n01.fits.mag')
	if os.path.exists(datedir+filter+'140110Pers'+standard_name2+'n01.fits.coo'):
		os.remove(datedir+filter+'140110Pers'+standard_name2+'n01.fits.coo')
	if os.path.exists(datedir+filter+'140110Pers'+standard_name2+'n01.fits.mag'):
		os.remove(datedir+filter+'140110Pers'+standard_name2+'n01.fits.mag')
	if os.path.exists(datedir+filter+'140110Pers'+standard_name3+'n01.fits.coo'):
		os.remove(datedir+filter+'140110Pers'+standard_name3+'n01.fits.coo')
	if os.path.exists(datedir+filter+'140110Pers'+standard_name3+'n01.fits.mag'):
		os.remove(datedir+filter+'140110Pers'+standard_name3+'n01.fits.mag')
	if os.path.exists(datedir+filter+'140110Pers'+standard_name4+'n01.fits.coo'):
		os.remove(datedir+filter+'140110Pers'+standard_name4+'n01.fits.coo')
	if os.path.exists(datedir+filter+'140110Pers'+standard_name4+'n01.fits.mag'):
		os.remove(datedir+filter+'140110Pers'+standard_name4+'n01.fits.mag')
	if os.path.exists(datedir+filter+'140110Pers'+standard_name5+'n01.fits.mag'):
		os.remove(datedir+filter+'140110Pers'+standard_name5+'n01.fits.mag')
	if os.path.exists(datedir+filter+'140110Pers'+standard_name6+'n02.fits.mag'):
		os.remove(datedir+filter+'140110Pers'+standard_name6+'n02.fits.mag')
		
#	print('Doing photomtery on the standard star 1')
#	standard_star(datedir+filter+'140124Pers'+standard_name1+'n01.fits', 0)

#	print('Doing photomtery on the standard star 2')
#	standard_star(datedir+filter+'140124Pers'+standard_name2+'n01.fits', 1)

#	print('Doing photomtery on the standard star 3')
#	standard_star(datedir+filter+'140124Pers'+standard_name3+'n01.fits', 2)

#	print('Doing photomtery on the standard star 4')
#	standard_star(datedir+filter+'140124Pers'+standard_name4+'n01.fits', 3)

#	print('Doing photomtery on the standard star 5')
#	standard_star(datedir+filter+'140124Pers'+standard_name5+'n01.fits', 4)

#	print('Doing photomtery on the standard star 6')
#	standard_star(datedir+filter+'140124Pers'+standard_name6+'n02.fits', 0)

	print('Determining the number of stars on the image:')
	stars, image, agtergrond, brt_src, canstrs = stardetection(datedir+runnumber+filter+fitsname)

	print(brt_src)

	print('Writing the position of the bightest sources to a coordindate file:')
	brt_srcfile = open(datedir+runnumber+filter+fitsname+'.brt_src.coo','w')
	for s in range(0,canstrs-1):
		print(brt_src[s,0],brt_src[s,1])
		brt_srcfile.write(str(brt_src[s,0])+'\t'+str(brt_src[s,1])+'\n')    
	brt_srcfile.close()    

	OverWrite = 0#raw_input('Performing photometry on '+str(np.shape(stars[0]))+' stars. If you want to overwrite the parameter file for a clean run press 1. If you want to use the old files press 0:')

	print()
	photometry(stars,image,agtergrond,brt_src)
	aper_cor(datedir+runnumber+filter+fitsname)

	print('Matching the 2mass sources to their corresponding NIR sources:')
	match_2mass(datedir+runnumber, filter+fitsname)
	
#print('Calibrating the standard star extinction coefficients from the custom standard star catalog PERSSON.dat')
#calibration(datedir+runnumber, [standard_name1, standard_name2, standard_name3, standard_name4, standard_name5, standard_name6])

print('Calibrating the NIR magnitudes and colours to their corresponding 2mass sources:')
calibrate_2mass(datedir+runnumber, fitsname)

#print('Constructing correlation observation file and calibrating the colours:')
#invert_calibration(datedir+runnumber, fitsname)
	    
#	if (filter == 'j') or (filter =='h'):
#		plt.clf()

#                plt.cla()
                
                
#close(fig)                


print("Done!")
