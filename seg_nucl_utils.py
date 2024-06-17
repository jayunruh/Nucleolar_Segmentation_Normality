#this file has the functions for nucleolar segmentation
import numpy as np
import tifffile
import scipy.ndimage as ndi
import pandas as pd

#this is the structuring element for flood fill labeling
#full structure means diagonal pixels will get connected
fs=np.ones((3,3))

def poorMansRollingBall(img,rollballrad=100,smstdev=5):
    '''
    we don't have an easy to use rolling ball option in python so subtract a smoothed 
    minimum filtered image instead
    '''
    fimg=img.astype(float)
    return fimg-ndi.minimum_filter(ndi.gaussian_filter(fimg,sigma=smstdev),size=rollballrad*2)

def labelClearEdges(mask,border=2):
    '''
    clear objects within border distance of the edge and label them (return nobjects as well)
    '''
    labels=ndi.label(mask,structure=fs)[0]
    edge=np.concatenate([labels[:border,:].flat,labels[-border:,:].flat,
                         labels[:,:border].flat,labels[:,-border:].flat])
    edgevals=np.unique(edge)
    filtered=labels.copy()
    for i in range(len(edgevals)):
        if(edgevals[i]!=0):
            filtered[labels==edgevals[i]]=0.0
    return ndi.label(filtered>0.5,structure=fs)

def filterObjects(labels,nobj,minsize=4,maxsize=-1,fillholes=True):
    '''
    filter out large and small objects from a labeled image (return nobjects as well)
    '''
    objareas=ndi.sum(labels>0.5,labels,range(1,nobj+1))
    filtered=labels.copy()
    for i in range(len(objareas)):
        if(objareas[i]<minsize):
            filtered[labels==(i+1)]=0.0
        elif(maxsize>0 and objareas[i]>maxsize):
            filtered[labels==(i+1)]=0.0
    filled=ndi.binary_fill_holes(filtered>0.5,structure=fs)
    return ndi.label(filled>0.5,structure=fs)

def segmentNuclei(dapiimg,rollballrad=100,smstdev=5,nucthresh=0.1,minsize=1000,maxsize=4000):
    '''
    segment the nuclei and return preprocessed, unfiltered mask,filtered labels and number of nuclei
    '''
    dapisub=poorMansRollingBall(dapiimg)
    #threshold at fraction of max
    dapimask=dapisub>(nucthresh*dapisub.max())
    #eliminate the edge objects and label them
    dapilabels,nnuclei=labelClearEdges(dapimask)
    #filter out large and small nuclei
    dapilabels,nnuclei=filterObjects(dapilabels,nnuclei,minsize=minsize,maxsize=maxsize)
    return dapisub,dapimask,dapilabels,nnuclei

def segmentNucleoli(nuclimg,thirdimg,dapilabels,nnuclei,
                    rollballrad=15.0,rbsigma=1.0,sigma=0.7,nuclthresh=0.4):
    '''
    segment the nucleoli and return preprocessed, preprocessed third, labeled, and nnuceoli
    note that "third" is a third image (not nuceoli or nuclei) that needs measured
    '''
    #start with background subtraction
    nucleolisub=poorMansRollingBall(nuclimg,rollballrad=15.0,smstdev=1.0)
    thirdsub=poorMansRollingBall(thirdimg,rollballrad=15.0,smstdev=1.0)
    #and a filter
    nucleoli=ndi.gaussian_filter(nucleolisub,sigma=0.7)
    third=ndi.gaussian_filter(thirdsub,sigma=0.7)
    #get the min and max values for each nucleus
    minvals=ndi.minimum(nucleoli,labels=dapilabels,index=range(1,nnuclei+1))
    maxvals=ndi.maximum(nucleoli,labels=dapilabels,index=range(1,nnuclei+1))
    threshvals=minvals+nuclthresh*(maxvals-minvals)
    #make the threshold mask
    threshlevels=dapilabels.copy()
    for i in range(len(threshvals)):
        threshlevels[dapilabels==(i+1)]=threshvals[i]
    #apply the mask in areas where there are nuclei
    nucleolimask=nucleoli>threshlevels
    nucleolimask[threshlevels==0.0]=0
    nucleolilabels,nnucleoli=ndi.label(nucleolimask,structure=fs)
    #finally filter out the small nucleoli
    nucleolilabels,nnucleoli=filterObjects(nucleolilabels,nnucleoli,minsize=4,fillholes=False)
    return nucleoli,third,nucleolilabels,nnucleoli

def measureAll(dapilabels,nnuclei,nucleolilabels,nnucleoli,nucleoli,third):
    '''
    measure all of the nucleoli and associate them with their parent nucleus
    '''
    nuclear_means=ndi.mean(nucleoli,dapilabels,range(1,nnuclei+1))
    nuclear_stds=ndi.standard_deviation(nucleoli,dapilabels,range(1,nnuclei+1))
    #use the ndimage sum function with a boolean image to get the area
    nuclear_areas=ndi.sum(dapilabels>0.5,dapilabels,range(1,nnuclei+1))
    nucleolar_means=ndi.mean(nucleoli,nucleolilabels,range(1,nnucleoli+1))
    nucleolar_stds=ndi.standard_deviation(nucleoli,nucleolilabels,range(1,nnucleoli+1))
    nucleolar_areas=ndi.sum(nucleolilabels>0.5,nucleolilabels,range(1,nnucleoli+1))
    #measure the dapi label "intensity" in each nucleolus to get it's parent label
    nucleolar_ids=ndi.mean(dapilabels,nucleolilabels,range(1,nnucleoli+1))
    #now measure the third image
    third_means=ndi.mean(third,nucleolilabels,range(1,nnucleoli+1))
    third_stds=ndi.standard_deviation(third,nucleolilabels,range(1,nnucleoli+1))
    thirdnuc_means=ndi.mean(third,dapilabels,range(1,nnuclei+1))
    thirdnuc_stds=ndi.standard_deviation(third,dapilabels,range(1,nnuclei+1))
    #count the number of nucleoli per nucleus
    nucleolar_count=[0]*nnuclei
    for i in range(nnucleoli):
        nid=int(np.round(nucleolar_ids[i]))-1
        nucleolar_count[nid]+=1
    #now get the measurements for every nucleolus as a dictionary
    #nuclear measurements: nucleolar_id, nucid, area, avg, stdev, nucleolar count; nucleolar measurements: area, avg, stdev,thirdnucavg,thirdnucstdev,thirdavg,thirdstdev
    measurements=[]
    for i in range(nnucleoli):
        nid=int(np.round(nucleolar_ids[i]))-1
        mdict={'id':i,'nuclear_id':nid+1,'nuclear_area':nuclear_areas[nid],'nuclear_avg':nuclear_means[nid],
               'nuclear_stdev':nuclear_stds[nid],'number_nucleoli':nucleolar_count[nid],'nucleolar_area':nucleolar_areas[i],
               'nucleolar_avg':nucleolar_means[i],'nucleolar_stdev':nucleolar_stds[i],'third_nucavg':thirdnuc_means[nid],
               'third_nucstdev':thirdnuc_stds[nid],'third_nuclavg':third_means[i],'third_nuclstdev':third_stds[i]}
        measurements.append(mdict)
    #return as a dataframe
    return pd.DataFrame(measurements)