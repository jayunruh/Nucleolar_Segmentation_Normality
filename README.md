# Nucleolar_Segmentation_Normality
This repository contains python versions of tools for segmenting nucleoli and measuring their normality as described by Potapoval et. al. (https://elifesciences.org/reviewed-preprints/88799).  Python doesn't reproduce the exact numbers that Fiji does (see https://github.com/jayunruh/Jay_Plugins/blob/master/segment_nucleoli_jru_v4.java) because of slight algorithmic differences but the numbers are statistically very similar and I would expect trends to be the essentially identical.  The dependencies here should be limited to numpy, tifffile, matplotlib, scipy, and pandas which I would recommend installing via conda-forge.  The code takes 3 images: dapi, a nucleolar marker (UBTF in our case) and a nucleolar granular component marker (e.g. Nucleolin).  Normality measures the ratio of the granular component nucleolar enrichment to the UBTF enrichment.  High normality indicates normal enrichment of nucleolin in the nucleolus whereas low normality indicates nucleolar dispersal which is often a response to stress or ribosomal transcriptional disruption.

Here is a listing of the files and their function:

 * segment_nucleoli_demo.ipynb: a walkthrough demonstration of segmenting and measuring a single set of images.
 * seg_nucl_utils.py: the functions from segment_nucleoli_demo in stand alone form.
 * measure_nucleolar_normality.ipynb: demonstrates calculating the normality score from the individual measurements.
 * batch_segment_nucleoli.ipynb: demonstrates high throughtput measurements of all images in example_data.
 * The java_outputs folder has the mask, nucleolar ROIs and measurements for one image from java.