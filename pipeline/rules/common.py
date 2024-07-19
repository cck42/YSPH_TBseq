import pandas as pd
import os

def getrundate(wildcards):
	rundate=wildcards.path.split('runs/')[1].split('_')[0]

def getsamplename(paths):
	sample=os.path.basename(wildcards.path).split('_')[1]
	return sample

def datereformat(wildcards):
	newdate=pd.to_datetime(wildcards.date).strftime('%Y%m%d')
	return newdate

def getr1(wildcards):
	samplepaths=pd.read_table(config['seq_data'],dtype=str,header=None)
	matchedpath=samplepaths.loc[matchedpath[0]==wildcards.sample]
