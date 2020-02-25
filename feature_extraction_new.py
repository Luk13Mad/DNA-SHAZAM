import numpy as np
import pandas as pd
import argparse
import os
import sys
import itertools as itt
import concurrent.futures as cf
import time


parser = argparse.ArgumentParser(description='feature extraction')
parser.add_argument('-in',"--input", type=str,help='input file, full path',required=True)
parser.add_argument('-out',"--output", type=str,help='output directory',default=os.getcwd())
parser.add_argument('-w',"--window", type=int,help='number of columns in matrix, default 100',default=100)
parser.add_argument('-s',"--sub", type=int,help='number of rows/columns for submatrices (nxn), default 4',default=4)


args = vars(parser.parse_args())

'''
script extracts feature vector from input. input is file with shannon entropy
'''

def read_input_to_df(file):
    '''
    reads .csv file with entropy data of kmers. .csv-file generated by Shannon_entropy.py
    labels columns 
    replace NA with 0
    
    Parameters
    ----------
    file: input file with data

    Returns
    -------
    dataframe
    '''
    try:
        df=pd.read_csv(file, sep ="\t", header = None, compression='gzip')
    except FileNotFoundError:
        sys.exit("input file not found")

    #replaces NaN with 0
    df.columns=["kmer","entropy","length","window"]
    df["entropy"]=df["entropy"].fillna(0)
    return df

def split_kmer_length(df):
    '''
    returns sub dataframe where only one kmer-length is represented but all windows are present
    alphabetically orders each window by name of kmer
    
    Parameters
    ----------
    df: dataframe

    Returns
    -------
    dataframe
    '''
    m=max(df["length"])
    n=max(df["window"])
    for i in range(3,m+1):
        x=df.loc[df["length"]==i].copy(deep=True)
        yield x
            
def shorten(df,window):
    '''
    shortens dataframe so number of windows is multiple of matrix columns with minimal loss
    
    Parameters
    ----------
    df: dataframe

    Returns
    -------
    dataframe
    '''
    m=len(set(df["window"]))
    overhang=m%window
    new=m-overhang
    df_short=df[df.window<=new-1]
    print("{0} windows cut from df kmer-length = {1}".format(overhang,df["length"].iloc[0]),file=sys.stderr)
    return df_short

def iterate_df(df,window):
    '''
    returns subsets of shortened dataframe sorted by window e.g. window 0 to 99
    
    Parameters
    ----------
    df: dataframe

    Returns
    -------
    dataframe
    '''
    m=len(set(df["window"]))
    
    if (m%window)!=0:
        sys.exit("needs shortened dataframe")
    print("{0} matrices will be generated".format(m/window),file=sys.stderr)
    for i in range(0,m,window):
        df_lower=df[df["window"]>=i]
        df_twoside=df_lower[df_lower["window"]<(i+window)]
        yield df_twoside
    
def generate_matrix(df,length,window):
    '''
    returns matrix with entropy, each column is one time-window
    
    Parameters
    ----------
    df: dataframe
    length:size of kmers
    window:size of window
    
    Returns
    -------
    matrix
    '''
    minimum=min(df["window"])
    maximum=max(df["window"])
    lines=4**int(length)
    matrix=np.zeros((lines,window))
    for cnt,i in enumerate(range(minimum,maximum+1)):
        matrix[:,cnt]=df[df["window"]==i]["entropy"]
    return matrix

def process_matrix(matrix,length,window,sub):
    '''
    processes matrix which contains entropy values for 100 windows
    sets all values which are <mean(matrix) to zero
    divides matrix into 4x4-submatrix and calculates sum for these submatrices
    saves sum to vector
    
    Parameters
    ----------
    matrix: matrix from generate matrix
    length:size of kmers
    window:size of window

    Returns
    -------
    vector
    '''
    mean=matrix.mean()
    tmp1=np.where(matrix < mean,0,matrix)
    vector=[]
    lines=4**int(length)
    for i, j in itt.product(list(range(0,lines,sub)),list(range(0,window,sub))):
        vector.append(np.sum(tmp1[i:i+sub,j:j+sub]))
    return vector

def processing(df):
    '''
    runs process for given input
    
    Parameters
    ----------
    df: dataframe from split_kmer_length

    Returns
    -------
    dict
    '''
    dfs=shorten(df,args["window"])
    check_input(args["window"],dfs["length"].iloc[0],args["sub"])
    intermediate={}
    for cnt,sub_df in enumerate(iterate_df(dfs,args["window"])):
        a1=generate_matrix(sub_df,sub_df["length"].iloc[0],args["window"])
        a2=process_matrix(a1,sub_df["length"].iloc[0],args["window"],args["sub"])
        tmp2=np.argsort(a2)[-6:]
        tmp2=int("".join(map(str,tmp2)))
        intermediate.update({cnt:tmp2})
    df_inter=pd.DataFrame.from_dict(intermediate,columns=['feat'],orient="index")
    return df_inter
    
def check_input(window,length,sub):
    '''
    checks if matrix generation can be done without crash
    
    Parameters
    ----------
    length:size of kmers
    window:size of window

    Returns
    -------
    '''
    lines=4**int(length)
    if (lines%sub)!=0:
        sys.exit("chosen kmer length doesn't work for matrix generation")
    elif (window%sub)!=0:
        sys.exit("chosen window length doesn't work for matrix generation")
    
def print_results(results,out):
    '''
    prints results from processing() to file
    
    Parameters
    ----------
    results: dict with results
    out: output directory

    Returns
    -------
    '''
    printdf=pd.DataFrame()
    for df in results:
        printdf=printdf.append(df)
            
    printdf.to_csv("{0}/{1}.featvec.csv.gz".format(out,args["input"].split("/")[-1]),sep="\t",header=False,index=True, compression='gzip')

#main loop    
result={}
df=read_input_to_df(args["input"])
workers=int(max(df["length"]))
arguments=[]
length=[]
for f in split_kmer_length(df):
    arguments.append(f)
    length.append(f["length"].iloc[0])

with cf.ProcessPoolExecutor(max_workers=workers) as executor:
    result=tuple(executor.map(processing,arguments))
    
for r,l in zip(result,length):
    r["length"]=np.repeat(l,len(r))

print_results(result,args["output"])
print("Done")