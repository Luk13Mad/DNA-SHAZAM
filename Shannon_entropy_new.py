from Bio import SeqIO
import numpy as np
import pandas as pd
import argparse
import os
import itertools as itt
import re
import sys
import concurrent.futures as cf


parser = argparse.ArgumentParser(description='generate kmers')
parser.add_argument('-in',"--input", type=str,help='input FASTA-file, full path',required=True)
parser.add_argument('-out',"--output", type=str,help='output directory',default=os.getcwd())
parser.add_argument('-w',"--window", type=int,help='scan window size, default 150',default=150)
parser.add_argument('-k',"--kmer", type=int,help='max length of kmers to be generated e.g. if you want max kmer-size=13 pass -k 14',default=4)
parser.add_argument('-np',"--nproc", type=int,help='number of sub-rpocesses to be started',default=1)

args = vars(parser.parse_args())


    
'''
Script generates shannon entropy of kmers in a sequence

output format 
"window" "length" "kmers" "entropy" "origin"
'''

def generate_possible_kmers(kmer):
    '''
    generates alle possible kmers of specified length
    
    Parameters
    ----------
    kmer: length of kmer for which alle possible variations to be generated

    Returns
    -------
    dict
    '''
    bases=['A','C','G','T']
    possible_kmers={}
    for i in range(3,kmer):
        tmp=[''.join(p) for p in itt.product(bases, repeat=i)]
        possible_kmers.update({i:sorted(tmp)})
    return possible_kmers

def window(fseq,window_size):
    '''
    slide reading window over string
    
    Parameters
    ----------
    fseq: a string
    window_size:int 
                size of sliding window, default 150

    Returns
    -------
    iterator 
    '''
    cnt=0
    #max number of windows for given window_size
    max = int(len(fseq)) - int(window_size) + 1
    if window_size > len(fseq):
        sys.exit("Error: window bigger than string \nError in \"window(fseq, window_size)\" ") 
    while max > cnt:
        for i in range(max):
            yield fseq[i:i+window_size]
            cnt +=1

def get_position(windows,k):
    '''
    returns position of each possible kmer if found
    
    Parameters
    ----------
    seq: a string
    k: length of the kmer

    Returns
    -------
    dict
    '''
    position_raw={}
    for w in possible_kmers[k]:
        tmp=[]
        ite=re.finditer(r"(?={})".format(w),windows)
        for i in ite:
            tmp.append(int(i.start())+1)
        position_raw.update({w:tmp})
    return position_raw

def calc_alpha(position):
    '''
    returns alpha value for each found kmer
    
    Parameters
    ----------
    position: a dict with the positions of each found kmer
    
    Returns
    -------
    dict
    '''
    alpha={}
    for w in position.keys():
        p=position[w]
        if len(p)!=0:
            al=[]
            al=[1/(x-y) for x,y in zip(p[1:],p)]
            al=[1/p[0]]+al
            alpha.update({w:al})
        else:
            alpha.update({w:[]})
    return alpha
    
def calc_beta(alpha):
    '''
    returns beta value for each found kmer
    
    Parameters
    ----------
    alpha: a dict with the alpha values of each found kmer

    Returns
    -------
    dict
    '''
    beta={}
    for w in alpha.keys():
        a=alpha[w]
        if len(a)!=0:
            be=list(np.cumsum(a))
            beta.update({w:be})
        else:
            beta.update({w:[]})
    return beta
    
def calc_q(beta):
    '''
    returns q value for each found kmer
    
    Parameters
    ----------
    beta: a dict with the beta values of each found kmer

    Returns
    -------
    dict
    '''
    q={}
    for w in beta.keys():
        x=np.array(beta[w])
        if len(x)!=0:
            beta_sum=np.sum(x)
            y=x/beta_sum
            q.update({w:list(y)})
        else:
            q.update({w:[]})
    return q
    
def inner(x):
    return (x*np.log2(x))
    
def shannon_entropy(q):
    '''
    returns H value for each found kmer
    
    Parameters
    ----------
    q: a dict with the q values of each found kmer

    Returns
    -------
    dict
    '''
    H={}
    for w in q.keys():
        v=q[w]
        if len(v)!=0:
            t=list(map(inner,v))
            tmp=-np.sum(t)
            if tmp==0:
                tmp=1
            H.update({w:[tmp,len(v)]})
        else:
            H.update({w:[]})
    return H
    
def rescale(H):
    '''
    returns normalized to scale [0,1] H value for each found kmer
    normalization done to range H can take [0,log2(k)] where k is the number of elements of the discrete probability distribution
    
    Parameters
    ----------
    H: a dict with the H value for each found kmer
    q: a dict with the q values of each found kmer

    Returns
    -------
    dict
    '''
    for w in H.keys():
        if H[w]:
            min=0
            max=np.log2(H[w][1])
            if max !=0:
                H[w][0]=(H[w][0]-min)/(max-min)
            else:
                H[w][0]=1
        else:
            H[w]="NA"
    return H

def process(H):
    '''
    gets rid of not needed value
    
    Parameters
    ----------
    H: a dict with the H value for each found kmer
    
    Returns
    -------
    dict
    '''
    for w in H.keys():
        if H[w]!="NA":
            H[w]=H[w][0]
    return H
    
def print_result(df,name):
    '''
    prints results to file
    
    Parameters
    ----------
    df: dataframe with data to print
    name: name of the sequence in question
    '''
    for index,row in df.iterrows():
        for k in range(3,args["kmer"]):
            df_res=pd.DataFrame.from_dict(row[k],columns=['entropy'],orient='index')
            df_res["length"]=np.repeat(len(list(row[k].keys())[0]),len(row[k]))
            df_res["window"]=np.repeat(index,len(row[k]))
            df_res.to_csv("{1}/{0}.csv.gz".format(name,args["output"]),sep="\t",header=False,index=True, compression='gzip',mode="a+")
    print("{0} windows generated for read {1}".format(max(df_res["window"])+1,name),file=sys.stderr)
    
def main(raw):
    read=str(raw.seq).strip()
    name=raw.id
    tmp=[]
    for w in window(read,args["window"]):
        tmp.append(w.strip())
    Windows=pd.Series(tmp).rename(name)
    del tmp
    df=pd.DataFrame({"windows":Windows})
    for k in range(3,args["kmer"]):
        temp1=Windows.apply(get_position,k=k)
        temp2=temp1.apply(calc_alpha)
        del temp1
        temp3=temp2.apply(calc_beta)
        del temp2
        temp4=temp3.apply(calc_q)
        del temp3
        temp5=temp4.apply(shannon_entropy)
        del temp4
        temp6=temp5.apply(rescale)
        del temp5
        temp7=temp6.apply(process)
        del temp6
        df[k]=temp7.values
        del temp7
    print_result(df,name)
    
possible_kmers=generate_possible_kmers(args["kmer"])
generator_reads=SeqIO.parse(args["input"],"fasta")#generator returning reads

with cf.ProcessPoolExecutor(max_workers=args["nproc"]) as executor:#starting multiple sub-process
    for r in executor.map(main,generator_reads):
        pass

print("Done")