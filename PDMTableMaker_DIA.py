# Module metadata variables
__author__ = "Cristina Amparo Devesa Arbiol"
__credits__ = ["Cristina Amparo Devesa Arbiol", "Jose Rodriguez", "Samuel Lozano-Juarez", "Jesus Vazquez"]
__license__ = "Creative Commons Attribution-NonCommercial-NoDerivs 4.0 Unported License https://creativecommons.org/licenses/by-nc-nd/4.0/"
__version__ = "0.1.0"
__maintainer__ = "Jose Rodriguez"
__email__ = "cristinaamparo.devesa@cnic.es;jmrodriguezc@cnic.es"
__status__ = "Development"

# Import modules
import re
from Bio import SeqIO
import numpy as np
import configparser
import argparse
import os
import logging
import sys
from collections import OrderedDict  

###################
# Local functions #
###################
def get_fasta_report(file):
    '''
    Create the FASTA report
    '''
    def _create_key_id(rec):
        if (rec.startswith("sp") or rec.startswith("tr") or rec.startswith("cRAP_sp")) and "|" in rec:
            return rec.split("|")[1]
        else:
            return rec
    indb = SeqIO.index(file, "fasta", key_function=_create_key_id)
    return indb


def Obtain_n(MasterProtein, dicc_fasta,clean_seq,m,npos) : 

    MasterProtein = MasterProtein.strip(" ")
    
    final_q_pos = ""
    #The fasta sequence corresponding to this identifier is saved 
    for iden in dicc_fasta:
        if MasterProtein == iden:
            result=str(dicc_fasta[iden].seq.upper()).replace("X","L")
            break
    
    

    
    pattern=re.compile(clean_seq.replace("L","l").replace("I","[IL]").replace("l","[IL]")) # Problems that may exist with leucine and isoleucine are solved
    
    dicc_seqs={}
    pos = 0
    
    # The corresponding fasta sequence is rigorously scrutinized so that no chance is missed  
    listab=[]
    listae = []
    listan = []
    while True:
        match = pattern.search(result, pos)


        if not match:
            break
        s = match.start()
        e = match.end()
        s = s+1
        if s-1 == 0:
            pos1 = 0
        else:
            pos1 = s-2
            
        try:
            p2 = result[e+1]
            pos2 = e+1
        except:
            pos2 = e
        s=s-1
        listab.append(str(s+1))
        listae.append(str(e))
        
        final_pos = s+1
        if final_q_pos == "":
            final_q_pos = str(final_pos+m-1)
            initial_q_pos = final_pos
        else: 
            final_q_pos = final_q_pos+";"+str(final_pos+m-1)
            initial_q_pos = final_pos
            
        
        pos = e #  Move forward in text for the next search


    for position in listab: 
        listan.append(str(int(position)+npos-1))

    b_e_list =  [val for pair in zip(listab, listae) for val in pair] 
    new_b_e_list = []
    for i in range(len(b_e_list)): 
        if float(i)%2 ==0: 
            new_b_e_list.append(b_e_list[i]+"-"+b_e_list[i+1])
    initial_final =";".join(new_b_e_list )
    initial = ";".join(listab)
    final = ";".join(listae)
    ns = ";".join(listan)
    multiple_n = "NO"
    if initial.find(";")!=-1:
        multiple_n = "YES"

    return initial, final, ns, multiple_n,initial_final

def ListMaker(df,seq,counts,dicc_fasta,MasterProtein):
    
    """
    ListMaker returns the frequency of the scan, the amino acid in the first position, cleaned sequence and m and l positions.
    """

    if seq.find("#")!=-1: 
        npos = 0
        clean_seq = seq[seq.find("#")+1:] # Clean sequence is obtained.
        aa = "U" # Aminoacid first position
        freq=int(counts.loc[seq]) # Frecuency of the scan
        m = 0 # Position starting from N-term
        l = -1 # Position starting from C-term
        n = 0 # Protein position
        pd = clean_seq+":"+seq[seq.find("[")+1:seq.find("]")] # Amino acid positions and DM
        dqna = seq[seq.find("[")+1:seq.find("]")]
        b = Obtain_n(MasterProtein,dicc_fasta,clean_seq,m,npos)[0] # Protein positions
        e = Obtain_n(MasterProtein,dicc_fasta,clean_seq,m,npos)[1]
        n = b 
        b_e = Obtain_n(MasterProtein,dicc_fasta,clean_seq,m,npos)[4]
        multiple_n = Obtain_n(MasterProtein,dicc_fasta,clean_seq,m,npos)[3]
        if multiple_n=="YES":
            first_b = b.split(";")[0]
            first_n = first_b
        else:
            first_b = b 
            first_n = b 
                               
    elif seq.find(":")!=-1:
        npos = 0 
        clean_seq = seq[:seq.find(":")] # Clean sequence is obtained.
        aa = "U" # Aminoacid first position
        freq=int(counts.loc[seq]) # Frecuency of the scan
        m = -1 # Position starting from N-term
        l = 0 # Position starting from C-term
        pd = clean_seq+":"+seq[seq.find(":")+1]  # peptide and DM
        n = len(clean_seq)+1 # Protein position
        dqna = seq[seq.find("[")+1:seq.find("]")]
        b = Obtain_n(MasterProtein,dicc_fasta,clean_seq,m,npos)[0] # Protein positions
        e = Obtain_n(MasterProtein,dicc_fasta,clean_seq,m,npos)[1]
        n = b 
        b_e = Obtain_n(MasterProtein,dicc_fasta,clean_seq,m,npos)[4]
        multiple_n = Obtain_n(MasterProtein,dicc_fasta,clean_seq,m,npos)[3]
        if multiple_n=="YES":
            fisrt_b = b.split(";")[0]
            first_n = first_b
        else:
            first_b = b 
            first_n = b 
                               
    elif seq.find("_") != -1:
        clean_seq = seq[:seq.find("_")]
        npos = 0 
        aa = "U" # Aminoacid first position
        freq=int(counts.loc[seq]) # Frecuency of the scan
        m = np.nan # Position starting from N-term
        l = np.nan # Position starting from C-ter
        pd = (seq[:seq.find("_")])+":"+(seq[seq.find("_")+1:]) # Peptide and DM
        dqna = (seq[seq.find("_")+1:])
        b = Obtain_n(MasterProtein,dicc_fasta,clean_seq,m,npos)[0] # Protein positions
        e = Obtain_n(MasterProtein,dicc_fasta,clean_seq,m,npos)[1]
        n = b
        b_e = Obtain_n(MasterProtein,dicc_fasta,clean_seq,m,npos)[4]
        multiple_n = Obtain_n(MasterProtein,dicc_fasta,clean_seq,m,npos)[3]
        if multiple_n=="YES":
            first_b = b.split(";")[0]
            first_n = first_b
        else:
            first_b = b 
            first_n = b 
                                 
    elif seq.find("(") != -1:
        clean_seq = seq[:seq.find("(")]+seq[seq.find(")")+1:] # Clean sequence is obtained.
        npos = seq.find("(")
        aa = seq[seq.find("(")-1] # Aminoacid first position
        freq=int(counts.loc[seq]) # Frecuency of the scan
        m = int(seq.find("(")) # Position starting from N-term
        l = int(len(clean_seq)-m+1) # Position starting from C-term
        pd = clean_seq+":"+seq[seq.find("(")+1:seq.find(")")] # peptide and DM
        dqna = seq[seq.find("(")+1:seq.find(")")] 
        n = Obtain_n(MasterProtein,dicc_fasta,clean_seq,m,npos)[2] # Protein positions
        b = Obtain_n(MasterProtein,dicc_fasta,clean_seq,m,npos)[0] # Protein positions
        e = Obtain_n(MasterProtein,dicc_fasta,clean_seq,m,npos)[1]
        multiple_n = Obtain_n(MasterProtein,dicc_fasta,clean_seq,m,npos)[3]
        b_e = Obtain_n(MasterProtein,dicc_fasta,clean_seq,m,npos)[4]
        if multiple_n=="YES":
            first_n = n.split(";")[0]
            first_b = b.split(";")[0]
        else:
            first_n = n 
            first_b = b 
            
    elif seq.find("[") != -1:
        clean_seq = seq[:seq.find("[")]+seq[seq.find("]")+1:] # Clean sequence is obtained.
        npos = seq.find("[")
        aa = seq[seq.find("[")-1] # Aminoacid first position
        freq=int(counts.loc[seq]) # Frecuency of the scan
        m = int(seq.find("[")) # Position starting from N-term
        l = int(len(clean_seq)-m+1) # Position starting from C-term
        pd = clean_seq+":"+seq[seq.find("[")+1:seq.find("]")] # peptide and DM
        dqna = seq[seq.find("[")+1:seq.find("]")] 
        n = Obtain_n(MasterProtein,dicc_fasta,clean_seq,m,npos)[2] # Protein positions
        b = Obtain_n(MasterProtein,dicc_fasta,clean_seq,m,npos)[0] # Protein positions
        e = Obtain_n(MasterProtein,dicc_fasta,clean_seq,m,npos)[1]
        multiple_n = Obtain_n(MasterProtein,dicc_fasta,clean_seq,m,npos)[3]
        b_e = Obtain_n(MasterProtein,dicc_fasta,clean_seq,m,npos)[4]
        if multiple_n=="YES":
            first_n = n.split(";")[0]
            first_b = b.split(";")[0]
        else:
            first_n = n 
            first_b = b
            
    else:
        clean_seq = seq # Clean sequence is obtained.
        npos = 0
        aa = "U" # Aminoacid first position
        freq=int(counts.loc[seq]) # Frecuency of the scan
        m = np.nan # Position starting from N-term
        l = np.nan # Position starting from C-term
        pd = clean_seq # peptide and DM
        dqna = seq[seq.find("(")+1:seq.find(")")] 
        b = Obtain_n(MasterProtein,dicc_fasta,clean_seq,m,npos)[0] # Protein positions
        e = Obtain_n(MasterProtein,dicc_fasta,clean_seq,m,npos)[1]
        n = b
        multiple_n = Obtain_n(MasterProtein,dicc_fasta,clean_seq,m,npos)[3]
        b_e = Obtain_n(MasterProtein,dicc_fasta,clean_seq,m,npos)[4]
        if multiple_n=="YES":
            first_b = b.split(";")[0]
            first_n = first_b
        else:
            first_n = n 
            first_b = b 
        

    return clean_seq,seq,aa,freq,m,l,pd,n,dqna,b,e,first_b,first_n,b_e

def get_d(val):
    """
    Extract the deltamass value from a peptide sequence
    """
    #search the deltamass between brackets [ ]
    match = re.search(r'\[([^\[\]]+)\]', val)
    if match:
        return match.group(1)
    #search between parentheses ( )
    match = re.search(r'\(([^\(\)]+)\)', val)
    if match:
        return match.group(1)
    #search if the mod is labile _
    if '_' in val:
        return val.split('_', 1)[1]
    #if  there are not previous match, return 0 as deltamass
    return '0'

def get_misscleavages(string):
    """
    Get the number of missed cleavages in a peptide sequence.
    KK, KR, KP, RR, RK and RP are not considered as missed cleavages.
    """
    #keep only capital letters, deleting lower cases and special characters
    seq = ''.join([c for c in string if c.isupper()])
    count = 0
    #iterate through the sequence
    for i, aa in enumerate(seq):
        if aa in ('K', 'R'):
            if i < len(seq) - 1 and seq[i+1] not in ('K', 'R', 'P'):
                count += 1
    return count

def build_pgm(row, seq_column_name = "pdm"):
   """
   Build the pgm sequence from the modified sequence and group column
   """
   seq = re.sub(r'[^A-Z]', '', str(row[seq_column_name]))
   g = row["g"]
   pgm = f"{seq};{g}"
   if g != "NM":
       extra = ""
       seq_val = str(row[seq_column_name])
       if "[" in seq_val:
           extra = str(seq_val.index("["))
       elif "(" in seq_val:
           extra = str(seq_val.index("("))
       elif "_" in seq_val:
           extra = ""
       pgm += f";{extra}"
   return pgm

##################
# Main functions #
##################
def main(file,infile1,fastafile):
    
    """
    Reading configuration file
    """
    import pandas as pd
    config = configparser.ConfigParser(inline_comment_prefixes='#')
    config.read(file)
    logging.info("Reading PDMTableMaker configuration file")  
    
    seq_column_name = config["PDMTableMaker_Parameters"].get("sequence_column_name") # Sequence with DM column name
    MasterProtein_column_name =  config["PDMTableMaker_Parameters"].get("MasterProtein_column_name") # Master protien name column name 
    Outfile_suffix =  config["PDMTableMaker_Parameters"].get("Outfile_suffix") # Chosen suffix for output file

    logging.info("Processing input file")
    
    dicc_fasta = get_fasta_report(fastafile)
    output = infile1
    df = pd.read_csv(infile1, sep="\t", float_precision='high',low_memory=False)
    
    df_withPGM = df.copy()
    df_withPGM["d"] = df_withPGM[seq_column_name].apply(get_d)
    df_withPGM["g"] = df_withPGM["d"].apply(lambda x: x if x != "0" else "NM")
    df_withPGM["pgm"] = df_withPGM.apply(build_pgm, args = (seq_column_name,),axis=1)
    df_withPGM.to_csv(output[:-4]+"_plus.tsv", index=False, sep='\t', encoding='utf-8')

    counts = df[seq_column_name].value_counts() # The number of times that the species appear is saved in the variable count
    df2 = pd.DataFrame(columns=["p","pdm","q","qFreq","pFreq","pd","d","Missing_Cleavage","a","m","n","l","b","e","first_b","first_n","c","k","qdna","qna","qk","qc","qdnaFreq","qnaFreq","qkFreq","qcFreq","A","M","L","N","qdNA","m_left", "m_right"],dtype=float) # Dataframe 2 is created with the aim of 

    cont = 0
    seqlist = [] # In this list it will be saved the sequences already analyzed
    dic_pd_M = {}
    dic_b_e = {}
    dic_b_e_mn = {}
    dic_qna_p = {}
    dic_qna_p_mn = {}
    dic_qdna_freq = {}
    dic_qna_freq = {}
    dic_q_freq = {}
    dic_p_freq = {}
    for index, row in df.iterrows():
        if row[seq_column_name] not in seqlist:

            seqlist.append(row[seq_column_name])
            p,seq,aa,ScanFreq,m,l,pd,n,dqna,b,e,first_b,first_n,b_e= ListMaker(df,row[seq_column_name],counts,dicc_fasta,row[MasterProtein_column_name])              
            d = get_d(row[seq_column_name])
            pdm = row[seq_column_name]
            qId = row[MasterProtein_column_name]

            df2.loc[cont,"p"] = p
            df2.loc[cont,"q"] = qId
            df2.loc[cont,"pdm"] = pdm
            df2.loc[cont,"pd"] = pd
            df2.loc[cont,"d"] = d
            df2.loc[cont,"a"] = aa
            df2.loc[cont,"m"] = m
            df2.loc[cont,"l"] = l
            df2.loc[cont,"n"] = n
            df2.loc[cont,"b"] = b
            df2.loc[cont,"e"] = e
            df2.loc[cont,"first_b"] = first_b
            df2.loc[cont,"first_n"] = first_n
            df2.loc[cont,"Missing_Cleavage"] = get_misscleavages(row[seq_column_name])

            df2.loc[cont,"ScanFreq"] =  int(ScanFreq)

            if row[seq_column_name].find("_")==-1: 
                qdna = qId+":"+str(dqna)+":"+str(n)+":"+aa
                df2.loc[cont,"qdna"] = qdna
                qna = qId+":"+str(n)+":"+aa
                df2.loc[cont,"qna"] = qna
            else:
                qdna = qId+":"+str(dqna)+":"+str(n)+":"+"U"
                df2.loc[cont,"qdna"] = qdna
                qna = qId+":"+str(n)+":"+"U"
                df2.loc[cont,"qna"] = qna

            cont = cont+1
            # The M and L positions with maximun ScanFreq are saved in a dictionary
            if row[seq_column_name].find("_")==-1 and row[seq_column_name].find("#")==-1 and row[seq_column_name].find(":")==-1:
                if p not in dic_pd_M.keys():
                    dic_pd_M[p] = {}
                if d not in dic_pd_M[p].keys():
                    dic_pd_M[p][d] = {}
                    dic_pd_M[p][d] = ScanFreq,m,l,n,qdna,qna
                if ScanFreq > dic_pd_M[p][d][0]:
                    dic_pd_M[p][d] = ScanFreq,m,l,n,qdna,qna

            # Just one b value 
            if  b.find(";")==-1:
                if qId not in dic_b_e.keys() : 
                    dic_b_e[qId]={}

                if p not in dic_b_e[qId].keys() :
                    dic_b_e[qId][p]=int(b),int(e)

            # More than one b value 
            if b.find(";")!=-1:  
                if qId not in dic_b_e_mn.keys(): 
                    dic_b_e_mn[qId]={}

                if p not in dic_b_e_mn[qId].keys() :
                    dic_b_e_mn[qId][p]={}
                if pdm not in dic_b_e_mn[qId][p].keys():
                    dic_b_e_mn[qId][p][pdm]=b_e


            # Just one n value 
            if n.find(";")==-1:
                if qId not in  dic_qna_p.keys() and n.find(";")==-1:
                    dic_qna_p[qId]={}
                if p not in  dic_qna_p[qId].keys():

                    dic_qna_p[qId][p]=""

                    dic_qna_p[qId][p]=dic_qna_p[qId][p]+str(n)
                else:
                    dic_qna_p[qId][p]=dic_qna_p[qId][p]+","+str(n)


            # More than one n value 
            if n.find(";")!=-1:
                if qId not in  dic_qna_p_mn.keys() and n.find(";")!=-1:
                    dic_qna_p_mn[qId]={}
                if p not in  dic_qna_p_mn[qId].keys():

                    dic_qna_p_mn[qId][p]={}

                if pdm not in  dic_qna_p_mn[qId][p].keys():
                    dic_qna_p_mn[qId][p][pdm]=n

                
            

            if qId not in  dic_q_freq.keys():
                dic_q_freq[qId]={}
                dic_q_freq[qId]=int(ScanFreq)

            else:
                dic_q_freq[qId]=  dic_q_freq[qId]+int(ScanFreq)


            if p not in  dic_p_freq.keys():
                dic_p_freq[p]=int(ScanFreq)

            else:
                dic_p_freq[p]=  dic_p_freq[p]+int(ScanFreq)

            if qdna not in  dic_qdna_freq.keys():
                dic_qdna_freq[qdna]=int(ScanFreq)

            else:
                dic_qdna_freq[qdna]=  dic_qdna_freq[qdna]+int(ScanFreq)
            if qna not in  dic_qna_freq.keys():
                dic_qna_freq[qna]=int(ScanFreq)

            else:
                dic_qna_freq[qna]=  dic_qna_freq[qna]+int(ScanFreq)
                    


    dic_qna_p_sort = {}           
    for q in dic_qna_p:
        dic_qna_p_sort[q]= {}

        
        for p in dic_qna_p[q]:
            
            listqna = dic_qna_p[q][p].split(",")
            listqna = list(map(int, listqna))
            listqna = sorted(listqna)

            longitud = len(listqna)
            if longitud ==1: 
                nb= listqna[0]
                ne = listqna[0]
            else: 
                nb = listqna[0]
                ne = listqna[longitud-1]
            dic_qna_p_sort[q][p]= nb,ne


    dic_qna_cluster={}
    for q in dic_qna_p_sort:
        min_nb = 2000000000
        max_ne = 0 
        lista= []

        sorted_q_b_e_qna = OrderedDict(sorted(dic_qna_p_sort[q].items(), key=lambda x: x[1]))

        p_number = 0 
        dic_qna_cluster[q]={}
        longitud = len(dic_qna_p_sort[q].values())

        number_of_p = 0
        for p in sorted_q_b_e_qna: 

            number_of_p = number_of_p+1

                
            if sorted_q_b_e_qna[p][0]<=min_nb:
                min_nb = sorted_q_b_e_qna[p][0]
                max_ne = sorted_q_b_e_qna[p][1]
                clusterb = min_nb
                clustere = max_ne
                p_number = p_number+1
                lista.append(p)
            elif sorted_q_b_e_qna[p][0]<=max_ne:
                
                p_number = p_number+1
                lista.append(p)
                
                if sorted_q_b_e_qna[p][1]>=max_ne: 
                    clustere = sorted_q_b_e_qna[p][1]
                    max_ne = sorted_q_b_e_qna[p][1]

            
            else: 
                for element in lista: 
                    dic_qna_cluster[q][element]=str(clusterb)+"_"+str(clustere)

                lista = []
                new_cluster= "yes"
                min_nb = 2000000000
                max_ne = 0
                p_number = 0
                if sorted_q_b_e_qna[p][0]<=min_nb:
                    min_nb = sorted_q_b_e_qna[p][0]
                    max_ne = sorted_q_b_e_qna[p][1]
                    clusterb = min_nb
                    clustere = max_ne
                    p_number = p_number+1
                    lista.append(p)
                if sorted_q_b_e_qna[p][0]<=max_ne:
                    p_number = p_number+1
                    lista.append(p)

                    if sorted_q_b_e_qna[p][1]>=max_ne: 
                        clustere = sorted_q_b_e_qna[p][1]
                        max_ne = sorted_q_b_e_qna[p][1]
            if number_of_p == longitud: 
                for element in lista: 
                    if clusterb ==1:
                        clusterb="1"

                    dic_qna_cluster[q][element]=str(clusterb)+"_"+str(clustere)

    dic_cluster={}
    for q in dic_b_e:
        min_b = 2000000000
        max_e = 0 
        lista= []
        #sortedDict = sorted(dic_b_e[q].values())
        sorted_q_b_e = OrderedDict(sorted(dic_b_e[q].items(), key=lambda x: x[1]))

        p_number = 0 
        dic_cluster[q]={}
        longitud = len(dic_b_e[q].values())
        number_of_p = 0
        for p in sorted_q_b_e:

            number_of_p = number_of_p+1

                
            if sorted_q_b_e[p][0]<=min_b:
                min_b = sorted_q_b_e[p][0]
                max_e = sorted_q_b_e[p][1]
                clusterb = min_b
                clustere = max_e
                p_number = p_number+1
                lista.append(p)
            elif sorted_q_b_e[p][0]<=max_e:
                p_number = p_number+1
                lista.append(p)
                
                if sorted_q_b_e[p][1]>=max_e: 
                    clustere = sorted_q_b_e[p][1]
                    max_e = sorted_q_b_e[p][1]

                    
            
            else: 
                for element in lista: 
                    dic_cluster[q][element]=str(clusterb)+"_"+str(clustere)

                lista = []
                new_cluster= "yes"
                min_b = 2000000000
                max_e = 0
                p_number = 0
                if sorted_q_b_e[p][0]<=min_b:
                    min_b = sorted_q_b_e[p][0]
                    max_e = sorted_q_b_e[p][1]
                    clusterb = min_b
                    clustere = max_e
                    p_number = p_number+1
                    lista.append(p)
                if sorted_q_b_e[p][0]<=max_e:
                    p_number = p_number+1
                    lista.append(p)

                    if sorted_q_b_e[p][1]>=max_e: 
                        clustere = sorted_q_b_e[p][1]
                        max_e = sorted_q_b_e[p][1]
            if number_of_p == longitud: 
                for element in lista: 
                    if clusterb ==1:
                        clusterb="1"

                    dic_cluster[q][element]=str(clusterb)+"_"+str(clustere)

                
    dic_qk_freq = {}
    dic_qc_freq = {}
    cont2 = 0
    for index, row in df2.iterrows(): # A M and  L columns are added 
        d = row["d"] 
        df2.loc[cont2,"m_right"]= 0
        df2.loc[cont2,"m_left"]= 0
        try:
            if row["b"].find(";")!=-1: 
                c = dic_b_e_mn[row["q"]][row["p"]][row["pdm"]]
            else: 

                c = dic_cluster[row["q"]][row["p"]]

            df2.loc[cont2,"c"]= c.replace("ene","1")
            qc = row["q"]+":"+c 
            df2.loc[cont2,"qc"]= qc
        
            if row["n"].find(";")!=-1: 
                k = dic_qna_p_mn[row["q"]][row["p"]][row["pdm"]]
            else: 

                k = dic_qna_cluster[row["q"]][row["p"]]

            df2.loc[cont2,"k"]= k
            qk = row["q"]+":"+k
            df2.loc[cont2,"qk"]= qk

            if qk not in  dic_qk_freq.keys():
                dic_qk_freq[qk]=row["ScanFreq"]
            else:
                 dic_qk_freq[qk]=dic_qk_freq[qk]+row["ScanFreq"]

            if qc not in  dic_qc_freq.keys():
                dic_qc_freq[qc]=row["ScanFreq"]
            else:
                 dic_qc_freq[qc]=dic_qc_freq[qc]+row["ScanFreq"]

            df2.loc[cont2,"qdnaFreq"]= dic_qdna_freq[row["qdna"]]
            df2.loc[cont2,"qnaFreq"]= dic_qna_freq[row["qna"]]
            df2.loc[cont2,"qFreq"]=  dic_q_freq[row["q"]]
            df2.loc[cont2,"pFreq"]= dic_p_freq[row["p"]]
           
            if row["pdm"].find("(") != -1:
                M = dic_pd_M[row["p"]][row["d"]][1]
                L = dic_pd_M[row["p"]][row["d"]][2]
                N = dic_pd_M[row["p"]][row["d"]][3]
                qdNA = dic_pd_M[row["p"]][row["d"]][4]
                qNA = dic_pd_M[row["p"]][row["d"]][5]
                
                df2.loc[cont2,"M"] = M
                df2.loc[cont2,"A"]=row["p"][M-1]
                df2.loc[cont2,"L"]= L
                df2.loc[cont2,"N"]= N
                df2.loc[cont2,"qdNA"]= qdNA
                df2.loc[cont2,"qNA"]= qNA
            elif row["pdm"].find("[") != -1:
                M = dic_pd_M[row["p"]][row["d"]][1]
                L = dic_pd_M[row["p"]][row["d"]][2]
                N = dic_pd_M[row["p"]][row["d"]][3]
                qdNA = dic_pd_M[row["p"]][row["d"]][4]
                qNA = dic_pd_M[row["p"]][row["d"]][5]
                
                df2.loc[cont2,"M"] = M
                df2.loc[cont2,"A"]=row["p"][M-1]
                df2.loc[cont2,"L"]= L
                df2.loc[cont2,"N"]= N
                df2.loc[cont2,"qdNA"]= qdNA
                df2.loc[cont2,"qNA"]= qNA
            else:
                df2.loc[cont2,"M"] = np.nan
                df2.loc[cont2,"A"] = "U"
                df2.loc[cont2,"L"] = np.nan
                df2.loc[cont2,"N"]= np.nan
                df2.loc[cont2,"qdNA"]= row["q"]+":"+str(round(row["d"],6))+"::"+"U"
                df2.loc[cont2,"qNA"]= row["q"]+"::"+"U"
            cont2 = cont2+1
        except: 
            cont2 = cont2+1
            pass

    cont3 = 0
    for index, row in df2.iterrows():
        try:
      
            df2.loc[cont3,"qkFreq"]= dic_qk_freq[row["qk"]]
            df2.loc[cont3,"qcFreq"]= dic_qc_freq[row["qc"]]
            cont3 = cont3+1
        except: 
            cont3 = cont3+1
            pass


    logging.info("Writing output file")
    name = output[:-4]
    #add g column (keep d value if not 0, otherwise NM)
    df2["g"] = df2["d"].apply(lambda x: x if x != "0" else "NM")
    #add pgm column
    df2["pgm"] = df2.apply(build_pgm, axis=1)
    #create pgmFreq as the mean of the intensities of the samples for that precursor
    import pandas as pd #because someone thought that calling "pd" a variable was a good idea LOL
    cols_excluir = ["Proteotypic", "Precursor.Charge"]
    cols_numericas = df.select_dtypes(include="number").columns.difference(cols_excluir)
    df2["pgmFreq"] = df2["pdm"].apply(
        lambda x: (
            np.log2(
                df.loc[df[seq_column_name] == x, cols_numericas]
                  .stack()
                  .mean()
            )
            if df.loc[df[seq_column_name] == x, cols_numericas].stack().mean() > 0
            else 0
        )
    )
    df2.to_csv(name+Outfile_suffix+"_DIAversion.txt", index=False, sep='\t', encoding='utf-8')
    logging.info('end script')

if __name__ == '__main__':
    
    # multiprocessing.freeze_support()

    # parse arguments
    parser = argparse.ArgumentParser(
        description='PDMTableMaker',
        epilog='''
        Example:
            python PDMTbleMaker.py
        ''')
      
    # default DM0Solver configuration file
    defaultconfig = os.path.join(os.path.dirname(__file__), "config/Solver.ini")
    
    parser.add_argument('-i', '--infile', required=True, help='Input file with paths of the file(s) ')
    parser.add_argument('-f', '--fastafile', required=True, help='Path to input file')
    parser.add_argument('-c', '--config', default=defaultconfig, help='Path to custom config.ini file')
    parser.add_argument('-v', dest='verbose', action='store_true', help="Increase output verbosity")
    
    args = parser.parse_args()
    # logging debug level. By default, info level
    log_file = outfile = args.infile[:-4] + 'PDMTable_log.txt'
    log_file_debug = outfile = args.infile[:-4] + 'PDMTable_log_debug.txt'
    if args.verbose:
        logging.basicConfig(level=logging.DEBUG,
                            format='%(asctime)s - %(levelname)s - %(message)s',
                            datefmt='%m/%d/%Y %I:%M:%S %p',
                            handlers=[logging.FileHandler(log_file_debug),
                                      logging.StreamHandler()])
    else:
        logging.basicConfig(level=logging.INFO,
                            format='%(asctime)s - %(levelname)s - %(message)s',
                            datefmt='%m/%d/%Y %I:%M:%S %p',
                            handlers=[logging.FileHandler(log_file),
                                      logging.StreamHandler()])

    
    

    
    
    #start main function
    logging.info('start script: '+"{0}".format(" ".join([x for x in sys.argv])))
    

        
    main(args.config, args.infile, args.fastafile)