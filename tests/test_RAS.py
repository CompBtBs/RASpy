# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import pandas as pd
import cobra as cb
from cobra import Model,Reaction,Metabolite
import pytest
import numpy as np
import sys
from scanpy import AnnData
from ras import RAS_computation as rc
sys.path.insert(0, '../raspy')

def test_gpr1():
    
    ### check the computation of gpr


    model = Model('example_model')
    
    #single gene gpr
    reaction = Reaction('A')
    reaction.gene_reaction_rule = 'GeneA'
    model.add_reactions([reaction])
    
    #or gpr
    reaction = Reaction('B')
    reaction.gene_reaction_rule = 'GeneA or GeneB'
    model.add_reactions([reaction])
    
    #and gpr
    reaction = Reaction('C')
    reaction.gene_reaction_rule = 'GeneB and GeneC'
    model.add_reactions([reaction])
    
    #complex gpr 1
    reaction = Reaction('D')
    reaction.gene_reaction_rule = '(GeneA or GeneB) and GeneC'
    model.add_reactions([reaction])
 
    #complex gpr 2
    reaction = Reaction('E')
    reaction.gene_reaction_rule = '( GeneA or GeneB ) and (GeneC and GeneD)'
    model.add_reactions([reaction])    
 
    ### create count matrix
    num_cells=5
    data=pd.DataFrame(index=["GeneA","GeneB","GeneC","GeneD"],columns=["Cell"+str(i) for i in range(num_cells)])    
    data.iloc[:,:]=np.random.randint(low=0,high=10,size=(4,num_cells))
    adata=AnnData(data.T)
         
    ras_object=rc(adata,model)
    ras_adata=ras_object.compute()
    
    df=ras_adata.to_df().T
    
    
    assert df.shape==(5,num_cells)
    
    assert all(df.loc["A"]==data.loc["GeneA"])
    assert all(df.loc["B"]==data.loc["GeneA"]+data.loc["GeneB"])
    
    assert all(df.loc["C"]==np.min([data.loc["GeneB"],data.loc["GeneC"]],axis=0))
    
    val=data.loc["GeneA"]+data.loc["GeneB"]

    assert all(df.loc["D"]==np.min([val,data.loc["GeneC"]],axis=0))


    val1=data.loc["GeneA"]+data.loc["GeneB"]
    val2=np.min([data.loc["GeneC"],data.loc["GeneD"]],axis=0)

    assert all(df.loc["E"]==np.min([val1,val2],axis=0))
    
def test_gpr2():
    
    ### Check nan in the GPR


    model = Model('example_model')
    
    #single gene gpr
    reaction = Reaction('A')
    reaction.gene_reaction_rule = 'GeneA'
    model.add_reactions([reaction])
    
    #or gpr
    reaction = Reaction('B')
    reaction.gene_reaction_rule = 'GeneA or GeneE'
    model.add_reactions([reaction])
    
    #and gpr
    reaction = Reaction('C')
    reaction.gene_reaction_rule = 'GeneB and GeneE'
    model.add_reactions([reaction])
    

    ### create count matrix
    num_cells=5
    data=pd.DataFrame(index=["GeneA","GeneB","GeneC","GeneD"],columns=["Cell"+str(i) for i in range(num_cells)])    
    data.iloc[:,:]=np.random.randint(low=0,high=10,size=(4,num_cells))
    adata=AnnData(data.T)
         
    ras_object=rc(adata,model)
    ras_adata=ras_object.compute()
    
    df=ras_adata.to_df().T
    

    assert df.shape==(3,num_cells)
    
    assert all(df.loc["A"]==data.loc["GeneA"])
    assert all(df.loc["B"]==data.loc["GeneA"])
    assert all(df.loc["C"]==data.loc["GeneB"])
    
    
    ras_adata=ras_object.compute(and_expression=np.min,
                                 or_expression=np.sum,
                                 drop_na_rows=False)
    
    df=ras_adata.to_df().T

    assert df.shape==(3,num_cells)
    assert all(df.loc["B"].isna()) 
    assert all(df.loc["C"].isna()) 


def test_gpr3():
    
    ### check the different operation to assemble the gpr


    model = Model('example_model')
    

    #or gpr
    reaction = Reaction('B')
    reaction.gene_reaction_rule = 'GeneA or GeneB'
    model.add_reactions([reaction])
    
    #and gpr
    reaction = Reaction('C')
    reaction.gene_reaction_rule = 'GeneB and GeneC'
    model.add_reactions([reaction])
    
 
    ### create count matrix
    num_cells=5
    data=pd.DataFrame(index=["GeneA","GeneB","GeneC","GeneD"],columns=["Cell"+str(i) for i in range(num_cells)])    
    data.iloc[:,:]=np.random.randint(low=0,high=10,size=(4,num_cells))
    adata=AnnData(data.T)
         
    ras_object=rc(adata,model)
    ras_adata=ras_object.compute(and_expression=np.nanmean,
                                 or_expression=np.nanmean)
    
    df=ras_adata.to_df().T

    assert df.shape==(2,num_cells)    
    assert all(df.loc["B"]==(data.loc["GeneA"]+data.loc["GeneB"])/2)    
    assert all(df.loc["C"]==(data.loc["GeneB"]+data.loc["GeneC"])/2)   
    