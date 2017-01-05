# -*- coding: utf-8 -*-
"""
Created on Tue Jan  3 13:49:29 2017

@author: jmorel

"""
def dilution_recc(chop_ppm,p_conc,chop_ng='',sample_id=1,target=[1,3,12,40,100],target_v=50,vol_min=5,vol_max=200,drop_null=True):
    """
    This function builds a list of reccomended dilution levels for a given sample.
    Predominatley used for the CHOP assay, it can also be used for any general
    application to take the guesswork out of dilutions.
    
    Example: You want to dilute a sample at 345ppm to a stock of 1,2,12,40, 
    100 ppm. What is the best way to do it without extensive pen/paper work?
    
    Dilution levels are incremented in whole number fractions. Ie a 1:3 has the
    next higher dilution level of 1:4. There are no fractional dilution factors.
    
    Parameters:
    -----
    chop_ppm: float
        Ppm of CHOP in starting sample. CHOP ppm is measured on a part CHOP
        to part protein of interest.
    p_conc: float
        Concentration of protein of interest. Typically a mAb product
    chop_ng: float, optional 
        Concentration of CHOP in ng/mL as measured by the standard ELISA avaiable
        from Cygnus or Perkin Elmer. CHOP in ng/mL is product of the CHOP ppm
        and the p_concentration in mg/mL
    sample_id: int, optional, default=1
        Sample ID to be recored in the resultant dataframe. Sample IDs should be
        sequential and/or unique
    target: list, optional, default=[1,3,12,40,100]
        List of target dilution factors
    target_v: float, optional, default=50
        Target volume of sample stock to be created in uL
    vol_min: float, optional, default=5
        Minimum volume in uL desired to be used for dilution series. This is
        the volume the user/machine would pipette
    vol_max: float, optional, default=200
        Maximum volume in uL desired to be used for dilution series. This is
        the volume the user/machine would pipette
    drop_null: bool, optional, default=True
        Indicator to/not drop values outside of the vol_min, vol_max ranges. True
        indicates volumes outside of range are discarded, False means the volumes
        are not dropped
    
    Output
    -----
    PD Dataframe with info
    """    
    import pandas as pd
    
    chop_ng=chop_ppm*p_conc
    dilution_table=pd.DataFrame(columns=['sample_id','chop_ppm','p_conc','chop_ng','target','dil_factor','approx_chop','sample_v','buffer_v'])

    n=0
    for i in target:
       if chop_ng%i >= 0.1*i:
           df = 1+chop_ng//(i)
       else:
           df = chop_ng//i
       if df == 0:
           break;
        
       apc = chop_ng/df
       sample_v = target_v*apc/chop_ng
       buffer_v = target_v-sample_v
       dilution_table.loc[n,['sample_id','chop_ppm','p_conc','chop_ng','target','dil_factor','approx_chop','sample_v','buffer_v']]=sample_id,chop_ppm,p_conc,chop_ng,i,df,apc,sample_v,buffer_v
       n+=1
       
    if drop_null == True:
       dilution_table=dilution_table[:][(vol_max>=dilution_table['sample_v']) & (dilution_table['sample_v']>=vol_min)]
       dilution_table=dilution_table[:][(vol_max>=dilution_table['buffer_v']) & (dilution_table['buffer_v']>=vol_min)]
    
    if len(dilution_table) == 0:
        error_warning='Intermediate dilution series may be required.'
        print('Error in sample number %s due to volume restrictions.\n' %sample_id)
        print(error_warning+'\n')
    
    return dilution_table