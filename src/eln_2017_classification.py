def eln_2017_classification(df_merge):
    """
    This function builds the ELN 2017 classification of AML for the dataframe given in input. The classification is a new column labeled 'eln_2017'

    It supposes that the following are columns of the input:
        - t_8_21, inv_16, t_16_16, CEBPA_bi, NPM1, FLT3_ITD t_9_11, TP53, t_6_9, MLL, t_9_22, inv_3, t_3_3,
        minus5, del_5q, minus7, minus17, RUNX1 (or hotspots starting with RUNX1), ASXL1 (ibid), complex


    :param df_final:
    :return: pandas.DataFrame, copy of df_final with 'eln_2017' column
    """
    import pandas as pd
    import numpy as np
    df = df_merge.copy(deep=True)
    df.loc[:, 'eln_2017'] = np.nan
    df.loc[(df.eln_2017.isnull()) & (df.inv_16 == 1), 'eln_2017'] = 'favorable'
    df.loc[(df.eln_2017.isnull()) & (df.t_15_17 == 1), 'eln_2017'] = 'favorable'
    df.loc[(df.eln_2017.isnull()) & (df.t_9_11 == 1), 'eln_2017'] = 'intermediate'
    df.loc[(df.eln_2017.isnull()) & (df.complex == 1), 'eln_2017'] = 'adverse'
    df.loc[(df.eln_2017.isnull()) & (df.TP53 == 1), 'eln_2017'] = 'adverse'
    df.loc[(df.eln_2017.isnull()) & (df.t_6_9 == 1), 'eln_2017'] = 'adverse'
    if 'MLL_PTD' in df.columns.values:
        df.loc[(df.eln_2017.isnull()) & (df.MLL_PTD == 1), 'eln_2017'] = 'adverse'
    else:
        df.loc[(df.eln_2017.isnull()) & (df.MLL == 1), 'eln_2017'] = 'adverse'
    df.loc[(df.eln_2017.isnull()) & (df.t_9_22 == 1), 'eln_2017'] = 'adverse'
    df.loc[(df.eln_2017.isnull()) & (df.inv_3 == 1), 'eln_2017'] = 'adverse'
    df.loc[(df.eln_2017.isnull()) & (
                (df.minus5 == 1) | (df.del_5q == 1) | (df.minus7 == 1) | (df.minus17 == 1)), 'eln_2017'] = 'adverse'
    df.loc[(df.eln_2017.isnull()) & (df.NPM1 == 0) & (df.ITD == 1), 'eln_2017'] = 'adverse'  #yanis change for my dataset FLT3_ITD to new ITD

    df.loc[(df.eln_2017.isnull()) & (df.t_8_21 == 1), 'eln_2017'] = 'favorable'
    #df.loc[(df.eln_2017.isnull()) & (df.inv_16 == 1), 'eln_2017'] = 'favorable'
    df.loc[(df.eln_2017.isnull()) & (df.CEBPA_bi == 1), 'eln_2017'] = 'favorable'
    df.loc[(df.eln_2017.isnull()) & (df.NPM1 == 1) & (df.ITD == 0), 'eln_2017'] = 'favorable' #yanis change for my dataset FLT3_ITD to new ITD

    # As we distinguish between hotspot mutations,  several columns may encode
    # for a gene mutation

    for col in [c for c in df.columns if c.startswith('RUNX1')]:
        df.loc[(df.eln_2017.isnull()) & (df[col] == 1), 'eln_2017'] = 'adverse'
    for col in [c for c in df.columns if c.startswith('ASXL1')]:
        df.loc[(df.eln_2017.isnull()) & (df[col] == 1), 'eln_2017'] = 'adverse'

    for col in [x for x in df.columns if (x.startswith('minus') or x.startswith('del'))]:
        df.loc[(df.eln_2017.isnull()) & (df[col] == 1), 'eln_2017'] = 'adverse'

    df.loc[(df.eln_2017.isnull()) & (df.NPM1 == 0) & (df.ITD == 0), 'eln_2017'] = 'intermediate'   #yanis change for my dataset FLT3_ITD to new ITD
    df.loc[(df.eln_2017.isnull()) & (df.NPM1 == 1) & (df.ITD == 1), 'eln_2017'] = 'intermediate'
    
    df.loc[(df.eln_2017.isnull()) & (df.t_v_11 == 1), 'eln_2017'] = 'adverse'  # latest change
    
    df.loc[df.eln_2017.isnull(), 'eln_2017'] = 'intermediate'

    return df.eln_2017
