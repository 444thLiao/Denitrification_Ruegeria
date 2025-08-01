import pandas as pd
S2S = {'YTW':'NM1',
       'PC':'SM10',
       #'BI': 'PM11',
       'PI': 'MM17',
       'NP': 'MM19',
       'LC': 'SM1',
       'SW': 'SM4',
       'SI': 'PM3'
       }
    
# merged_f = []
# for f in glob('/mnt/ivy/thliao/project/coral_ruegeria/EPD_data/raw/*.csv'):
#     df = pd.read_csv('../EPDdata/env_data',sep=',')
#     merged_f.append(df)
water_df = pd.read_csv("D:\\Desktop\\nan papers\\EPDdata\\All_stations.tsv",sep='\t')
water_df.loc[:,'Dates'] = pd.to_datetime(water_df['Dates'])

def r_symbol(i):
    if '<' in str(i):
        i = i.replace('<','').strip()
    elif '>' in str(i):
        i = i.replace('>','').strip()
    return i
water_df = water_df.applymap(r_symbol)
#water_df.to_csv('/mnt/ivy/thliao/project/coral_ruegeria/EPD_data/All_stations.tsv',sep='\t')

subwater_df = water_df.loc[water_df['Station'].isin(S2S.values()),:]

from datetime import datetime
subwater_df = subwater_df.loc[(subwater_df.Dates >= datetime(2010, 5, 1)),:]
subwater_df = subwater_df.loc[~subwater_df.isna().any(axis=1),:]
subwater_df = subwater_df.reindex(columns=[ 
'Station', 'Total Nitrogen (mg/L)','Total Inorganic Nitrogen (mg/L)','Nitrite Nitrogen (mg/L)','Nitrate Nitrogen (mg/L)','Ammonia Nitrogen (mg/L)', 'Total Phosphorus (mg/L)','Turbidity (NTU)',  'Salinity (psu)', 'pH',])

# for i in subwater_df.columns[5:]:
#     subwater_df.loc[:,i] = subwater_df[i].astype(float)
mean_df = subwater_df.groupby('Station').mean()
mean_df.insert(0, 'EPD Station',list(mean_df.index))
rev_S2S = {v:k for k,v in S2S.items()}
mean_df.index = [rev_S2S[_] for _ in mean_df.index]
mean_df = mean_df.reindex(index=['YTW',
'PC',
'SW',
'LC',
'NP',
'SI',
'PI',])


mean_df.loc[:,'Total Nitrogen (mM/L)'] = mean_df['Total Nitrogen (mg/L)']/14
mean_df.loc[:,'Total Inorganic Nitrogen (mM/L)'] = mean_df['Total Inorganic Nitrogen (mg/L)']/14
mean_df.loc[:,'Nitrite Nitrogen (mM/L)'] = mean_df['Nitrite Nitrogen (mg/L)']/46
mean_df.loc[:,'Nitrate Nitrogen (mM/L)'] = mean_df['Nitrate Nitrogen (mg/L)']/62
mean_df.loc[:,'Ammonia Nitrogen (mM/L)'] = mean_df['Ammonia Nitrogen (mg/L)']/17
mean_df.loc[:,'Total Phosphorus (mM/L)'] = mean_df['Total Phosphorus (mg/L)']/30
mean_df.loc[:,'N:P ratio'] = mean_df['Total Nitrogen (mM/L)']/mean_df['Total Phosphorus (mM/L)']
mean_df.to_excel('../figuresandtables/Table S1(mMversion).xlsx')


