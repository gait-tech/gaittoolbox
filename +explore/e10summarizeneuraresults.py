import pandas as pd
import numpy as np

import plotly.plotly as py
import plotly.graph_objs as go
import plotly

result_mat_path = '../neura/explore/results2.csv'
summary_output_path = '../neura/explore/summary.csv'

def get_action_type(s):
    return s.split('-')[-2]
df = pd.read_csv(result_mat_path)
df['type'] = df['name'].apply(lambda s: s.split('-')[-2])

df['anklePosMean'] = np.mean(df[['LTIO_1', 'LTIO_2', 'LTIO_3', 'RTIO_1', 'RTIO_2', 'RTIO_3']], axis=1)
df['anklePosStd'] = np.mean(df[['LTIOStd_1', 'LTIOStd_2', 'LTIOStd_3', 'RTIOStd_1', 'RTIOStd_2', 'RTIOStd_3']], axis=1)
df['ankleLPosMean'] = np.mean(df[['LTIO_1', 'LTIO_2', 'LTIO_3']], axis=1)
df['ankleLPosStd'] = np.mean(df[['LTIOStd_1', 'LTIOStd_2', 'LTIOStd_3']], axis=1)
df['ankleRPosMean'] = np.mean(df[['RTIO_1', 'RTIO_2', 'RTIO_3']], axis=1)
df['ankleRPosStd'] = np.mean(df[['RTIOStd_1', 'RTIOStd_2', 'RTIOStd_3']], axis=1)
df['kneePosMean'] = np.mean(df[['LFEO_1', 'LFEO_2', 'LFEO_3', 'RFEO_1', 'RFEO_2', 'RFEO_3']], axis=1)
df['kneePosStd'] = np.mean(df[['LFEOStd_1', 'LFEOStd_2', 'LFEOStd_3', 'RFEOStd_1', 'RFEOStd_2', 'RFEOStd_3']], axis=1)
df['kneeOriYMean'] = np.mean(df[['qLKNE_2', 'qRKNE_2']], axis=1)
df['kneeOriYStd'] = np.mean(df[['qLKNEStd_2', 'qRKNEStd_2']], axis=1)
df['kneeLOriYMean'] = np.mean(df[['qLKNE_2']], axis=1)
df['kneeLOriYStd'] = np.mean(df[['qLKNEStd_2']], axis=1)
df['kneeROriYMean'] = np.mean(df[['qRKNE_2']], axis=1)
df['kneeROriYStd'] = np.mean(df[['qRKNEStd_2']], axis=1)

df_result = df.groupby(['type', 'label']).mean()
# df_result.to_csv(summary_output_path)
print(df_result.head())

traces_name = {'nonlin-fmincon': 'Nssv+M02+C135', 'nonlin-SCKF': 'Nssv+M02+C152', 
               'lin-world-frame': 'Nssv+M02+C202', 'lin-pelv-frame': 'Nssv+M02+C203'}

target_params = ['anklePos', 'ankleLPos', 'ankleRPos', 'kneeOriY', 'kneeLOriY', 'kneeROriY']

for target in target_params:
    traces = []
    for k in traces_name:
        v = traces_name[k]
        idx = (df_result.index.get_level_values('label') == v)
        df_buf = df_result.ix[idx]
   
        trace = go.Bar(
            x=df_buf['{}Mean'.format(target)].index.get_level_values('type').values,
            y=df_buf['{}Mean'.format(target)].values,
            name=k,
            error_y=dict(
                type='data',
                array=df_buf['{}Std'.format(target)].values,
                visible=True
            )
        )
        traces.append(trace)

    layout = go.Layout(
        title='{}'.format(target),
        xaxis=dict(
            title='Types of motion',
        ),
        yaxis=dict(
            title='rmse (m)',
        ),
        barmode='group'
    )
    fig = go.Figure(data=traces, layout=layout)
    plotly.offline.plot(fig, filename='{}-bar.html'.format(target))
