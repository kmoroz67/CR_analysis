import pandas as pd
import plotly.express as px

def top_clonotypes_plot(save_path,
                        samples,
                        time_points,
                        chain,
                        seq_type,
                        material_type,
                        cdr3xx,
                        top,
                        column='all',
                        legend=True):
    DF = pd.DataFrame()
    cdr3 = []
    if 'Ig' in chain:
        ch = 'isotype'
        title = '{} {}<br>Isotype: {}'.format(seq_type, material_type, chain)
    else:
        ch = 'chain'
        title = '{} {}<br>Chain: {}'.format(seq_type, material_type,
                                            chain.capitalize())
    samples_dict = dict(zip(time_points, samples))
    for time_point in time_points:
        df = pd.read_csv(samples_dict[time_point],
                         sep='\t',
                         usecols=[cdr3xx, ch, 'fraction'])
        if chain == 'light':
            df = df.query("chain=='kappa' | chain=='lambda'").reset_index(
                drop=True)
        else:
            df = df[df[ch] == chain].reset_index(drop=True)
        df = df.groupby(cdr3xx)['fraction'].apply(sum).reset_index(
            name='Cumulative abundance')
        df.sort_values(by=['Cumulative abundance'],
                       ascending=False,
                       inplace=True)
        df['Timepoint'] = [time_point] * df.shape[0]
        DF = pd.concat([DF, df], ignore_index=True)
        if column == 'all':
            cdr3 += list(df[cdr3xx][:top])
            cdr3 = list(set(cdr3))
        else:
            if time_point == column:
                cdr3 = list(df[cdr3xx][:top])
    DF.rename(columns={cdr3xx: cdr3xx.upper()[:-2]}, inplace=True)
    cdr3xx = cdr3xx.upper()[:-2]
    DF = DF[DF[cdr3xx].isin(cdr3)]
    for time_point in time_points:
        DF_time_point = DF[DF['Timepoint'] == time_point]
        if DF_time_point.shape[0] < len(cdr3):
            cdr3_zero = [
                x for x in cdr3 if x not in list(DF_time_point[cdr3xx])
            ]
            df_zero = pd.DataFrame({
                cdr3xx: cdr3_zero,
                'Cumulative abundance': [0] * len(cdr3_zero),
                'Timepoint': [time_point] * len(cdr3_zero)
            })
            DF = pd.concat([DF, df_zero], ignore_index=True)
    if legend:
        width = 920
        height = 800
    else:
        width = 800
        height = 800
    fig = px.area(DF,
                  x="Timepoint",
                  y="Cumulative abundance",
                  color=cdr3xx,
                  title=title,
                  width=width,
                  height=height)
    fig.update_layout(showlegend=legend,
                      font=dict(family="Source Sans Pro", size=16))
    fig.write_image('{}{}_{}_top_{}_timepoints_{}_{}.svg'.format(
        save_path, material_type, seq_type, top, column, chain))
    

# data example
save_path = '/home/akosenkov/Metrics/'
legend = True
top = 5
cdr3xx = 'cdr3aa'
chain = 'light'
samples = [
    '/uftp/Blood/db_calc_pipeline/BG001082_december_sum/B_Cells_BULK_April_BCR_qc_pass.txt',
    '/uftp/Blood/db_calc_pipeline/BG001082_december_sum/B_Cells_BULK_August_BCR_qc_pass.txt',
    '/uftp/Blood/db_calc_pipeline/BG001082_december_sum/B_Cells_BULK_BCR_qc_pass.txt'
]
time_points = ['April', 'August', 'December']
seq_type = 'BULK'
material_type = 'PBMC'
# column = 'April'
column = 'all'

# launch example
fig = top_clonotypes_plot(save_path,
                          samples,
                          time_points,
                          chain,
                          seq_type,
                          material_type,
                          cdr3xx,
                          top,
                          column,
                          legend=True)