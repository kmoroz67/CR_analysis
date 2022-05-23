import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def tumor_clon_fraction_plot(samples,
                             time_points,
                             tumor_clone_cdr3,
                             material_type,
                             seq_type,
                             save_path,
                             max_y_lim=None,
                             cdr3xx='cdr3aa',
                             figsize=(3, 3),
                             DPI=150,
                             yticklables=True):
#     colors = ['#299BED', '#0F57C4', '#000096']
    colors = ['#D2E1F4', '#A6C4EA', '#88AADD', '#6D94CE']
    barplot_palette = dict(zip(time_points, colors))
    #     samples_dict = dict(zip(time_points, samples))
    for_table = {'Time': time_points, 'Fraction': []}
    for sample in samples:
        df = pd.read_csv(sample, sep='\t')
        df = df[df[cdr3xx] == tumor_clone_cdr3].reset_index(drop=True)
        if len(df) == 0:
            for_table['Fraction'].append(0)
        else:
            for_table['Fraction'].append(sum(df['fraction']))
    DF = pd.DataFrame(for_table)
    fig, ax = plt.subplots(1, 1, facecolor='white', figsize=figsize, dpi=DPI)
    sns.barplot(data=DF, x='Time', y='Fraction', palette=barplot_palette)
    plt.title('{} {}\n{}'.format(material_type, seq_type, tumor_clone_cdr3))
    if max_y_lim == None:
        y_ax_list = ax.get_yticks()
        ax.set_ylim([0, max(ax.get_yticks())])
    else:
        ax.set_ylim([0, max_y_lim])
    ax.set_xlabel('')
    ax.grid(False)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=60)
    if not yticklables:
        ax.set_ylabel('')
        ax.set_yticklabels([])
        ax.get_yaxis().set_ticks([])
        plt.tick_params(left=False)
    plt.savefig('{}tumor_clone_{}_{}_{}.svg'.format(save_path, material_type, seq_type, tumor_clone_cdr3), bbox_inches='tight')


def tumor_clon_fraction_multyplot(samples_lists,
                                  time_points,
                                  tumor_clone_cdr3,
                                  material_types,
                                  seq_types,
                                  save_path,
                                  cdr3xx='cdr3aa'):
    yticklables = [False] * len(samples_lists)
    yticklables[0] = True
    all_samples = sum(samples_lists, [])
    tumor_fraction = []
    for sample in all_samples:
        df = pd.read_csv(sample, sep='\t')
        tumor_fraction.append(
            sum(df[df[cdr3xx] == tumor_clone_cdr3]['fraction']))
    max_y_lim_df = pd.DataFrame({'x': ['x'], 'y': [max(tumor_fraction)]})
    fig, ax = plt.subplots(1, 1)
    sns.barplot(data=max_y_lim_df, x='x', y='y', ax=ax)
    max_y_lim = max(ax.get_yticks())
    del ax
    del fig
    for i in range(len(samples_lists)):
        tumor_clon_fraction_plot(samples_lists[i],
                                 time_points,
                                 tumor_clone_cdr3,
                                 material_types[i],
                                 seq_types[i],
                                 save_path,
                                 max_y_lim,
                                 cdr3xx=cdr3xx,
                                 figsize=(3, 3),
                                 DPI=150,
                                 yticklables=yticklables[i])


# samples_lists = [
#     [
#         '/uftp/Blood/db_calc_pipeline/BG001082_december_sum/PBMC_BULK_April_BCR_qc_pass.txt',
#         '/uftp/Blood/db_calc_pipeline/BG001082_december_sum/PBMC_BULK_August_BCR_qc_pass.txt',
#         '/uftp/Blood/db_calc_pipeline/BG001082_december_sum/PBMC_BULK_BCR_qc_pass.txt'
#     ],
#     [
#         '/uftp/Blood/db_calc_pipeline/BG001082_december_sum/B_Cells_BULK_April_BCR_qc_pass.txt',
#         '/uftp/Blood/db_calc_pipeline/BG001082_december_sum/B_Cells_BULK_August_BCR_qc_pass.txt',
#         '/uftp/Blood/db_calc_pipeline/BG001082_december_sum/Total_B_Cells_BULK_BCR_qc_pass.txt'
#     ]
# ]
# seq_types = ['BULK', 'BULK']
# material_types = ['PBMC', 'B Cells']
# time_points = ['April', 'August', 'December']
# tumor_clone_cdr3 = 'CQQHHSSPYTF'
# save_path = '/uftp/Blood/db_calc_pipeline/BG001082_december_sum/'

# tumor_clon_fraction_multyplot(samples_lists,
#                               time_points,
#                               tumor_clone_cdr3,
#                               material_types,
#                               seq_types,
#                               save_path,
#                               cdr3xx='cdr3aa')