import numpy as np
import pandas as pd
import os


def run(
        probes_to_fragments_path: str,
        samples_names: np.ndarray | list,
        wt_capture_eff_dir: str,
        sparse_mat_path: str,
        binned_contacts_path: str,
        statistics_path: str,
        output_dir: str):

    df_probes = pd.read_csv(probes_to_fragments_path, sep='\t', index_col=0).T
    df_sparse_mat = pd.read_csv(sparse_mat_path, header=0, sep="\t", names=['frag_a', 'frag_b', 'contacts'])
    df_binned_contacts = pd.read_csv(binned_contacts_path, header=0, sep="\t", index_col=0)
    df_stats = pd.read_csv(statistics_path, header=0, sep="\t", index_col=0)

    ref_wt: dict = \
        {k.lower(): pd.read_csv(wt_capture_eff_dir+k, sep='\t') for k in os.listdir(wt_capture_eff_dir)}

    probes_dsdna_frag = df_probes.loc[df_probes['type'] == 'ds']['frag_id'].values
    probes_dsdna_name = df_probes.loc[df_probes['type'] == 'ds'].columns.values

    probes_dsdnaneg_frag = df_probes.loc[df_probes['type'] == 'ds']['frag_id'].values
    probes_dsdnaneg_name = df_probes.loc[df_probes['type'] == 'ds_neg'].columns.values

    probes_ssdna_frag = df_probes.loc[df_probes['type'] == 'ds']['frag_id'].values
    probes_ssdna_name = df_probes.loc[df_probes['type'] == 'ss'].columns.values

    probes_name = np.concatenate((probes_dsdna_name, probes_dsdnaneg_name, probes_ssdna_name))
    probes_frag = np.concatenate((probes_dsdna_frag, probes_dsdnaneg_frag, probes_ssdna_frag))

    total_contacts = sum(df_sparse_mat["contacts"])
    pondered_dsdna_mut = np.empty(shape=np.shape(probes_dsdna_frag))
    pondered_dsdna_wt = np.empty(shape=np.shape(probes_dsdna_frag))
    df_stats = df_stats.assign(contacts_divided_by_total=(df_stats["total"] / total_contacts))

    df_wt2h = ref_wt['wt2h.tsv']
    df_wt2h = ref_wt['wt4h.tsv']

    # loop calculate average dsdna capture efficiency in mutant and WT, and the correction factor mut/WT
    counter = 0
    for frag in probes_dsdna_frag:
        pondered_dsdna_mut[counter] = df_stats.loc[df_stats['fragments'] == int(frag), 'total'].tolist()[0]
        pondered_dsdna_wt[counter] = df_wt2h.loc[df_wt2h['fragments'] == int(frag)]["contacts_over_total"].tolist()[0]
        counter += 1
    dsdna_correct_fact = np.mean(pondered_dsdna_mut) / np.mean(pondered_dsdna_wt)

    # Apply the dsdna capture efficiency to the fraction of total contacts and add it to the stats
    mut_efficiency = df_stats["total"] / np.mean(pondered_dsdna_mut)
    df_stats = df_stats.assign(dsdna_norm_capture_efficiency=mut_efficiency)
    df_stats = df_stats.assign(Capture_efficiency_norm_WT=mut_efficiency)

    # divide the dsdna-normalized contacts to the WT_capture efficiency for each probe
    # to get the correction factor for each probe
    for frag in df_stats["fragments"]:
        mut_value = df_stats.loc[df_stats['fragments'] == int(frag), "dsdna_norm_capture_efficiency"].tolist()[0]
        wt_value = df_wt2h.loc[df_wt2h["fragments"] == int(frag), "Capture_efficiency_WT"].tolist()[0]
        df_stats.loc[df_stats['fragments'] == int(frag), "Capture_efficiency_norm_WT"] = mut_value / wt_value

    df_binned_contacts_pondered = df_binned_contacts.copy(deep=True)

    for probe in df_stats['probes']:
        frag_id = df_probes.loc[probe, 'frag_id']
        if frag_id not in df_binned_contacts.columns.values:
            continue
        x = df_stats.loc[df_stats['probes'] == probe, 'Capture_efficiency_norm_WT'].tolist()[0]
        df_binned_contacts_pondered.loc[:, frag_id] = df_binned_contacts.loc[:, frag_id].multiply(x)

    pass

#     for binn in binning:
#         contacts = pd.read_csv(str(contacts_dir + binn + "/" + sample + "_" + binn + "_frequencies.tsv"), header=0,
#                                sep="\t")
#         header = contacts.loc[0:4, ]
#         contacts = contacts.drop(contacts.index[0:5], axis=0)
#         contacts = contacts.astype(listtype)
#         contacts_pondered = contacts
#         for probe in stats["fragments"]:
#             # print(probe)
#             # print(float(stats.loc[(stats['fragments'] == probe), "Capture_efficiency_norm_WT"]))
#             # contacts_pondered[str(probe)] = contacts[str(probe)]/float(stats.loc[(stats['fragments'] == probe), "Capture_efficiency_norm_WT"])
#             a = float(stats.loc[(stats['fragments'] == probe), "Capture_efficiency_norm_WT"])
#             # print(a)
#             # print((contacts_pondered.loc[:, str(probe)])*float(stats.loc[(stats['fragments'] == probe), "Capture_efficiency_norm_WT"]))
#             contacts_pondered.loc[:, str(probe)] = contacts.loc[:, str(probe)].multiply(a)
#             # print(float(stats.loc[(stats['fragments'] == probe), "Capture_efficiency_norm_WT"]))
#         contacts_pondered = pd.concat((header, contacts_pondered))
#         contacts_pondered.to_csv(
#             str(contacts_dir + binn + '/' + sample + '_' + binn + '_frequencies_pondered_over_WT.tsv'), sep="\t",
#             header=True, index=False)
#
# binning = ["0kb"]
#
# for sample in samples:
#     sparse_mat = pd.read_csv(str(sparse_mat_dir + sample + "_S288c_DSB_LY_Capture_artificial_cutsite_q30.txt"),
#                              header=0, sep="\t", names=['Frag_A', 'Frag_B', 'Contacts'])
#     stats = pd.read_csv(str(stats_dir + sample + "_global_statistics.tsv"), header=0, sep="\t")
#     stats_pondered = stats
#     total_contacts = sum(
#         sparse_mat["Contacts"])  # from sparse_mat: get total contacts from which probes enrichment is calculated
#     pondered_dsdna_mut = np.empty(shape=np.shape(probes_dsdna_pos))
#     pondered_dsdna_WT = np.empty(shape=np.shape(probes_dsdna_pos))
#     stats = stats.assign(contacts_divided_by_total=(stats["total"] / total_contacts))
#
#     # loop calculate average dsdna capture efficiency in mutant and WT, and the correction factor mut/WT
#     a = 0
#     for pos in probes_dsdna_pos:
#         pondered_dsdna_mut[a] = float(stats[(stats['fragments'] == pos)]["total"])
#         pondered_dsdna_WT[a] = float(ref_WT[(ref_WT['fragments'] == pos)]["contacts_over_total"])
#         a = a + 1
#     dsdna_correct_fact = np.mean(pondered_dsdna_mut) / np.mean(pondered_dsdna_WT)
#
#     # Apply the dsdna capture efficiency to the fraction of total contacts and add it to the stats
#     mut_efficiency = stats["total"] / np.mean(pondered_dsdna_mut)
#     stats = stats.assign(dsdna_norm_capture_efficiency=mut_efficiency)
#     stats = stats.assign(Capture_efficiency_norm_WT=mut_efficiency)
#
#     # divide the dsdna-normalized contacts to the WT_capture efficiency for each probe to get the correction factor for each probe
#     for probe in stats["fragments"]:
#         mut_value = float(stats.loc[(stats['fragments'] == probe), "dsdna_norm_capture_efficiency"])
#         WT_value = float(ref_WT.loc[(ref_WT["fragments"] == probe), "Capture_efficiency_WT"])
#         stats.loc[(stats['fragments'] == probe), "Capture_efficiency_norm_WT"] = float(mut_value / WT_value)
#     stats.to_csv(str(stats_dir + sample + "_global_statistics_dsdna_and_WT_pondered.tsv"), sep="\t", header=True,
#                  index=False)
#
#     # create a list that specifies the data type "float" for each probe loaded in the contacts dataframe
#     listtype = {'chr': 'str', 'positions': 'int64'}
#     for probe in stats['fragments']:
#         listtype[str(probe)] = 'float64'
#
#     for binn in binning:
#         contacts = pd.read_csv(str(contacts_dir + binn + "/" + sample + "_" + binn + "_frequencies.tsv"), header=0,
#                                sep="\t")
#         header = contacts.loc[0:4, ]
#         contacts = contacts.drop(contacts.index[0:5], axis=0)
#         contacts = contacts.astype(listtype)
#         contacts_pondered = contacts
#         for probe in stats["fragments"]:
#             # print(probe)
#             # print(float(stats.loc[(stats['fragments'] == probe), "Capture_efficiency_norm_WT"]))
#             # contacts_pondered[str(probe)] = contacts[str(probe)]/float(stats.loc[(stats['fragments'] == probe), "Capture_efficiency_norm_WT"])
#             a = float(stats.loc[(stats['fragments'] == probe), "Capture_efficiency_norm_WT"])
#             # print(a)
#             # print((contacts_pondered.loc[:, str(probe)])*float(stats.loc[(stats['fragments'] == probe), "Capture_efficiency_norm_WT"]))
#             contacts_pondered.loc[:, str(probe)] = contacts.loc[:, str(probe)].multiply(a)
#             # print(float(stats.loc[(stats['fragments'] == probe), "Capture_efficiency_norm_WT"]))
#         contacts_pondered = pd.concat((header, contacts_pondered))
#         contacts_pondered.to_csv(
#             str(contacts_dir + binn + '/' + sample + '_' + binn + '_frequencies_pondered_over_WT.tsv'), sep="\t",
#             header=True, index=False)
