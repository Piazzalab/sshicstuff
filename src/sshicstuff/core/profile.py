"""
This module contains functions to plot 4C-like profiles and to organize the contacts made by each probe with the genome.
"""

import os
import re
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.io as pio
from plotly.subplots import make_subplots

import sshicstuff.core.graph as graph
import sshicstuff.core.methods as methods
import sshicstuff.log as log

logger = log.logger


def plot_profiles(
    profile_contacts_path: str,
    oligo_capture_path: str,
    chr_coord_path: str,
    output_dir: str = None,
    extension: str = "pdf",
    rolling_window: int = 1,
    region: str = None,
    log_scale: bool = False,
    user_y_min: float = None,
    user_y_max: float = None,
    width: int = 1200,
    height: int = 600,
):
    """
    Plot 4C-like profiles.
    """

    profile_type = "contacts"
    if "frequencies" in profile_contacts_path:
        profile_type = "frequencies"

    df: pd.DataFrame = pd.read_csv(profile_contacts_path, sep="\t")
    frags_col = df.filter(regex=r"^\d+$|^\$").columns.to_list()
    df_oligo: pd.DataFrame = pd.read_csv(oligo_capture_path, sep=",")
    probes_to_frag = dict(
        zip(df_oligo["fragment"].astype(str), df_oligo["name"].astype(str))
    )
    df_coords = pd.read_csv(chr_coord_path, sep="\t")

    full_genome_size = df_coords.loc[:, "length"].sum()
    x_min = 0
    x_max = full_genome_size

    b = re.search(r"_(\d+)kb_profile_", profile_contacts_path).group(1)
    if b == 0:
        binsize = 0
        bin_suffix = "0kb"
    else:
        binsize = int(b) * 1000  # kb to bp
        bin_suffix = f"{b}kb"

    sample_name = os.path.basename(profile_contacts_path).split(f"_{bin_suffix}_")[0]

    # Create the output directory
    if not output_dir:
        output_dir = os.path.dirname(profile_contacts_path)

    output_dir = os.path.join(output_dir, "plots")
    output_dir = os.path.join(output_dir, bin_suffix)
    if log_scale:
        output_dir = os.path.join(output_dir, "log")
    else:
        output_dir = os.path.join(output_dir, "raw")
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    chr_list_unique = pd.unique(df.chr)

    # Set coordinates preference
    coord_mode = ["genomic", "unbinned"]
    x_col = "genome_start" if binsize == 0 else "genome_bins"
    if region:
        coord_mode[0] = "chromosomal"
        region_list = region.split("-")
        if len(region_list) == 1:
            chr_, start_, end_ = region_list[0], "", ""
        else:
            chr_, start_, end_ = region_list
        max_chr_region_len = df_coords.loc[df_coords.chr == chr_, "length"].values[0]
        x_col = "start" if binsize == 0 else "chr_bins"

        if not start_ and not end_:
            df = df[df["chr"] == chr_]
            df_coords = df_coords[df_coords["chr"] == chr_]
            x_min = 0
            x_max = max_chr_region_len
        else:
            start_, end_ = int(start_), int(end_)
            x_min = start_
            x_max = end_

        df = df[(df["chr"] == chr_) & (df[x_col] >= x_min) & (df[x_col] <= x_max)]
        x_col = "start" if binsize == 0 else "chr_bins"

    if binsize > 0:
        df_bins = graph.build_bins_template(df_coords, binsize)
        if region:
            df_bins = df_bins[df_bins["chr"] == chr_]

        coord_mode[1] = "binned"
        x_min = x_min // binsize * binsize
        x_max = x_max // binsize * binsize + binsize
        df_bins = df_bins[
            (df_bins["chr_bins"] >= x_min) & (df_bins["chr_bins"] <= x_max)
        ]

        if rolling_window > 1:
            for chr_ in chr_list_unique:
                df.loc[df["chr"] == chr_, frags_col] = (
                    df.loc[df["chr"] == chr_, frags_col]
                    .rolling(window=rolling_window, min_periods=1)
                    .mean()
                )

    y_min = float(user_y_min) if user_y_min else 0.0
    y_max = float(user_y_max) if user_y_max else df[frags_col].max().max()

    log_status = ""
    if log_scale:
        log_status = "log_"
        data = df[frags_col].values
        data[data == 0] = np.nan
        new_data = np.log10(data)
        y_min = np.nanmin(new_data) if not user_y_max else float(user_y_max)
        y_max = np.nanmax(new_data) if not user_y_min else float(user_y_min)
        df[frags_col] = new_data

    y_ticks = np.linspace(y_min, y_max, 5)
    y_tick_text = [f"{tick:.3f}" for tick in y_ticks]

    colors_rgba = methods.generate_colors("rgba", len(frags_col))

    if region:
        for ii_f, frag in enumerate(frags_col):
            fig = go.Figure()
            probe = probes_to_frag.get(frag, "")
            output_path = os.path.join(
                output_dir,
                f"{sample_name}_{frag}_{probe}_{profile_type}_{bin_suffix}_{log_status}{region}.{extension}",
            )

            fig.add_trace(
                go.Scattergl(
                    x=df[x_col],
                    y=df[frag],
                    name=frag,
                    mode="lines",
                    line=dict(width=1, color=colors_rgba[ii_f]),
                    marker=dict(size=4),
                    showlegend=False,
                )
            )

            fig.update_layout(
                title=f"{sample_name}",
                xaxis=dict(
                    title=dict(text=f"{region} coordinates", standoff=80),
                    tickformat="d",
                    range=[x_min, x_max],
                    showgrid=False,
                ),
                yaxis=dict(
                    title=f"{profile_type.capitalize()}{' - log' if log_scale else ''}",
                    tickvals=y_ticks,
                    ticktext=y_tick_text,
                    range=[y_min, y_max],
                    showgrid=False,
                ),
                xaxis_type="linear",
                xaxis_tickformat="d",
                yaxis_tickformat="%.4e" if log_scale else "%.4d",
                plot_bgcolor="white",
                paper_bgcolor="white",
                width=width,
                height=height,
            )

            pio.write_image(fig, output_path, engine="kaleido")

    else:
        df_10kb_tmp = graph.build_bins_template(df_coords=df_coords, bin_size=10000)
        colorbar, chr_ticks_pos = graph.colorbar_maker(df_10kb_tmp)

        # plot for each prob or group of probes
        for ii_f, frag in enumerate(frags_col):
            probe = probes_to_frag.get(frag, "")
            output_path = os.path.join(
                output_dir,
                f"{sample_name}_{frag}_{probe}_{profile_type}_{bin_suffix}_{log_status}.{extension}",
            )

            fig = make_subplots(
                rows=2,
                cols=1,
                row_heights=[0.94, 0.06],
                vertical_spacing=0.06,
                specs=[[{"type": "scatter"}], [{"type": "bar"}]],
            )

            fig.add_trace(
                go.Scatter(
                    x=df[x_col],
                    y=df[frag],
                    mode="lines",
                    line=dict(width=1, color=colors_rgba[ii_f]),
                    marker=dict(size=4),
                    showlegend=False,
                ),
                row=1,
                col=1,
            )

            fig.add_trace(colorbar, row=2, col=1)
            title = f" {sample_name} <br> {frag} {'- ' + probe}"
            fig.update_layout(
                title=title,
                xaxis=dict(
                    title=dict(text="Genomic coordinates", standoff=80),
                    tickformat="d",
                    range=[x_min, x_max],
                    showgrid=False,
                ),
                xaxis2=dict(
                    tickmode="array",
                    tickvals=chr_ticks_pos,
                    ticktext=df["chr"].unique(),
                    tickfont=dict(size=12),
                ),
                yaxis=dict(
                    title=f"{profile_type.capitalize()}{' - log' if log_scale else ''}",
                    tickvals=y_ticks,
                    ticktext=y_tick_text,
                    range=[y_min, y_max],
                    showgrid=False,
                ),
                yaxis2=dict(
                    showticklabels=False,
                ),
                xaxis_showgrid=False,
                yaxis_showgrid=False,
                xaxis_type="linear",
                xaxis_tickformat="d",
                yaxis_tickformat="%.4e" if log_scale else "%.4d",
                xaxis_range=[x_min, x_max],
                yaxis_range=[y_min, y_max],
                hovermode="closest",
                plot_bgcolor="white",
                paper_bgcolor="white",
                width=width,
                height=height,
            )

            pio.write_image(fig, output_path, engine="kaleido")


def profile_contacts(
    filtered_table_path: str,
    oligo_capture_with_frag_path: str,
    chromosomes_coord_path: str,
    normalize: bool = False,
    output_path: str = None,
    additional_groups_path: str = None,
    force: bool = False,
):
    """
    Organize the contacts made by each probe with the genome and save the results as two .tsv files:
    one for contacts and one for frequencies.

    Parameters
    ----------
    filtered_table_path : str
        Path to the filtered table (sshictuff filter script output).
    oligo_capture_with_frag_path : str
        Path to the oligo capture file (table .csv or .tsv for oligo capture information).
        Must be the file with the fragments associated made with the 'associate' command.
    chromosomes_coord_path : str
        Path to the chromosomes coordinates file containing the length of each chromosome arms.
    normalize : bool
        Normalize the contacts by the total number of contacts.
    output_path : str
        Path to the output directory.
    additional_groups_path : str
        Path to the additional groups file (table .csv or .tsv for additional groups information).
    force : bool
        Force the overwriting of the output file if the file exists.
    """

    methods.check_if_exists(filtered_table_path)
    methods.check_if_exists(oligo_capture_with_frag_path)
    methods.check_if_exists(chromosomes_coord_path)

    if not output_path:
        output_path = filtered_table_path.replace(
            "filtered.tsv", "0kb_profile_contacts.tsv"
        )

    basedir = os.path.dirname(output_path)
    if not os.path.exists(basedir):
        os.makedirs(basedir)

    if os.path.exists(output_path) and not force:
        logger.warning("Output file already exists: %s", output_path)
        logger.warning("Use the --force / -F flag to overwrite the existing file.")
        return

    chr_coord_delim = "," if chromosomes_coord_path.endswith(".csv") else "\t"
    df_coords: pd.DataFrame = pd.read_csv(
        chromosomes_coord_path, sep=chr_coord_delim, index_col=None
    )
    df_chr_len = df_coords[["chr", "length"]]
    chr_list = list(df_chr_len["chr"].unique())
    df_chr_len["chr_start"] = df_chr_len["length"].shift().fillna(0).astype("int64")
    df_chr_len["cumu_start"] = df_chr_len["chr_start"].cumsum()

    oligo_delim = "," if oligo_capture_with_frag_path.endswith(".csv") else "\t"
    df_oligo: pd.DataFrame = pd.read_csv(oligo_capture_with_frag_path, sep=oligo_delim)
    probes = df_oligo["name"].to_list()
    fragments = df_oligo["fragment"].astype(str).to_list()

    df: pd.DataFrame = pd.read_csv(filtered_table_path, sep="\t")
    df_contacts: pd.DataFrame = pd.DataFrame(columns=["chr", "start", "sizes"])
    df_contacts: pd.DataFrame = df_contacts.astype(
        dtype={"chr": str, "start": int, "sizes": int}
    )

    for x in ["a", "b"]:
        y = methods.frag2(x)
        df2 = df[~pd.isna(df["name_" + x])]

        for probe in probes:
            if probe not in pd.unique(df2["name_" + x]):
                tmp = pd.DataFrame(
                    {
                        "chr": [np.nan],
                        "start": [np.nan],
                        "sizes": [np.nan],
                        probe: [np.nan],
                    }
                )

            else:
                df3 = df2[df2["name_" + x] == probe]
                tmp = pd.DataFrame(
                    {
                        "chr": df3["chr_" + y],
                        "start": df3["start_" + y],
                        "sizes": df3["size_" + y],
                        probe: df3["contacts"],
                    }
                )

            df_contacts = pd.concat([df_contacts, tmp])

    group = df_contacts.groupby(by=["chr", "start", "sizes"], as_index=False)
    df_contacts: pd.DataFrame = group.sum()
    df_contacts = methods.sort_by_chr(df_contacts, chr_list, "chr", "start")
    df_contacts.index = range(len(df_contacts))

    for probe, frag in zip(probes, fragments):
        df_contacts.rename(columns={probe: frag}, inplace=True)

    df_contacts: pd.DataFrame = df_contacts.loc[:, ~df_contacts.columns.duplicated()]

    df_merged: pd.DataFrame = df_contacts.merge(df_chr_len, on="chr")
    df_merged["genome_start"] = df_merged["cumu_start"] + df_merged["start"]
    df_contacts.insert(3, "genome_start", df_merged["genome_start"])

    if normalize:
        df_frequencies = df_contacts.copy(deep=True)
        for frag in fragments:
            frag_sum = df_frequencies[frag].sum()
            if frag_sum > 0:
                df_frequencies[frag] /= frag_sum

    if additional_groups_path:
        df_additional: pd.DataFrame = pd.read_csv(additional_groups_path, sep="\t")
        probes_to_fragments = dict(zip(probes, fragments))
        methods.make_groups_of_probes(df_additional, df_contacts, probes_to_fragments)
        if normalize:
            methods.make_groups_of_probes(
                df_additional, df_frequencies, probes_to_fragments
            )

    df_contacts.to_csv(output_path, sep="\t", index=False)
    if normalize:
        df_frequencies.to_csv(
            output_path.replace("contacts", "frequencies"), sep="\t", index=False
        )


def profile_probes_only(
    filtered_table_path: str,
    oligo_capture_with_frag_path: str,
    output_path: str = None,
    force: bool = False,
):

    methods.check_if_exists(filtered_table_path)
    methods.check_if_exists(oligo_capture_with_frag_path)

    if not output_path:
        output_path = filtered_table_path.replace(
            "filtered.tsv", "probes_vs_probes_profile.tsv"
        )

    basedir = os.path.dirname(output_path)
    if not os.path.exists(basedir):
        os.makedirs(basedir)

    if os.path.exists(output_path) and not force:
        logger.warning("Output file already exists: %s", output_path)
        logger.warning("Use the --force / -F flag to overwrite the existing file.")
        return

    oligo_delim = "," if oligo_capture_with_frag_path.endswith(".csv") else "\t"
    df_oligo: pd.DataFrame = pd.read_csv(oligo_capture_with_frag_path, sep=oligo_delim)
    probes = df_oligo["name"].to_list()
    fragments = df_oligo["fragment"].astype(int).to_list()

    df: pd.DataFrame = pd.read_csv(filtered_table_path, sep="\t")
    df2 = df[df["frag_a"].isin(fragments) & df["frag_b"].isin(fragments)]

    square_matrix = np.zeros((len(probes), len(probes)))
    for i, probe_a in enumerate(probes):
        frag_a = df_oligo.loc[df_oligo["name"] == probe_a, "fragment"].values[0]
        for j, probe_b in enumerate(probes):
            frag_b = df_oligo.loc[df_oligo["name"] == probe_b, "fragment"].values[0]
            contacts = df2[(df2["frag_a"] == frag_a) & (df2["frag_b"] == frag_b)][
                "contacts"
            ].sum()
            square_matrix[i, j] = contacts

    # normalize the matrix
    total_contacts = square_matrix.sum()
    if total_contacts > 0:
        square_matrix /= total_contacts

    # make the matrix symmetric
    square_matrix = np.maximum(square_matrix, square_matrix.T)
    df_contacts = pd.DataFrame(square_matrix, columns=probes, index=probes)
    df_contacts.to_csv(output_path, sep="\t", index=True)


def rebin_profile(
    contacts_unbinned_path: str,
    chromosomes_coord_path: str,
    bin_size: int,
    output_path: str = None,
    force: bool = False,
) -> None:
    """
    Rebin the contacts from the unbinned contacts file to the binned contacts file.

    Parameters
    ----------
    contacts_unbinned_path : str
        Path to the unbinned contacts file.
    chromosomes_coord_path : str
        Path to the chromosomes coordinates file containing the length of each chromosome arms.
    bin_size : int
        Size of the bins (resolution).
    output_path : str
        Path to the output file.
    force : bool
        Force the overwriting of the output file if the file exists.
    """

    methods.check_if_exists(contacts_unbinned_path)
    methods.check_if_exists(chromosomes_coord_path)

    bin_suffix = methods.get_bin_suffix(bin_size)
    if not output_path:
        output_path = contacts_unbinned_path.replace(
            "0kb_profile", f"{bin_suffix}_profile"
        )

    if os.path.exists(output_path) and not force:
        logger.warning("[Rebin] : Output file already exists: %s", output_path)
        logger.warning(
            "[Rebin] : Use the --force / -F flag to overwrite the existing file."
        )
        return

    df = pd.read_csv(contacts_unbinned_path, sep="\t")
    coord_delim = "," if chromosomes_coord_path.endswith(".csv") else "\t"
    df_coords: pd.DataFrame = pd.read_csv(
        chromosomes_coord_path, sep=coord_delim, index_col=None
    )

    chr_sizes = dict(zip(df_coords.chr, df_coords.length))
    chr_list, chr_bins = [], []

    for c, l in chr_sizes.items():
        chr_list.append([c] * (l // bin_size + 1))
        chr_bins.append(np.arange(0, (l // bin_size + 1) * bin_size, bin_size))

    chr_list = np.concatenate(chr_list)
    chr_bins = np.concatenate(chr_bins)

    df_template = pd.DataFrame(
        {
            "chr": chr_list,
            "chr_bins": chr_bins,
            "genome_bins": np.arange(0, len(chr_bins) * bin_size, bin_size),
        }
    )

    df["end"] = df["start"] + df["sizes"]
    df["start_bin"] = df["start"] // bin_size * bin_size
    df["end_bin"] = df["end"] // bin_size * bin_size
    df.drop(columns=["genome_start"], inplace=True)

    df_cross_bins = df[df["start_bin"] != df["end_bin"]].copy()
    df_in_bin = df.drop(df_cross_bins.index)
    df_in_bin["chr_bins"] = df_in_bin["start_bin"]

    df_cross_bins_a = df_cross_bins.copy()
    df_cross_bins_b = df_cross_bins.copy()
    df_cross_bins_a["chr_bins"] = df_cross_bins["start_bin"]
    df_cross_bins_b["chr_bins"] = df_cross_bins["end_bin"]

    fragments_columns = df.filter(regex=r"^\d+$|^\$").columns.to_list()

    correction_factors = (
        df_cross_bins_b["end"] - df_cross_bins_b["chr_bins"]
    ) / df_cross_bins_b["sizes"]
    for c in fragments_columns:
        df_cross_bins_a[c] *= 1 - correction_factors
        df_cross_bins_b[c] *= correction_factors

    df_binned = pd.concat([df_cross_bins_a, df_cross_bins_b, df_in_bin])
    df_binned.drop(columns=["start_bin", "end_bin"], inplace=True)

    df_binned = df_binned.groupby(["chr", "chr_bins"]).sum().reset_index()
    df_binned = methods.sort_by_chr(df_binned, chr_list, "chr_bins")
    df_binned = pd.merge(df_template, df_binned, on=["chr", "chr_bins"], how="left")
    df_binned.drop(columns=["start", "end", "sizes"], inplace=True)
    df_binned.fillna(0, inplace=True)

    df_binned.to_csv(output_path, sep="\t", index=False)
