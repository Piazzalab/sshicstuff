import os
import pandas as pd
from sshicstuff.core import methods
from sshicstuff.log import logger

def filter_contacts(sparse_mat_path, oligo_capture_path, fragments_list_path, output_path=None, force=False):
    """
    Filter sparse Hi-C contact matrix to retain only contacts where at least one fragment has an associated oligo.

    Parameters:
    - sparse_mat_path (str): Path to the Hi-C sparse matrix file.
    - oligo_capture_path (str): Path to the oligo capture file.
    - fragments_list_path (str): Path to the fragment list file.
    - output_path (str, optional): Output path for the filtered matrix.
    - force (bool, optional): Overwrite existing output file if True.
    """
    if not output_path:
        output_path = sparse_mat_path.replace(".txt", "_filtered.tsv")

    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    if not force and os.path.exists(output_path):
        logger.warning("[Filter] : Output file already exists: %s", output_path)
        logger.warning("[Filter] : Use the --force / -F flag to overwrite.")
        return

    for path in [sparse_mat_path, oligo_capture_path, fragments_list_path]:
        methods.check_if_exists(path)

    methods.check_file_extension(sparse_mat_path, ".txt")
    methods.check_file_extension(oligo_capture_path, [".csv", ".tsv"])
    methods.check_file_extension(fragments_list_path, ".txt")

    df_fragments = load_fragments(fragments_list_path)
    df_oligos = load_oligos(oligo_capture_path)
    df_contacts = load_contacts(sparse_mat_path)

    df_oligo_fragments = associate_oligos_to_fragments(df_fragments, df_oligos)
    df_joined = pd.concat([
        join_contacts("a", df_contacts, df_fragments, df_oligo_fragments),
        join_contacts("b", df_contacts, df_fragments, df_oligo_fragments)
    ])

    df_joined.drop(columns=["frag"], inplace=True)
    df_joined.sort_values(by=["frag_a", "frag_b", "start_a", "start_b"], inplace=True)
    df_joined.reset_index(drop=True).to_csv(output_path, sep="\t", index=False)

    logger.info("[Filter] : Filtered contacts saved to %s", output_path)

def load_fragments(path):
    """
    Load and format the fragments file.

    Parameters:
    - path (str): Path to the fragments file.

    Returns:
    - pd.DataFrame: Processed fragments DataFrame.
    """
    df = pd.read_csv(path, sep="\t")
    return pd.DataFrame({
        "frag": range(len(df)),
        "chr": df["chrom"],
        "start": df["start_pos"],
        "end": df["end_pos"],
        "size": df["size"],
        "gc_content": df["gc_content"],
    })

def load_oligos(path):
    """
    Load and format the oligo capture file.

    Parameters:
    - path (str): Path to the oligo capture file.

    Returns:
    - pd.DataFrame: Processed oligos DataFrame.
    """
    sep = "," if path.endswith(".csv") else "\t"
    df = pd.read_csv(path, sep=sep)
    df.columns = [c.lower() for c in df.columns]
    df = df[["chr", "start", "end", "name", "type", "sequence"]]
    return df.sort_values(by=["chr", "start"]).reset_index(drop=True)

def load_contacts(path):
    """
    Load sparse matrix of fragment contacts.

    Parameters:
    - path (str): Path to the contact matrix.

    Returns:
    - pd.DataFrame: Processed contacts DataFrame.
    """
    df = pd.read_csv(path, sep="\t", header=None, skiprows=1)
    df.columns = ["frag_a", "frag_b", "contacts"]
    return df

def associate_oligos_to_fragments(fragments, oligos):
    """
    Assign each oligo to its corresponding fragment based on position.

    Parameters:
    - fragments (pd.DataFrame): Fragment annotations.
    - oligos (pd.DataFrame): Oligo annotations.

    Returns:
    - pd.DataFrame: Oligos with matched fragment information.
    """
    oligo_mid = ((oligos["end"] - oligos["start"] - 1) // 2 + oligos["start"] - 1).astype(int)
    starts = []

    for i, row in oligos.iterrows():
        frag_matches = fragments[
            (fragments["chr"] == row["chr"]) &
            (fragments["start"] <= oligo_mid[i]) &
            (fragments["end"] >= oligo_mid[i])
        ]
        start = frag_matches.iloc[-1]["start"] if row["chr"] == "chr_artificial" else frag_matches.iloc[0]["start"]
        starts.append(start)

    oligos["start"] = starts
    merged = fragments.merge(oligos.drop(columns="end"), on=["chr", "start"])
    return merged

def join_contacts(which, contacts, fragments, oligo_frags):
    """
    Join contacts with oligos on one side and fragment metadata on the other.

    Parameters:
    - which (str): Either "a" or "b", indicating which fragment to join with oligos.
    - contacts (pd.DataFrame): Hi-C contact pairs.
    - fragments (pd.DataFrame): Fragment metadata.
    - oligo_frags (pd.DataFrame): Oligos with fragment association.

    Returns:
    - pd.DataFrame: Enriched contact table.
    """
    frag_col = f"frag_{which}"
    y = methods.frag2(which)  # swap "a" <-> "b"

    joined = contacts.merge(oligo_frags, left_on=frag_col, right_on="frag")

    fragments_renamed = fragments.rename(columns={"frag": f"frag_{y}_meta"})
    joined = joined.merge(
        fragments_renamed,
        left_on=f"frag_{y}",
        right_on=f"frag_{y}_meta",
        suffixes=(f"_{which}", f"_{y}")
    )

    joined.drop(columns=[f"frag_{y}_meta"], inplace=True)

    for col in ["type", "name", "sequence"]:
        if col in joined.columns:
            joined.rename(columns={col: f"{col}_{which}"}, inplace=True)

    return joined
