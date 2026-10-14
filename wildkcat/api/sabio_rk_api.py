import time, re
import requests
import logging
import pandas as pd 
from io import StringIO
from functools import lru_cache 

BASE = "https://sabiork.h-its.org/export-api/sabio/kinlaw-entry/json"
PAGE_SIZE = 1000
MIN_SECONDS_BETWEEN_REQUESTS = 1.1   # stays below 60 requests per minute

# --- Sabio-RK API ---


# Thank to Xiaoming Hu for the implementation 


def fetch_entries(q):
    """Yield every entry matching the Solr query q, following all pages."""
    page = 1
    while True:
        r = requests.post(BASE, json={"q": q, "page": page, "pageSize": PAGE_SIZE})
        if r.status_code == 429:          # rate limit hit: wait and retry
            time.sleep(5)
            continue
        r.raise_for_status()
        body = r.json()
        yield from body["data"]
        if page >= body["meta"]["total_pages"]:
            return
        page += 1
        time.sleep(MIN_SECONDS_BETWEEN_REQUESTS)


def kcat_rows(entry):
    """Yield one flat dict per kcat parameter of an entry, with KEGG IDs."""
    links = entry["external_links"]
    kegg = {l["id"]: l["value"] for l in links["compound"] if l["key"] == "KeggCompoundID"}
    kegg_rxn = next((l["value"] for l in links["reaction"] if l["key"] == "KeggReactionID"), None)

    def species(role):
        return [(s["compound"]["name"], kegg.get(s["compound"]["id"]))
                for s in entry["reaction"]["species"] if s["role"] == role]

    enz = entry["enzyme_description"]
    cond = entry["experimental_conditions"]
    base = {
        "EntryID": entry["id"],
        "ECNumber": enz["ec_number"],
        "KeggReactionID": kegg_rxn,
        "Reaction": entry["reaction"]["equation"],
        "Substrates": species("Substrate"),
        "Products": species("Product"),
        "UniProt": sorted({p["uniprot_id"] for p in enz["proteins"] if p["uniprot_id"]}),
        "Organism": entry["general"]["organism"]["name"],
        "EnzymeVariant": enz["wildtype"],
        # "MutantSpec": enz["mutant_spec"],
        "Temperature": cond["envvar_temperature"]["start_value"],
        "pH": cond["envvar_ph"]["start_value"],
        # "PubMedID": entry["publication"].get("pubmed_id"),
    }
    for p in entry["kineticlaw"]["parameter"]:
        if p["parameter_type"]["name"] != "kcat" or p["start_value"] is None:
            continue
        yield {**base, "Value": p["start_value"], "Unit": p["unit"]["name"]}


def send_sabio_query(q) -> pd.DataFrame:
    """Return all kcat records for the Solr query q."""
    df = pd.DataFrame([row for e in fetch_entries(q) for row in kcat_rows(e)]) 
    return df 


@lru_cache(maxsize=None)
def get_turnover_number_sabio(ec_number) -> pd.DataFrame:
    """
    Retrieve turnover number (kcat) data from SABIO-RK for a given EC number.

    Parameters:
        ec_number (str): Enzyme Commission number.

    Returns:
        pd.DataFrame: DataFrame containing kcat entries for a given EC.
    """
    query = f'ParameterType:"kcat" AND ECNumber:"{ec_number}"'

    df_entries = send_sabio_query(query)

    if df_entries.empty:
        logging.warning("%s: No data found for the query in SABIO-RK." % f"{ec_number}")
        return pd.DataFrame()  # Return empty DataFrame if no data found

    df_entries = format_sabio_output(df_entries)

    return df_entries


@lru_cache(maxsize=None)
def get_enzyme_sabio(uniprot_id) -> pd.DataFrame:
    """
    Retrieve enzyme data from SABIO-RK for a given UniProtKB accession.

    Parameters:
        uniprot_id (str): UniProtKB accession.
    
    Returns:
        pd.DataFrame: DataFrame containing kcat entries for a given UniProt.
    """
    query = f'ParameterType:"kcat" AND UniProtKB_AC:"{uniprot_id}"'

    df_entries = send_sabio_query(query)

    if df_entries.empty:
        logging.warning("%s: No data found for the query in SABIO-RK." % f"{uniprot_id}")
        return pd.DataFrame()  # Return empty DataFrame if no data found

    df_entries = format_sabio_output(df_entries)

    return df_entries


def simplify_uniprot(value) -> str:
    """
    Extract UniProtKB accessions from a raw SABIO-RK UniProt field.

    Parameters:
        value (str | list | tuple | Any): Raw field, e.g. ['(P08200)(O29610)'].
    Returns:
        str: ';'-separated unique accessions (order preserved), or '' if none found.
    """
    if isinstance(value, str):
        value = [value]
    elif not isinstance(value, (list, tuple)):
        return ""

    accessions = (
        ac
        for item in value if isinstance(item, str)
        for ac in re.compile(
            r"\b(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})\b")
            .findall(item)
    )
    return ";".join(dict.fromkeys(accessions))


def convert_kcat_to_per_second(df: pd.DataFrame) -> pd.DataFrame:
    """
    Convert kcat values to s^(-1).

    Parameters:
        df (pd.DataFrame): DataFrame with 'value' and 'unit' columns.
    Returns:
        pd.DataFrame: Copy of df with values converted to s^(-1). Entries with a missing value or an unsupported unit are set to NaN / None.
    """
    df = df.copy()
    df["value"] = pd.to_numeric(df["value"], errors="coerce")
    factors = df["unit"].map({
        "s^(-1)": 1,
        "min^(-1)": 1 / 60,
        "h^(-1)": 1 / 3600,
    })

    unsupported = df.loc[factors.isna() & df["value"].notna(), "unit"]
    for unit, n in unsupported.value_counts(dropna=False).items():
        logging.warning("Unsupported kcat unit %r (%d entries): value set to NaN.", unit, n)

    df["value"] = df["value"] * factors
    df["unit"] = "s^(-1)"
    df.loc[df["value"].isna(), "unit"] = None
    return df


def normalize_enzyme_variant(df: pd.DataFrame) -> pd.DataFrame:
    """
    Add an 'EnzymeVariant' column derived from 'MutantSpec'.

    Parameters:
        df (pd.DataFrame): DataFrame with a 'MutantSpec' column.
    Returns:
        pd.DataFrame: Copy of df with 'EnzymeVariant' set to 'wildtype', 'mutant' or None (unspecified).
    """
    df = df.copy()
    spec = df["EnzymeVariant"].fillna("").astype(str).str.strip().str.lower()

    is_missing = spec.isin({"", "nan", "none", "null", "n/a", "na", "-"})
    is_wildtype = spec.str.fullmatch(re.compile(r"wild[\s\-_]?type"))

    df["EnzymeVariant"] = "mutant"
    df.loc[is_wildtype, "EnzymeVariant"] = "wildtype"
    df.loc[is_missing, "EnzymeVariant"] = None
    return df


def _join_entries(entries, index: int) -> str:
    """
    Join one field of (name, KEGG ID) tuples into a ';'-separated string.

    Parameters:
        entries (list[tuple] | Any): List of (name, kegg_id) tuples; non-list values yield ''.
        index (int): Tuple position to extract (0: name, 1: KEGG ID).
    Returns:
        str: ';'-separated values, with '' for missing items.
    """
    if not isinstance(entries, list):
        return ""
    return ";".join("" if e[index] is None else str(e[index]) for e in entries)


def format_sabio_output(df: pd.DataFrame) -> pd.DataFrame: 
    """
    Format raw SABIO-RK entries (as returned by send_sabio_query) into the standard kcat output schema.

    Parameters:
        df (pd.DataFrame): Raw SABIO-RK entries.
    Returns:
        pd.DataFrame: SABIO-RK entries, with substrates/products flattened to ';'-separated strings and kcat values converted to s^(-1).
    """

    # print(df.columns)

    df = normalize_enzyme_variant(df)
    # Flatten substrates and products
    for src, dst in [("Substrates", "Substrate"), ("Products", "Product")]:
        df[dst] = df[src].apply(_join_entries, args=(0,))
        df[f"{dst}_KEGG"] = df[src].apply(_join_entries, args=(1,))

    df["UniProtKB_AC"] = df["UniProt"].apply(simplify_uniprot)

    df = df.rename(columns={"Value": "value", "Unit": "unit"})
    df["db"] = "sabio_rk"
    df = df[["ECNumber", "KeggReactionID", "Reaction", "Substrate", "Substrate_KEGG", "Product", "Product_KEGG", "UniProtKB_AC",
             "Organism", "EnzymeVariant", "Temperature", "pH", "value", "unit", "db"
             ]]

    return convert_kcat_to_per_second(df).reset_index(drop=True)
    
    
if __name__ == "__main__":
    # Test : Send a request to SABIO-RK API
    # entries = get_enzyme_sabio(uniprot_id="P00722")

    # for elm in entries: 
    #     print(elm)

    df = get_turnover_number_sabio(ec_number="1.1.1.42")
    df.to_csv("in_progress/sabio_rk_new.tsv", sep='\t', index=False)
    # print(df)
    # df.to_csv("in_progress/sabio_rk_test.tsv", sep='\t', index=False)