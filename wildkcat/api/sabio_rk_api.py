import requests
import logging
import pandas as pd 
from io import StringIO
from functools import lru_cache 


# --- Sabio-RK API ---

# This should work in theory, but the API rate is problematic and isn't taken into account here. 
# Note: Requesting the kcat parameters using Sabio RK API doesn't just return the kcats (it also returns the KMs); furthermore, it also returns rows where the kcat is None. No comment. 


def send_sabio_query(query: str) -> pd.DataFrame:
    """
    Retrieve all kcat entries matching a query string from the SABIO-RK Export API. 

    Parameters:
        query (str): Query string (the `q` parameter).

    Returns:
        TODO
    """

    entries = []
    page = 1

    while True:
        params = {"q": query, "page": page, "pageSize": 1000}
        response = requests.get(
            "https://sabiork.h-its.org/export-api/sabio/kinlaw-entry/json",
            params=params,
        )
        response.raise_for_status()
        payload = response.json()

        entries.extend(payload.get("data", []))

        total_pages = payload.get("meta", {}).get("total_pages", page)
        if page >= total_pages:
            break
        page += 1

    return pd.DataFrame(entries)


@lru_cache(maxsize=None)
def get_turnover_number_sabio(ec_number) -> pd.DataFrame:
    """
    Retrieve turnover number (kcat) data from SABIO-RK for a given EC number.

    Parameters:
        ec_number (str): Enzyme Commission number.

    Returns:
        pd.DataFrame: DataFrame containing kcat entries for a given EC.
    """
    # -- Retrieve entryIDs --
    query = f'ParameterType:"kcat" AND ECNumber:"{ec_number}"'
    
    # Make GET request
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

    # -- Retrieve entryIDs --
    query = f'ParameterType:"kcat" AND UniProtKB_AC:"{uniprot_id}"'

    entries = send_sabio_query(query)
    if not entries:
        logging.warning('%s: No data found for the query in SABIO-RK.' % f"{uniprot_id}")
        return pd.DataFrame()  # Return empty DataFrame if no data found

    df_entries = format_sabio_output(df_entries)
    
    return df_entries


# old version
# def query_sabio(entryIDs) -> pd.DataFrame:
#     """
#     Retrieve SABIO-RK entries for given entry IDs.

#     Parameters:
#         entryIDs (list): List of SABIO-RK entry IDs.

#     Returns:
#         pd.DataFrame: DataFrame containing SABIO-RK entries.
#     """
#     parameters = 'https://sabiork.h-its.org/entry/exportToExcelCustomizable'

#     data_field = {'entryIDs[]': entryIDs}
#     # Possible fields to retrieve:
#     # EntryID, Reaction, Buffer, ECNumber, CellularLocation, UniProtKB_AC, Tissue, Enzyme Variant, Enzymename, Organism
#     # Temperature, pH, Activator, Cofactor, Inhibitor, KeggReactionID, KineticMechanismType, Other Modifier, Parameter,
#     # Pathway, Product, PubMedID, Publication, Rate Equation, SabioReactionID, Substrate
#     query = {'format':'tsv', 'fields[]':['EntryID', 'ECNumber', 'KeggReactionID', 'Reaction', 'Substrate', 'Product', 
#                                          'UniProtKB_AC', 'Organism', 'Enzyme Variant', 'Temperature', 'pH', 
#                                          'Parameter']}

#     # Make POST request
#     request = requests.post(parameters, params=query, data=data_field)
#     request.raise_for_status()

#     # Format the response into a DataFrame
#     df = pd.read_csv(StringIO(request.text), sep='\t')
#     df = df[df['parameter.name'].str.lower() == 'kcat'].reset_index(drop=True) # Keep only kcat parameters
#     # Convert Temperature and pH to numeric, coercing errors to NaN
#     df['Temperature'] = pd.to_numeric(df['Temperature'], errors='coerce')
#     df['pH'] = pd.to_numeric(df['pH'], errors='coerce')
#     # Drop unnecessary columns
#     df.drop(columns=['EntryID', 'parameter.name', 'parameter.type', 'parameter.associatedSpecies', 
#                      'parameter.endValue', 'parameter.standardDeviation'], inplace=True, errors='ignore')
#     # Drop duplicates based on normalized Substrate and Product sets
#     df["Substrate_set"] = df["Substrate"].fillna("").str.split(";").apply(lambda x: tuple(sorted(s.strip() for s in x if s.strip())))
#     df["Product_set"] = df["Product"].fillna("").str.split(";").apply(lambda x: tuple(sorted(s.strip() for s in x if s.strip())))
#     dedup_cols = [col for col in df.columns if col not in ["Substrate", "Product"]]
#     df = df.drop_duplicates(subset=dedup_cols + ["Substrate_set", "Product_set"], keep="first")
#     df = df.drop(columns=["Substrate_set", "Product_set"])
#     # Rename columns for consistency
#     df.rename(columns={
#         'ECNumber': 'ECNumber',
#         'KeggReactionID': 'KeggReactionID',
#         'Reaction': 'Reaction',
#         'Substrate': 'Substrate',
#         'Product': 'Product',
#         'UniProtKB_AC': 'UniProtKB_AC',
#         'Organism': 'Organism',
#         'Enzyme Variant': 'EnzymeVariant',
#         'Temperature': 'Temperature',
#         'pH': 'pH',
#         'parameter.startValue': 'value',
#         'parameter.unit': 'unit'
#     }, inplace=True)
#     # Add a column for the db
#     df['db'] = 'sabio_rk'
#     return df


def format_sabio_output(entries: pd.DataFrame) -> pd.DataFrame: 
    """
    Flatten raw SABIO-RK entries (as returned by send_sabio_query) into one
    row per kcat parameter.
    
    Parameters:
        entries (pd.DataFrame): Raw entries from send_sabio_query.
    
    Returns:
        pd.DataFrame: TODO
    """
    
    # keep only kcat parameters
    rows = []

    for entry in entries.to_dict("records"):
        general = entry.get("general") or {}
        reaction = entry.get("reaction") or {}
        kineticlaw = entry.get("kineticlaw") or {}
        enzyme_desc = entry.get("enzyme_description") or {}
        exp_cond = entry.get("experimental_conditions") or {}
        external_links = entry.get("external_links") or {}

        species = reaction.get("species") or []
        substrate = ";".join(
            s["compound"]["name"] for s in species if s.get("role") == "Substrate"
        )
        product = ";".join(
            s["compound"]["name"] for s in species if s.get("role") == "Product"
        )

        kegg_reaction_id = ";".join(dict.fromkeys(
            link["value"] for link in external_links.get("reaction") or []
            if link.get("key") == "KeggReactionID"
        ))
        uniprot_ac = ";".join(dict.fromkeys(
            link["value"] for link in external_links.get("kinlaw_entry") or []
            if link.get("key") == "UniProtKB_AC"
        ))

        enzyme_description = str(enzyme_desc.get("wildtype") or "").lower()
        if "wildtype" in enzyme_description:
            enzyme_variant = "wildtype"
        elif "mutant" in enzyme_description:
            enzyme_variant = "mutant"
        else:
            enzyme_variant = ""

        for param in kineticlaw.get("parameter") or []:

            if str(param.get("name", "")).lower() != "kcat":
                continue

            value = param.get("start_value")
            if value is None:
                continue  # skip kcat rows with no reported value

            unit = param.get("unit") or {}
            rows.append({
                "EntryID": entry.get("id"),
                "ECNumber": enzyme_desc.get("ec_number"),
                "KeggReactionID": kegg_reaction_id,
                "Reaction": reaction.get("equation"),
                "Substrate": substrate,
                "Product": product,
                "UniProtKB_AC": uniprot_ac,
                "Organism": (general.get("organism") or {}).get("name"),
                "EnzymeVariant": enzyme_variant,
                "Temperature": (exp_cond.get("envvar_temperature") or {}).get("start_value"),
                "pH": (exp_cond.get("envvar_ph") or {}).get("start_value"),
                "value_ori": value,
                "unit_ori": unit.get("name"),
                "value": param.get("n_start_value"),
                "unit": unit.get("n_name"),
            })

    df = pd.DataFrame(rows)
    if df.empty:
        return df

    # remove value_ori/unit_ori cols
    df = df.drop(["value_ori", "unit_ori"], axis=1)

    df["Temperature"] = pd.to_numeric(df["Temperature"], errors="coerce")
    df["pH"] = pd.to_numeric(df["pH"], errors="coerce")
    df["db"] = "sabio_rk"

    # EntryID excluded on purpose: different SABIO entries can report the
    # exact same kcat, and those should still collapse to one row.
    dedup_cols = [c for c in df.columns if c != "EntryID"]
    df = df.drop_duplicates(subset=dedup_cols, keep="first")

    # rm EntryID
    df = df.drop(columns=["EntryID"])
    return df

    
if __name__ == "__main__":
    # Test : Send a request to SABIO-RK API
    # entries = get_enzyme_sabio(uniprot_id="P00722")

    # for elm in entries: 
    #     print(elm)

    df = get_turnover_number_sabio(ec_number="1.1.1.42")
    df.to_csv("in_progress/sabio_rk_new.tsv", sep='\t', index=False)
    # print(df)
    # df.to_csv("in_progress/sabio_rk_test.tsv", sep='\t', index=False)