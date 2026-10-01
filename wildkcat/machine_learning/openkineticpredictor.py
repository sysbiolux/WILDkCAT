import re, os, sys, time
import logging
import requests
import itertools 
import pandas as pd 
from tqdm import tqdm
from dotenv import load_dotenv
from functools import lru_cache

from wildkcat.api.api_utilities import safe_requests_get, retry_api
from wildkcat.api.uniprot_api import convert_uniprot_to_sequence, identify_catalytic_enzyme


load_dotenv()


# --- API ---


def convert_kegg_compound_to_sid(kegg_compound_id) -> str | None:
    """
    Convert the KEGG compound ID to the PubChem Substance ID (SID).

    Parameters:
        kegg_compound_id (str): KEGG compound ID.

    Returns:
        str: The PubChem SID if found, otherwise None.
    """
    url = f"https://rest.kegg.jp/conv/pubchem/compound:{kegg_compound_id}"
    safe_get_with_retry = retry_api()(safe_requests_get)
    response = safe_get_with_retry(url)

    if response is None:
        return None

    if response.status_code != 200:
        return None

    match = re.search(r'pubchem:\s*(\d+)', response.text)
    sid = match.group(1) if match else None
    return sid


def convert_sid_to_cid(sid) -> int | None:
    """
    Converts a PubChem Substance ID (SID) to the corresponding Compound ID (CID).

    Parameters:
        sid (str): PubChem Substance ID.

    Returns:
        int or None: The corresponding PubChem Compound ID (CID), or None if not found.
    """
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/substance/sid/{sid}/cids/JSON"
    safe_get_with_retry = retry_api()(safe_requests_get)
    response = safe_get_with_retry(url)
    
    if response is None:
        return None

    if response.status_code == 200:
        try:
            cid = response.json()['InformationList']['Information'][0]['CID'][0]
        except (KeyError, IndexError):
            cid = None
    return cid


def convert_cid_to_smiles(cid) -> list | None:    
    """
    Converts a PubChem Compound ID (CID) to its corresponding SMILES representation.

    Parameters:
        cid (str): PubChem Compound ID.

    Returns:
       list or None: A list of SMILES strings if found, otherwise None.
    """
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/smiles/txt"
    try:
        safe_get_with_retry = retry_api()(safe_requests_get)
        response = safe_get_with_retry(url)

        if response is None:
            return None

        response.raise_for_status()
        smiles = response.text.strip().split('\n')
        return smiles
    except:
        return None


@lru_cache(maxsize=None)
def convert_kegg_to_smiles(kegg_compound_id) -> list | None:
    """
    Convert the KEGG compound ID to the PubChem Compound ID (CID).

    Parameters:
        kegg_compound_id (str): KEGG compound ID.

    Returns:
        list or None: A list of SMILES strings if found, otherwise None.
    """
    sid = convert_kegg_compound_to_sid(kegg_compound_id)
    if sid is None:
        logging.warning('%s: Failed to retrieve SID for KEGG compound ID' % (kegg_compound_id))
        return None
    cid = convert_sid_to_cid(sid)
    if cid is None:
        logging.warning('%s: Failed to retrieve CID for KEGG compound ID' % (kegg_compound_id))
        return None
    smiles = convert_cid_to_smiles(cid)
    if smiles is None:
        logging.warning('%s: Failed to retrieve SMILES for KEGG compound ID' % (kegg_compound_id))
        return None
    return smiles


# --- Create Input file ---


def create_prediction_input_file(kcat_df, method):
    """
    Generate input file for OpenKineticsPredictor and a mapping of substrate KEGG IDs to SMILES.

    Parameters: 
        kcat_df (pd.DataFrame): Input DataFrame containing kcat information.

    Returns:
        input_df (pd.DataFrame): DataFrame for OpenKineticsPredictor input.
        substrates_to_smiles (dict): Mapping KEGG ID <-> SMILES.
        uniprot_to_sequence (dict): Mapping UniProt ID <-> Protein Sequence.
    """
    
    def resolve_smiles(kegg_ids, names, cofactor):
        smiles_out = []
        for name, kid in zip(names, kegg_ids):
            if not kid:
                continue
            if name.lower() in [c.lower() for c in cofactor]:
                nonlocal counter_cofactor
                counter_cofactor += 1  # TODO: Should we add a warning if no cofactor is found for a reaction? 
                continue
            smiles = convert_kegg_to_smiles(kid)
            if smiles is not None:
                s = smiles[0]  # If multiple SMILES, take the first one 
                smiles_out.append(s)
                substrates_to_smiles[kid] = s
        return smiles_out
    
    single_substrate_methods   = {"CataPro", "DLKcat", "EITLEM", "KinForm-H", "KinForm-L", 
                                  "MMISA-KM", "OmniESI", "UniKP"}
    multiple_substrate_methods = {"CatPred"}
    full_reactions_methods     = {"TurNup"}

    is_single = method in single_substrate_methods
    is_multi = method in multiple_substrate_methods
    is_full = method in full_reactions_methods

    content = []
    substrates_to_smiles = {}
    uniprot_to_sequence = {}

    counter_no_catalytic, counter_kegg_no_matching, counter_rxn_covered, counter_cofactor = 0, 0, 0, 0

    for _, row in tqdm(kcat_df.iterrows(), total=len(kcat_df), desc="Generating OpenKineticsPredictor input"):

        uniprot = row['uniprot']
        ec_code = row['ec_code']

        if len(uniprot.split(';')) > 1:       
            catalytic_enzyme = identify_catalytic_enzyme(uniprot, ec_code)
            if catalytic_enzyme is None or (";" in str(catalytic_enzyme)):
                counter_no_catalytic += 1
                continue
            else: 
                uniprot = catalytic_enzyme

        # If the number of KEGG Compound IDs is not matching the number of names  
        kegg_ids = [s for s in row["substrates_kegg"].split(";") if s]
        names    = row["substrates_name"].split(";")
        if len(kegg_ids) != len(names):
            logging.warning(
                f"Number of KEGG compounds IDs does not match number of names for {ec_code}: {uniprot}."
            )
            counter_kegg_no_matching += 1

        # Get sequence 
        sequence = convert_uniprot_to_sequence(uniprot) 
        if sequence is None:
            continue
        uniprot_to_sequence[uniprot] = sequence

        smiles_list = []
        names = row['substrates_name'].split(';')
        kegg_ids = row['substrates_kegg'].split(';')
        
        # Get the cofactor for the EC code
        # cofactor = get_cofactor(ec_code) 
        cofactor = [] # Rm cofactor
        # Get SMILES
        smiles_list  = resolve_smiles(kegg_ids, names, cofactor)

        if not smiles_list:
            continue
        
        if is_single:
            
            for smiles_str in smiles_list:
                content.append({
                    "Protein Sequence": sequence,
                    "Substrate": smiles_str,
                })
            counter_rxn_covered += 1

        elif is_multi:
            content.append({
                "Protein Sequence": sequence,
                "Substrate": ".".join(smiles_list),
            })
            counter_rxn_covered += 1
        
        elif is_full:
            substrates_str = ";".join(smiles_list)
            products_str = ""
            has_products_col = is_full and "products_kegg" in kcat_df.columns
            if has_products_col:
                product_kegg_ids = [p for p in str(row["products_kegg"]).split(";") if p]
            product_smiles   = []
            for kid in product_kegg_ids:
                smiles = convert_kegg_to_smiles(kid)
                if smiles is not None:
                    s = smiles[0]
                    product_smiles.append(s)
                    substrates_to_smiles[kid] = s   # products also cached
                    products_str = ";".join(product_smiles)

            content.append({
                "Protein Sequence": sequence,
                "Substrates": substrates_str,
                "Products": products_str
                })
                
            counter_rxn_covered += 1


    # Generate OpenKineticsPredictor input file
    input_df = pd.DataFrame(content)
    # Remove duplicates
    before_duplicates_filter = len(input_df)
    input_df = input_df.drop_duplicates().reset_index(drop=True)
    nb_lines_dropped = before_duplicates_filter - len(input_df)
    
    # Remove 'nan' values
    seq_col = "Protein Sequence"
    sub_col = "Substrates" if is_full else "Substrate"
    input_df = input_df.dropna(subset=[seq_col, sub_col])
    input_df = input_df[(input_df[seq_col].str.strip() != "") & (input_df[sub_col].str.strip() != "")].reset_index(drop=True)

    # Generate reverse mapping as TSV
    substrates_to_smiles_df = pd.DataFrame(list(substrates_to_smiles.items()), columns=['KEGGID', 'smiles'])
    uniprot_to_sequence_df = pd.DataFrame(list(uniprot_to_sequence.items()), columns=['UniProtID', 'sequence'])

    report_statistics = {
        "rxn_covered": counter_rxn_covered,
        "cofactor_identified": counter_cofactor,
        "no_catalytic": counter_no_catalytic,
        "kegg_no_matching": counter_kegg_no_matching,
        "duplicates_enzyme_substrates": nb_lines_dropped,
    }

    return input_df, substrates_to_smiles_df, uniprot_to_sequence_df, report_statistics


def send_kinetic_predictor(output_folder, method, canonicalize_substrates):
    """
    TODO: Write the docsting 
    """
    api_key = os.getenv("KP_API_KEY")
    HEADERS = {"Authorization": f"Bearer {api_key}"}

    with open(os.path.join(output_folder, "machine_learning/prediction_input.csv"), "rb") as f:
        resp = requests.post(
            "https://predictor.openkinetics.org/api/v1/submit/",
            headers=HEADERS,
            files={"file": f},
            data={
                "targets": '["kcat"]', #TODO: Add the possibility to choose between kcat and km prediction
                "methods": f'{{"kcat":"{method}"}}',
                "handleLongSequences": "truncate",
                "canonicalizeSubstrates": f'"{canonicalize_substrates}"',
            },
        )

    # TODO: Handle erros 
    # print(resp.json())

    job_id = resp.json().get("jobId")
    status = "pending"

    print(f"Prediction {job_id} in progress ", end="", flush=True)
    spinner = itertools.cycle(["|", "/", "-", "\\"])

    while status != "Completed":

        sys.stdout.write(next(spinner))
        sys.stdout.flush()
        time.sleep(0.2)
        sys.stdout.write("\b")

        resp = requests.get(
            f"https://predictor.openkinetics.org/api/v1/status/{job_id}/",
            headers=HEADERS,
        )
        status = resp.json().get("status")

    print(f"\r{next(spinner)} Running prediction {job_id})")
    print(f"Prediction {job_id} completed")

    resp_json = requests.get(
        f"https://predictor.openkinetics.org/api/v1/result/{job_id}/?format=json",
        headers=HEADERS,
    )

    prediction_output = resp_json.json()
    # with open("prediction_output_unformated.json", "w", encoding="utf-8") as f: # TODO: Rm 
    #     json.dump(prediction_output, f, ensure_ascii=False, indent=4)
    # prediction_output_formatted = format_prediction_output(prediction_output, output_folder) #TODO 
    prediction_output_formatted = format_prediction_output(prediction_output, output_folder, method) #TODO  

    return prediction_output_formatted


def format_prediction_output(result_json, output_folder, method):
    """
    Format the raw prediction output from the OpenKineticsPredictor API into a structured DataFrame,
    mapping sequences back to UniProt IDs and SMILES back to KEGG IDs.

    Handles three input/output formats depending on the method:
        - Single substrate: one SMILES per prediction row (column "Substrate").
        - Multiple substrates: substrates joined by '.' in a single "Substrate" column.
        - Full reaction: substrates and products as separate columns ("Substrates", "Products"),
          each containing SMILES separated by ';'.

    Parameters:
        result_json (dict): Raw JSON response from the prediction API.
        output_folder (str): Path to the output folder containing auxiliary mapping files.
        method (str): Prediction method used; determines the expected column format.

    Returns:
        pd.DataFrame: Formatted DataFrame with columns ['UniProtID', 'KEGGID', 'Kcat'].
    """
    single_substrate_methods   = {"CataPro", "DLKcat", "EITLEM", "KinForm-H", "KinForm-L",
                                  "MMISA-KM", "OmniESI", "UniKP"}
    multiple_substrate_methods = {"CatPred"}
    full_reactions_methods     = {"TurNup"}

    is_single = method in single_substrate_methods
    is_multi  = method in multiple_substrate_methods
    is_full   = method in full_reactions_methods

    output_path = os.path.join(output_folder, "machine_learning/prediction_")
    uniprot_to_sequence_df  = pd.read_csv((output_path + 'input_uniprot_to_sequence.tsv'),  sep='\t')
    substrates_to_smiles_df = pd.read_csv((output_path + 'input_substrates_to_smiles.tsv'), sep='\t')

    result_df = pd.DataFrame(result_json["data"])
    result_df = result_df[result_df["kcat (1/s)"].notna()]  # Remove invalid 

    if is_full:
        result_df = result_df.rename(columns={
            "kcat (1/s)":       "Kcat",
            "Protein Sequence": "sequence",
            "Substrates":       "substrates_smiles",
            "Products":         "products_smiles",
        })
    else:
        result_df = result_df.rename(columns={
            "kcat (1/s)":       "Kcat",
            "Protein Sequence": "sequence",
            "Substrate":        "substrate_smiles",
        })

    smiles_to_kegg = (
        substrates_to_smiles_df
        .groupby('smiles')['KEGGID']
        .apply(lambda x: ';'.join(sorted(set(x.dropna().astype(str)))))
        .to_dict()
    )

    # Helpers
    def resolve_kegg_ids_single(smiles_str):
        return smiles_to_kegg.get(smiles_str)

    def resolve_kegg_ids_multi(smiles_dot_joined):

        individual = smiles_dot_joined.split('.')
        kegg_ids = [smiles_to_kegg[s] for s in individual if s in smiles_to_kegg]
        return ';'.join(kegg_ids) if kegg_ids else None

    def resolve_kegg_ids_full(substrates_semicolon_joined):

        individual = substrates_semicolon_joined.split(';')
        kegg_ids = [smiles_to_kegg[s] for s in individual if s in smiles_to_kegg]
        return ';'.join(kegg_ids) if kegg_ids else None

    if is_single:
        result_df['KEGGID'] = result_df['substrate_smiles'].apply(resolve_kegg_ids_single)
    elif is_multi:
        result_df['KEGGID'] = result_df['substrate_smiles'].apply(resolve_kegg_ids_multi)
    elif is_full:
        result_df['KEGGID_Substrates'] = result_df['substrates_smiles'].apply(resolve_kegg_ids_full)
        result_df['KEGGID_Products'] = result_df['products_smiles'].apply(resolve_kegg_ids_full)

    result_df = pd.merge(result_df, uniprot_to_sequence_df, on='sequence', how='left')

    if is_full: 
        output_df = result_df[['UniProtID', 'KEGGID_Substrates', 'KEGGID_Products', 'Kcat']].copy()
    else: 
        output_df = result_df[['UniProtID', 'KEGGID', 'Kcat']].copy()

    output_df.to_csv((output_path + 'output_' + method + '.tsv'), sep='\t', index=False)
    return output_df

# --- Integrate predictions into kcat file ---


def integrate_predictions(kcat_df, predictions_df, method) -> pd.DataFrame:
    """
    Integrates predictions into an kcat file.
    If multiple values are provided for a single combination of EC, Enzyme, Substrate, the minimum value is taken.

    Parameters:
        kcat_df (pd.DataFrame): Input DataFrame containing kcat information.
        predictions_df (pd.DataFrame): DataFrame containing model predictions
        method (string): TODO

    Returns:
        pd.DataFrame: The input kcat_df with an additional column 'predicted_kcat_s' containing
            the integrated predicted kcat (s^-1) values.
    """
    single_substrate_methods   = {"CataPro", "DLKcat", "EITLEM", "KinForm-H", "KinForm-L",
                                  "MMISA-KM", "OmniESI", "UniKP"}
    multiple_substrate_methods = {"CatPred"}
    full_reactions_methods     = {"TurNup"}

    is_single = method in single_substrate_methods
    is_multi  = method in multiple_substrate_methods
    is_full   = method in full_reactions_methods

    def get_full_pred_kcat(row):
        key = (row['uniprot'], row['substrates_kegg'], row['products_kegg'])
        values = prediction_map.get(key, [])
        return min(values) if values else None

    def get_multi_pred_kcat(row):
        key = (row['uniprot'], row['substrates_kegg'])
        values = prediction_map.get(key, [])
        return min(values) if values else None
    
    def get_min_pred_kcat(row):
        """
        For a given kcat_df row, retrieve all predicted kcat values matching its
        (uniprot, kegg_id) pairs and return the minimum — or None if no match exists.
        """
        uniprot  = row['uniprot']
        kegg_ids = str(row['substrates_kegg']).split(';')
        kcat_values = [
            prediction_map.get((uniprot, kegg_id.strip()))
            for kegg_id in kegg_ids
            if (uniprot, kegg_id.strip()) in prediction_map
        ]
        return min(kcat_values) if kcat_values else None  # If multiple substrates, take the minimum kcat value

    if is_full: 
        prediction_map = (
            predictions_df
            .groupby(['UniProtID', 'KEGGID_Substrates', 'KEGGID_Products'])['Kcat']
            .apply(list)
            .to_dict()
        )
        kcat_df['predicted_kcat_s'] = kcat_df.apply(get_full_pred_kcat, axis=1)

    elif is_multi: 
        prediction_map = (
            predictions_df
            .groupby(['UniProtID', 'KEGGID'])['Kcat']
            .apply(list)
            .to_dict()
        )
        kcat_df['predicted_kcat_s'] = kcat_df.apply(get_multi_pred_kcat, axis=1)

    elif is_single: 
        prediction_map = predictions_df.set_index(['UniProtID', 'KEGGID'])['Kcat'].to_dict()
        kcat_df['predicted_kcat_s'] = kcat_df.apply(get_min_pred_kcat, axis=1)

    return kcat_df


# if __name__ == "__main__":
    # Test : Retrieve SMILES from KEGG ID
    # print(convert_kegg_to_smiles("C00008"))

    # Test : Retrieve Sequence from UniProt ID
    # print(convert_uniprot_to_sequence("P0A796"))

    # Test : Integrate CataPro predictions into kcat file
    # kcat_df = pd.read_csv("output/ecoli_kcat_sabio.tsv", sep='\t')
    # substrates_to_smiles = pd.read_csv('in_progress/ml_test/substrates_to_smiles.tsv', sep='\t')
    # integrate_catapro_predictions(kcat_df, substrates_to_smiles, "in_progress/ml_test/catapro_output.csv", "in_progress/ml_test/ecoli_kcat_catapro.tsv")