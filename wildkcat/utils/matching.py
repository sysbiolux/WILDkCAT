import re 
import logging
import pandas as pd 
import numpy as np
from typing import Optional, Tuple, Dict, Any

from ..utils.temperature import arrhenius_equation, calculate_ea
from ..utils.organism import closest_enz, closest_taxonomy


# TODO: Limit the Ea to the same pH ? 


# --- Main --- 


def find_best_match(kcat_dict, api_output, general_criteria) -> Tuple[float, Optional[Dict[str, Any]]]:
    """
    Finds the best matching enzyme entry from the provided API output based on: 
        - Kcat specific criteria: 
            * Substrate 
            * Catalytic enzyme(s)
        - General criteria : 
            * Organism
            * Temperature
            * pH

    This function filters out mutant enzyme variants, orders the remaining entries based on enzyme and organism similarity,
    and iteratively computes a score for each candidate to identify the best match. If a candidate requires an Arrhenius
    adjustment, the kcat value is recalculated accordingly.

    Parameters:
        kcat_dict (dict): Dictionary containing enzyme information.
        api_output (pd.DataFrame): DataFrame containing kcat entries and metadata from an API.
        general_criteria (dict): Dictionary specifying matching criteria.

    Returns:
        tuple:
            best_score (float): The lowest score found, representing the best match.
            best_candidate (dict or None): Dictionary of the best matching candidate's data, or None if no match is found.
    """

    # 1. Remove mutant enzymes
    api_output = api_output[api_output["EnzymeVariant"].isin(['wildtype', None])].copy()
    if api_output.empty:
        return 15, None

    # 2. Compute score and adjust kcat if needed
    scores = []
    adjusted_kcats, adjusted_temps, arrhenius_warning = [], [], []

    for _, row in api_output.iterrows():
        candidate_dict = row.to_dict()
        score, arrhenius = compute_score(kcat_dict, candidate_dict, general_criteria, api_output)
        candidate_dict['arrhenius'] = arrhenius
        if arrhenius != 0:
            kcat = arrhenius_equation(candidate_dict, api_output, general_criteria)
            if 10e-8 < kcat < 10e+8: 
                candidate_dict['value'] = kcat
                candidate_dict['Temperature'] = np.mean(general_criteria["Temperature"])
            # If the kcat value calculated using the Arrhenius is aberrant, use the non correct value instead
            else:
                logging.warning(f"{candidate_dict.get('ECNumber')}: Corrected kcat ({kcat:.0f} s-1) is outside the expected range of 10e-8, 10e+8.")
                temperature = float(candidate_dict['Temperature'])
                candidate_dict['arrhenius'] = 0
                if np.isnan(temperature): 
                    score += 1
                else: 
                    score += 2
        scores.append(score)
        adjusted_kcats.append(candidate_dict.get('value', row['value']))
        adjusted_temps.append(candidate_dict.get('Temperature', row['Temperature']))
        arrhenius_warning.append(candidate_dict.get('arrhenius'))

    api_output['score'] = scores
    api_output['adj_kcat'] = adjusted_kcats
    api_output['adj_temp'] = adjusted_temps
    api_output['arrhenius'] = arrhenius_warning

    api_output["score"] = pd.to_numeric(api_output["score"], errors="coerce").fillna(13)
    api_output["adj_kcat"] = pd.to_numeric(api_output["adj_kcat"], errors="coerce")

    # Initialize columns for tie-breaking
    api_output['id_perc'] = -1
    api_output['organism_score'] = np.inf

    # 3. Keep only best-score candidates
    min_score = api_output['score'].min()
    tied = api_output[api_output['score'] == min_score]

    # 4. Tie-breaking
    if len(tied) > 1:
        # Tie-break with enzyme identity
        tied = closest_enz(kcat_dict, tied)
        if not tied['id_perc'].isna().all():
            max_id = tied['id_perc'].max()
            tied = tied[tied['id_perc'] == max_id]

    if len(tied) > 1:
        # Tie-break with taxonomy
        tied = closest_taxonomy(general_criteria, tied)
        if not tied['organism_score'].isna().all():
            min_tax = tied['organism_score'].min()
            tied = tied[tied['organism_score'] == min_tax]

    if len(tied) > 1:
        # Tie-break with max kcat value
        max_kcat = tied['adj_kcat'].max()
        tied = tied[tied['adj_kcat'] == max_kcat]

    # 5. Select best candidate
    best_candidate = tied.iloc[0].to_dict()
    best_candidate['catalytic_enzyme'] = kcat_dict.get('catalytic_enzyme')
    best_score = best_candidate['score']

    # 6. Compute organism_score and id_perc if not present 
    if best_candidate['organism_score'] == np.inf:
        if best_candidate.get('Organism') == general_criteria['Organism']:
            best_candidate['organism_score'] = 0
        else:
            tmp_df = pd.DataFrame([best_candidate])
            taxonomy_score = closest_taxonomy(general_criteria, tmp_df).iloc[0]['organism_score']
            best_candidate['organism_score'] = taxonomy_score

    if best_candidate['id_perc'] == -1:
        catalytic = kcat_dict.get('catalytic_enzyme')
        catalytic_list = str(catalytic).split(";") if catalytic and pd.notna(catalytic) else []

        if best_candidate.get('UniProtKB_AC') in catalytic_list:
            best_candidate['id_perc'] = 100.0
        else:
            tmp_df = pd.DataFrame([best_candidate])
            best_candidate['id_perc'] = closest_enz(kcat_dict, tmp_df).iloc[0]['id_perc']

    return best_score, best_candidate


# --- Utils --- 

def _norm_name(s: str) -> str:
    """Normalize substrates names"""
    if s is None:
        return ""
    s = s.strip().lower()
    # Remove prefixes (d-, l-, d -, l -)
    s = re.sub(r'\b[dl]\s*-\s*', '', s)
    # unify hyphens/spaces
    s = s.replace('-', ' ')
    # compress spaces
    s = re.sub(r'\s+', ' ', s)
    return s

def _to_set(x) -> set:
    """Transform a string 'a; b; c' into a normalized set."""
    if x is None or (isinstance(x, float) and pd.isna(x)):
        return set()
    if isinstance(x, str):
        parts = [p for p in (t.strip() for t in x.split(';')) if p]
    else:
        parts = list(x)
    return { _norm_name(p) for p in parts if p }

def _any_intersection(a, b) -> bool:
    return len(_to_set(a) & _to_set(b)) > 0


# --- Check parameters ---


def check_catalytic_enzyme(candidate, kcat_dict): 
    """
    Checks whether the enzyme in a candidate entry matches the model's enzyme.
    Identifies the catalytic enzyme using UniProt API.
    """
    if pd.notna(kcat_dict['catalytic_enzyme']):
        catalytic_enzymes = str(kcat_dict['catalytic_enzyme']).split(";")
        if candidate["UniProtKB_AC"] in catalytic_enzymes:
            return 0
    return 3


def check_organism(candidate, general_criteria): 
    """
    Checks whether the organism in a candidate entry matches the expected organism.
    """
    if candidate["Organism"] == general_criteria["Organism"]:
        return 0
    return 2


def check_variant(candidate):
    """
    Checks whether the enzyme variant in a candidate entry is wildtype or unknown.
    """
    if candidate["EnzymeVariant"] == "wildtype":
        return 0
    else:  # Unknown
        return 1


def check_pH(candidate, general_criteria):
    """
    Checks whether the pH in a candidate entry matches the expected pH.
    """
    ph_min, ph_max = general_criteria["pH"]
    candidate_ph = candidate.get("pH", None)

    if pd.isna(candidate_ph):
        return 1
    elif ph_min <= candidate_ph <= ph_max:
        return 0
    else:  # Out of range
        return 2
    

def check_temperature(candidate, general_criteria, api_output, min_r2=0.8, expected_range=(50000, 150000)): 
    """
    Checks whether the temperature in a candidate entry matches the expected temperature.
    If the temperature is within the specified range is not met, verify if the Arrhenius equation can be applied.
    """

    temp_min, temp_max = general_criteria["Temperature"]
    candidate_temp = candidate.get("Temperature")

    if pd.isna(candidate_temp):
        return 1, 0
    
    if temp_min <= candidate_temp <= temp_max:
        return 0, 0

    # Try to find a correct the kcat value using the Arrhenius equation
    ph_min, ph_max = general_criteria["pH"]
    
    # Base filters
    filters = (
        api_output["pH"].between(ph_min, ph_max)
        & (api_output["UniProtKB_AC"] == candidate["UniProtKB_AC"])
        & api_output["Temperature"].notna()
        & api_output["value"].notna()
    )
    
    valid_idx = api_output.apply(
        lambda row: check_substrate(row.to_dict(), None, candidate) == 0,
        axis=1
        )
    
    filters = filters & valid_idx

    temps_dispo = api_output.loc[filters, "Temperature"].nunique()
    api_filtered = api_output.loc[filters, ["Temperature", "value"]].copy()

    # Convert temperatures to Kelvin
    api_filtered["Temperature"] = api_filtered["Temperature"] + 273.15

    if temps_dispo >= 2:
        ea, r2, n = calculate_ea(api_filtered)
        if r2 >= min_r2 and ea > 0:
            if not (expected_range[0] <= ea <= expected_range[1]):
                logging.warning(f"{candidate.get('ECNumber')}: Estimated Ea ({ea:.0f} J/mol) is outside the expected range {expected_range} J/mol.")
            # Go Arrhenius
            return 0, n
    return 2, 0


def check_substrate(entry, kcat_dict=None, candidate=None):
    """
    Checks whether the substrate in a candidate entry matches the model's substrates.
    """
    api = entry.get("db", candidate.get("db") if candidate else None)

    # Normalize names
    entry_subs = entry.get("Substrate", "")
    entry_prods = entry.get("Product", "")
    entry_kegg = entry.get("KeggReactionID")

    cand_subs = candidate.get("Substrate", "") if candidate else ""
    cand_prods = candidate.get("Product", "") if candidate else ""
    cand_kegg = candidate.get("KeggReactionID") if candidate else ""

    model_subs = (kcat_dict or {}).get("substrates_name", "")
    model_prods = (kcat_dict or {}).get("products_name", "")
    model_kegg = (kcat_dict or {}).get("rxn_kegg")

    if api == "sabio_rk":
        
        entry_kegg = None if pd.isna(entry_kegg) else entry_kegg
        model_kegg = None if pd.isna(model_kegg) else model_kegg
        cand_kegg  = None if pd.isna(cand_kegg) else cand_kegg
        
        if model_kegg and entry_kegg and _norm_name(model_kegg) == _norm_name(entry_kegg):
            if _any_intersection(entry_subs, model_subs) or _any_intersection(entry_prods, model_prods):
                return 0
        if cand_kegg and entry_kegg and _norm_name(cand_kegg) == _norm_name(entry_kegg):
            if _any_intersection(entry_subs, cand_subs) or _any_intersection(entry_prods, cand_prods):
                return 0
        base_subs = model_subs or cand_subs
        if _any_intersection(entry_subs, base_subs):
            return 0
        return 4

    elif api == "brenda":
        base_subs = model_subs or cand_subs
        if _any_intersection(entry_subs, base_subs):
            return 0
        return 4

    return 4


# --- Scoring ---


def compute_score(kcat_dict, candidate, general_criteria, api_output):
    """
    Compute a score for the candidate based on the kcat dictionary and general criteria.
    """
    score = 0
    # Check catalytic enzyme
    score += check_catalytic_enzyme(candidate, kcat_dict) # + 0 or 3
    # Check organism
    if score != 0: 
        score += check_organism(candidate, general_criteria) # + 0 or 2  
    # Check variant
    score += check_variant(candidate) # + 0, 1
    # Check pH
    score += check_pH(candidate, general_criteria) # + 0, 1 or 2
    # Check substrate 
    score += check_substrate(candidate, kcat_dict) # + 0 or 4
    # Check temperature 
    temperature_penalty, arrhenius = check_temperature(candidate, general_criteria, api_output) # + 0, 1 or 2
    score += temperature_penalty
    return score, arrhenius
