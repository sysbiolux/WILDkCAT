import os 
import datetime
import logging
import pandas as pd
import numpy as np

from ..machine_learning.catapro import create_catapro_input_file, integrate_catapro_predictions
from ..machine_learning.openkineticpredictor import create_prediction_input_file, send_kinetic_predictor, integrate_predictions
from ..utils.generate_reports import report_prediction_input
from ..utils.manage_warnings import DedupFilter


# TODO: Add a warning if there is the same SMILE for multiple KEGG IDs ?
# TODO: Simplify the format of the final output file (remove some columns) 


# --- Format ---


def format_output(kcat_df, limit_penalty_score):
    """
    Formats the kcat DataFrame by selecting and rounding kcat values based on matching score and prediction availability.

    Parameters:
        kcat_df (pandas.DataFrame) : Input DataFrame containing kcat-related columns and prediction results.
        limit_penalty_score (float) : Threshold for the matching score to determine whether to use predicted kcat values.

    Returns: 
        pandas.DataFrame : Formatted DataFrame with selected and rounded kcat values, reordered columns, and updated source information.
    """
    kcat_df = kcat_df.rename(columns={"kcat": "kcat_source", "db": "kcat_source_db"})

    def choose_row(row):
        if pd.notna(row["kcat_source"]):
            if row["penalty_score"] >= limit_penalty_score and pd.notna(row["catapro_predicted_kcat_s"]):
                return pd.Series([row["catapro_predicted_kcat_s"], "catapro"])
            else:
                return pd.Series([row["kcat_source"], row["kcat_source_db"]])
        else:
            if pd.notna(row["catapro_predicted_kcat_s"]):
                return pd.Series([row["catapro_predicted_kcat_s"], "catapro"])
            else:
                return pd.Series([np.nan, np.nan])

    # Add final kcat + db
    kcat_df[["kcat", "kcat_db"]] = kcat_df.apply(choose_row, axis=1)

    # Round numeric columns
    kcat_df["kcat"] = kcat_df["kcat"].round(4)
    kcat_df["kcat_source"] = kcat_df["kcat_source"].round(4)
    kcat_df["kcat_id_percent"] = kcat_df["kcat_id_percent"].round(2)
    kcat_df["catapro_predicted_kcat_s"] = kcat_df["catapro_predicted_kcat_s"].round(1)

    # Drop useless columns
    kcat_df = kcat_df.drop(columns=["kcat_source", "catapro_predicted_kcat_s", "kcat_source_db"])

    # Rename columns 

    kcat_df = kcat_df.rename(columns={
        "kcat_db": "db"
    })

    # If db = catapro then remove the content in the columns "penalty_score", "kcat_substrate", "kcat_organism", "kcat_enzyme", "kcat_temperature", "kcat_ph", "kcat_variant", "kcat_id_percent"
    kcat_df.loc[kcat_df['db'] == 'catapro', [
        "warning_arr",
        "penalty_score", "kcat_substrate", "kcat_organism", "kcat_enzyme", 
        "kcat_temperature", "kcat_ph", "kcat_variant", "kcat_id_percent", "kcat_organism_score"
        ]] = np.nan
    
    # Reorder columns
    kcat_df = kcat_df[[
        "rxn", "rxn_kegg", "ec_code", "ec_codes", "direction", 
        "substrates_name", "substrates_kegg", "products_name", "products_kegg", 
        "genes", "uniprot", "catalytic_enzyme", 
        "warning_ec", "warning_enz", "warning_arr", 
        "kcat", "db", 
        "penalty_score", "kcat_substrate", "kcat_organism", "kcat_enzyme", "kcat_temperature", "kcat_ph", "kcat_variant", "kcat_id_percent", "kcat_organism_score"
        ]]

    return kcat_df


def format_output_kxp(kcat_df, limit_penalty_score, method):
    """
    Formats the kcat DataFrame by selecting and rounding kcat values based on matching score and prediction availability.

    Parameters:
        kcat_df (pandas.DataFrame) : Input DataFrame containing kcat-related columns and prediction results.
        limit_penalty_score (float) : Threshold for the matching score to determine whether to use predicted kcat values.
        method (str): TODO 
    Returns: 
        pandas.DataFrame : Formatted DataFrame with selected and rounded kcat values, reordered columns, and updated source information.
    """
    kcat_df = kcat_df.rename(columns={"kcat": "kcat_source", "db": "kcat_source_db"})

    def choose_row(row):
        if pd.notna(row["kcat_source"]):
            if row["penalty_score"] >= limit_penalty_score and pd.notna(row["predicted_kcat_s"]):
                return pd.Series([row["predicted_kcat_s"], method])
            else:
                return pd.Series([row["kcat_source"], row["kcat_source_db"]])
        else:
            if pd.notna(row["predicted_kcat_s"]):
                return pd.Series([row["predicted_kcat_s"], method])
            else:
                return pd.Series([np.nan, np.nan])

    # Add final kcat + db
    kcat_df[["kcat", "kcat_db"]] = kcat_df.apply(choose_row, axis=1)

    # Round numeric columns
    kcat_df["kcat"] = kcat_df["kcat"].round(4)
    kcat_df["kcat_source"] = kcat_df["kcat_source"].round(4)
    kcat_df["kcat_id_percent"] = kcat_df["kcat_id_percent"].round(2)
    kcat_df["predicted_kcat_s"] = kcat_df["predicted_kcat_s"].round(1)

    # Drop useless columns
    kcat_df = kcat_df.drop(columns=["kcat_source", "predicted_kcat_s", "kcat_source_db"])

    # Rename columns 
    kcat_df = kcat_df.rename(columns={
        "kcat_db": "db"
    })

    # If db = method then remove the content in the columns "penalty_score", "kcat_substrate", "kcat_organism", "kcat_enzyme", "kcat_temperature", "kcat_ph", "kcat_variant", "kcat_id_percent"
    kcat_df.loc[kcat_df['db'] == method, [
        "penalty_score", "kcat_substrate", "kcat_organism", "kcat_enzyme", 
        "kcat_temperature", "kcat_ph", "kcat_variant", "kcat_id_percent", "kcat_organism_score"
        ]] = np.nan
    
    # Reorder columns
    kcat_df = kcat_df[[
        "rxn", "rxn_kegg", "ec_code", "ec_codes", "direction", 
        "substrates_name", "substrates_kegg", "products_name", "products_kegg", 
        "genes", "uniprot", "catalytic_enzyme", 
        "warning_ec", "warning_enz", "warning_arr", 
        "kcat", "db", 
        "penalty_score", "kcat_substrate", "kcat_organism", "kcat_enzyme", "kcat_temperature", "kcat_ph", "kcat_variant", "kcat_id_percent", "kcat_organism_score"
        ]]

    return kcat_df



# --- Main ---


def run_prediction_part1(output_folder: str,
                         limit_penalty_score: int, 
                         report: bool = True) -> None:
    """
    Processes kcat data file to generate input files for CataPro prediction.
    Optionally, it can produce a summary report of the processed data.

    Parameters:
        output_folder (str): Path to the output folder where the results will be saved.
        limit_penalty_score (int): Threshold for filtering entries based on matching score.
        report (bool, optional): Whether to generate a report using the retrieved data (default: True). 
    """
    # Intitialize logging
    os.makedirs(os.path.join(output_folder, "logs"), exist_ok=True)
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    filename = f"logs/prediction1_{timestamp}.log"
    logging.getLogger().addFilter(DedupFilter())
    logging.basicConfig(filename=os.path.join(output_folder, filename), encoding='utf-8', level=logging.INFO)

    # Run prediction part 1
    # Read the kcat file
    if not os.path.exists(output_folder):
        raise FileNotFoundError(f"The specified output folder '{output_folder}' does not exist.")
    
    kcat_file_path = os.path.join(output_folder, "kcat_retrieved.tsv")
    if not os.path.isfile(kcat_file_path):
        raise FileNotFoundError(f"The specified file '{kcat_file_path}' does not exist in the output folder. Please run the function 'run_retrieval()' first.")

    kcat_df = pd.read_csv(kcat_file_path, sep='\t')

    # Subset rows with no values or matching score above the limit
    kcat_df = kcat_df[(kcat_df['penalty_score'] >= limit_penalty_score) | (kcat_df['penalty_score'].isnull())]
    # Drop rows with no UniProt ID or no substrates_kegg
    before_duplicates_filter = len(kcat_df) - 1 
    kcat_df = kcat_df[kcat_df['uniprot'].notnull() & kcat_df['substrates_kegg'].notnull()]
    nb_missing_enzymes = before_duplicates_filter - len(kcat_df) + 1 
    
    # Generate CataPro input file
    catapro_input_df, substrates_to_smiles_df, report_statistics = create_catapro_input_file(kcat_df)

    # Save the CataPro input file and substrates to SMILES mapping
    os.makedirs(os.path.join(output_folder, "machine_learning"), exist_ok=True)
    output_path = os.path.join(output_folder, "machine_learning/catapro_input.csv")
    kcat_df.to_csv(output_path, sep='\t', index=False)
    catapro_input_df.to_csv(output_path, sep=',', index=True)
    substrates_to_smiles_df.to_csv(output_path.replace('.csv', '_substrates_to_smiles.tsv'), sep='\t', index=False)
    logging.info(f"Output saved to '{output_path}'")

    # Add statistics 
    report_statistics["missing_enzymes"] = nb_missing_enzymes

    nb_predicted_transport = 0
    for _, row in kcat_df.iterrows():         
        if row['substrates_kegg'] == row['products_kegg'] and row['penalty_score'] >= limit_penalty_score: 
            nb_predicted_transport += 1
    
    report_statistics['predicted_transport'] = nb_predicted_transport

    if report:
        report_prediction_input(catapro_input_df, report_statistics, output_folder)


def run_prediction_part2(output_folder: str,
                         catapro_predictions_path: str,
                         limit_penalty_score: int) -> None:
    """
    Runs the second part of the kcat prediction pipeline by integrating Catapro predictions,
    mapping substrates to SMILES and formatting the output. 
    
    Parameters:
        output_folder (str): Path to the output folder where the results will be saved.
        catapro_predictions_path (str): Path to the CataPro predictions CSV file.
        limit_penalty_score (float): Threshold for taking predictions over retrieved values.
    """ 
    # Intitialize logging
    os.makedirs(os.path.join(output_folder, "logs"), exist_ok=True)
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    filename = f"logs/prediction2_{timestamp}.log"
    logging.getLogger().addFilter(DedupFilter())
    logging.basicConfig(filename=os.path.join(output_folder, filename), encoding='utf-8', level=logging.INFO)

    # Run prediction part 2
    if not os.path.exists(output_folder):
        raise FileNotFoundError(f"The specified output folder '{output_folder}' does not exist.")
    kcat_file_path = os.path.join(output_folder, "kcat_retrieved.tsv")
    if not os.path.isfile(kcat_file_path):
        raise FileNotFoundError(f"The specified file '{kcat_file_path}' does not exist in the output folder. Please run the function 'run_extraction()' first.")
    kcat_df = pd.read_csv(kcat_file_path, sep='\t')
    substrates_to_smiles_path = os.path.join(output_folder, "machine_learning/catapro_input_substrates_to_smiles.tsv")
    substrates_to_smiles = pd.read_csv(substrates_to_smiles_path, sep='\t')
    catapro_predictions_df = pd.read_csv(catapro_predictions_path, sep=',')
    kcat_df = integrate_catapro_predictions(kcat_df, 
                                            substrates_to_smiles,
                                            catapro_predictions_df
                                            )

    # Save the output as a TSV file
    kcat_df = format_output(kcat_df, limit_penalty_score)

    # Add warning for the predicted transport reactions in the log file
    for _, row in kcat_df.iterrows():         
        if row['substrates_kegg'] == row['products_kegg'] and row['db'] == 'catapro': 
            logging.warning(f"kcat for transport reaction: '{row['rxn']}' was ML-predicted. Transport reaction predictions are less reliable due to sparse training data")
    
    output_path = os.path.join(output_folder, "kcat_full.tsv")
    kcat_df.to_csv(output_path, sep='\t', index=False)
    logging.info(f"Output saved to '{output_path}'")


def run_prediction(output_folder: str,
                   limit_penalty_score: int, 
                   method: str = "CataPro",
                   canonicalize_substrates: bool = True,
                   report: bool = True) -> None:
    """
    Processes kcat data file to generate input files for CataPro prediction.
    Optionally, it can produce a summary report of the processed data.

    Parameters:
        output_folder (str): Path to the output folder where the results will be saved.
        limit_penalty_score (int): Threshold for filtering entries based on matching score.
        method (str, optional): The method to use for prediction, check the full list of available methods from https://predictor.openkinetics.org/api/v1/methods/ (default: "CataPro").
        canonicalize_substrates (bool, optional): Whether to canonicalize substrates (default: True).
        report (bool, optional): Whether to generate a report using the retrieved data (default: True).
    """
    
    # TODO: Integrate the possibility to get back to a specific step if the code crashed 

    # Intitialize logging
    os.makedirs(os.path.join(output_folder, "logs"), exist_ok=True)
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    filename = f"logs/prediction1_{timestamp}.log"
    logging.getLogger().addFilter(DedupFilter())
    logging.basicConfig(filename=os.path.join(output_folder, filename), encoding='utf-8', level=logging.INFO)

    # 1. Create the prediction file 
    # Read the kcat file
    if not os.path.exists(output_folder):
        raise FileNotFoundError(f"The specified output folder '{output_folder}' does not exist.")
    
    kcat_file_path = os.path.join(output_folder, "kcat_retrieved.tsv")
    if not os.path.isfile(kcat_file_path):
        raise FileNotFoundError(f"The specified file '{kcat_file_path}' does not exist in the output folder. Please run the function 'run_retrieval()' first.")
    
    all_methods = {"CataPro", "CatPred", "DLKcat", "EITLEM", "KinForm-H", "KinForm-L", "MMISA-KM", "OmniESI", "UniKP",
                   "CatPred", "TurNup"}
    
    if method not in all_methods:
        raise ValueError(f"Unknown method '{method} - Valid options: {(all_methods)}")

    output_path = os.path.join(output_folder, "machine_learning/prediction_input.csv")
    if not os.path.isfile(output_path): # TODO: Add a warning message stating that the step 1 is already process
        kcat_df = pd.read_csv(kcat_file_path, sep='\t')

        # Subset rows with no values or matching score above the limit
        kcat_df = kcat_df[(kcat_df['penalty_score'] >= limit_penalty_score) | (kcat_df['penalty_score'].isnull())]
        # Drop rows with no UniProt ID or no substrates_kegg
        before_duplicates_filter = len(kcat_df) - 1 
        kcat_df = kcat_df[kcat_df['uniprot'].notnull() & kcat_df['substrates_kegg'].notnull()]
        nb_missing_enzymes = before_duplicates_filter - len(kcat_df) + 1 
        
        # Generate prediction input file
        prediction_input_df, substrates_to_smiles_df, uniprot_to_sequence_df, report_statistics = create_prediction_input_file(kcat_df, method)

        # Save the input file and substrates to SMILES mapping
        os.makedirs(os.path.join(output_folder, "machine_learning"), exist_ok=True)
        output_path = os.path.join(output_folder, "machine_learning/prediction_input.csv")
        kcat_df.to_csv(output_path, sep='\t', index=False)
        prediction_input_df.to_csv(output_path, sep=',', index=False)
        substrates_to_smiles_df.to_csv(output_path.replace('.csv', '_substrates_to_smiles.tsv'), sep='\t', index=False)
        uniprot_to_sequence_df.to_csv(output_path.replace('.csv', '_uniprot_to_sequence.tsv'), sep='\t', index=False)
        logging.info(f"Output saved to '{output_path}'")

        # Add statistics 
        report_statistics["missing_enzymes"] = nb_missing_enzymes
        nb_predicted_transport = 0
        for _, row in kcat_df.iterrows():         
            if row['substrates_kegg'] == row['products_kegg'] and row['penalty_score'] >= limit_penalty_score: 
                nb_predicted_transport += 1
        report_statistics['predicted_transport'] = nb_predicted_transport

        if report:
            report_prediction_input(prediction_input_df, report_statistics, output_folder)

    # 2. Prediction 
    prediction_output = send_kinetic_predictor(output_folder, method, canonicalize_substrates)

    # 3. Integrate prediction to kcat 
    kcat_df = pd.read_csv(kcat_file_path, sep='\t')
    kcat_df = integrate_predictions(kcat_df, prediction_output, method)
    
    # Save the output as a TSV file
    kcat_df = format_output_kxp(kcat_df, limit_penalty_score, method)

    # Add warning for the predicted transport reactions in the log file
    for _, row in kcat_df.iterrows():         
        if row['substrates_kegg'] == row['products_kegg'] and row['db'] == method: # 
            logging.warning(f"kcat for transport reaction: '{row['rxn']}' was ML-predicted. Transport reaction predictions are less reliable due to sparse training data")
    
    output_path = os.path.join(output_folder, "kcat_full.tsv")
    kcat_df.to_csv(output_path, sep='\t', index=False)
    logging.info(f"Output saved to '{output_path}'")


if __name__ == "__main__":

    # TEST 
    # run_prediction("in_progress/models/e_coli_core", 9, "CataPro")
    run_prediction("in_progress/models/e_coli_core", 9, "KinForm-L")
 
    # Test : Retrieve Sequence from UniProt ID
    # print(convert_uniprot_to_sequence("P0A796"))
    
    # Test : Main function

    # run_prediction_part1("output/ecoli_kcat_brenda.tsv", -1, "output/machine_learning/ecoli_catapro_input.csv")
    # run_prediction_part2("output/ecoli_kcat_brenda.tsv", 
    #                      "output/machine_learning/ecoli_catapro_output.csv", 
    #                      "output/machine_learning/ecoli_catapro_input_substrates_to_smiles.tsv", 
    #                      8, 
    #                      "output/ecoli_kcat_full.tsv")
    
    # run_prediction_part1("output/yeast_kcat_brenda.tsv", -1, "output/machine_learning/yeast_catapro_input.csv")
    # run_prediction_part2("output/yeast_kcat_brenda.tsv", 
    #                      "output/machine_learning/yeast_catapro_output.csv", 
    #                      "output/machine_learning/yeast_catapro_input_substrates_to_smiles.tsv", 
    #                      8, 
    #                      "output/yeast_kcat_full.tsv")