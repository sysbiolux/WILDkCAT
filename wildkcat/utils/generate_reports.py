import os
import io
import base64
import logging
import datetime
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import plotly.express as px
from matplotlib.ticker import LogFormatter, MaxNLocator
from io import BytesIO


def report_extraction(model, df, report_statistics, output_folder, shader=False) -> None:
    """
    Generates a detailed HTML report summarizing kcat extraction results from a metabolic model.

    Parameters: 
        model (cobra.Model): The metabolic model object containing reactions, metabolites, and genes.
        df (pandas.DataFrame): DataFrame containing data from the run_extraction function.
        report_statistics (dict): Dictionary with statistics about EC code assignment and extraction issues.
        output_folder (str): Path to the output folder where the report will be saved.
        shader (bool, optional): If True, includes a shader canvas background in the report. Default is False.

    Returns: 
        None: The function saves the generated HTML report to 'reports/extract_kcat_report.html'. 
    """
    # Model statistics
    nb_model_reactions = len(model.reactions)
    nb_model_metabolites = len(model.metabolites)
    nb_model_genes = len(model.genes)
    rxn_with_ec = 0
    unique_ec_codes = []

    for rxn in model.reactions:
        ec_code = rxn.annotation.get('ec-code')
        if ec_code:
            rxn_with_ec += 1
            if isinstance(ec_code, str):
                ec_code = [ec_code.strip()]
            elif isinstance(ec_code, list):
                ec_code = [x.strip() for x in ec_code if x.strip()]
            else:
                ec_code = []
            unique_ec_codes.extend(ec_code)
    
    nb_model_ec_codes = len(set(unique_ec_codes))

    # Kcat statistics
    nb_reactions = df['rxn'].nunique()
    nb_ec_codes = df.loc[df["warning_ec"].fillna("") == "", "ec_code"].nunique()


    # nb_missing_ec = report_statistics.get('nb_missing_ec', np.nan)
    nb_incomplete_ec = report_statistics.get('nb_incomplete_ec', np.nan)
    nb_transferred_ec = report_statistics.get('nb_transferred_ec', np.nan)
    nb_missing_gpr = report_statistics.get('nb_missing_gpr', np.nan)
    nb_missing_catalytic_enzyme = report_statistics.get('nb_missing_catalytic_enzyme', 0)
    # nb_multiple_catalytic_enzymes = report_statistics.get('nb_multiple_catalytic_enzymes', np.nan)
    nb_of_lines_dropped_no_ec_no_enzyme = report_statistics.get('nb_of_lines_dropped_no_ec_no_enzyme', np.nan)
    nb_of_reactions_dropped_no_ec_no_enzyme = report_statistics.get('nb_of_reactions_dropped_no_ec_no_enzyme', np.nan)

    rxn_coverage = 100.0 * nb_reactions / nb_model_reactions if nb_model_reactions else 0
    # percent_ec_retrieved = 100.0 * nb_ec_codes / nb_model_ec_codes if nb_model_ec_codes else 0

    rxn_ec_coverage = 100.0 * rxn_with_ec / nb_model_reactions if nb_model_ec_codes else 0

    # Pie Chart
    pie_data = {
        "Retrieved": nb_ec_codes,
        "Transferred": nb_transferred_ec, 
        "Incomplete": nb_incomplete_ec,
    }

    pie_data = {k: v for k, v in pie_data.items() if v > 0}

    fig = px.pie(
        names=list(pie_data.keys()),
        values=list(pie_data.values()),
        color_discrete_sequence=["#55bb55", "#ee9944", "#cc4455"]
    )
    fig.update_traces(textinfo="percent+label", textfont_size=16)
    fig.update_layout(
        title="",
        title_font=dict(size=30, color="black"),
        showlegend=True
    )

    pie_chart_html = fig.to_html(full_html=False, include_plotlyjs="cdn")

    # Time
    generated_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    # Html report
    html = f"""
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>Extract kcat Report</title>
        {report_style()}
    </head>
    <body>
        <header>
            <canvas id="shader-canvas"></canvas>
            <div class="overlay">
                <h1>Extract k<sub>cat</sub> Report</h1>
                <p>Generated on {generated_time}</p>
            </div>
        </header>

        <div class="container">
            <!-- Model Overview -->
            <div class="card">
                <h2>Model Overview</h2>
                <div class="stats-grid">
                    <div class="stat-box">
                        <h3>{model.id}</h3>
                        <p>Model ID</p>
                    </div>
                    <div class="stat-box">
                        <h3>{nb_model_reactions}</h3>
                        <p>Reactions</p>
                    </div>
                    <div class="stat-box">
                        <h3>{nb_model_metabolites}</h3>
                        <p>Metabolites</p>
                    </div>
                    <div class="stat-box">
                        <h3>{nb_model_genes}</h3>
                        <p>Genes</p>
                    </div>
                </div>
            </div>

            <!-- kcat Extraction Table -->
            <div class="card">
                <h2>k<sub>cat</sub> Extraction Statistics</h2>
                <table>
                    <tr>
                        <th>Metric</th>
                        <th>Value</th>
                        <th>Visualization</th>
                    </tr>
                    <tr>
                        <td>Reaction with k<sub>cat</sub> information</td>
                        <td>{nb_reactions} ({rxn_coverage:.1f}%)</td>
                        <td>
                            <div class="progress">
                                <div class="progress-bar-table" style="width:{rxn_coverage}%;"></div>
                            </div>
                        </td>
                    </tr>
                    <tr>
                        <td>Reactions with EC information</td>
                        <td>{rxn_with_ec} ({rxn_ec_coverage:.1f}%)</td>
                        <td>
                            <div class="progress">
                                <div class="progress-bar-table" style="width:{rxn_ec_coverage}%;"></div>
                            </div>
                        </td>
                    </tr>
                    <tr>
                        <td>Total k<sub>cat</sub> in output</td>
                        <td>{len(df)}</td>
                        <td>-</td>
                    </tr>
                </table>
            </div>

            <!-- EC Issues Table -->
            <div class="card">
                <h2>Quality Control</h2>

                <p style="text-align: justify">
                    During the extraction process, several issues may arise that can impact the 
                    quality of the retrieved k<sub>cat</sub> data.
                    <br>
                    When an EC number is incomplete or has been transferred, WILDkCAT performs the 
                    retrieval based only on the available enzyme information. Such cases are indicated as 
                    'incomplete' or 'transferred' in the 'warning_ec' column. 
                    These situations reduce the likelihood of finding alternative k<sub>cat</sub> values.
                    <br>
                    Additionally, WILDkCAT attempts to identify the catalytic enzyme associated with 
                    each reaction. If no Gene–Protein–Reaction rule is available, or if the catalytic 
                    enzyme cannot be found via the UniProt API, the entry is labeled as 
                    'none' in the 'warning_enz' column.
                    <br>
                    If both the EC number and the catalytic enzyme information are missing for a given 
                    reaction, the corresponding row is removed due to insufficient information to assign 
                    a k<sub>cat</sub> value.
                </p>

                <table>
                    <tr>
                        <th>Cases</th>
                        <th>Count</th>
                    </tr>
                    <tr>
                        <td>Transferred EC codes</td>
                        <td>{nb_transferred_ec}</td>
                    </tr>
                    <tr>
                        <td>Incomplete EC codes</td>
                        <td>{nb_incomplete_ec}</td>
                    </tr>
                    <tr>
                        <td>Number of reactions without catalytic enzyme</td>
                        <td>{nb_missing_gpr + nb_missing_catalytic_enzyme}</td>
                    </tr>
                    <tr>
                        <td>Number of reactions dropped due to inconsistent or absent EC codes and enzymes</td>
                        <td>{nb_of_reactions_dropped_no_ec_no_enzyme}</td>
                    </tr>
                    <tr>
                        <td>Number of k<sub>cat</sub> values dropped due to inconsistent or absent EC codes and enzymes</td>
                        <td>{nb_of_lines_dropped_no_ec_no_enzyme}</td>
                    </tr>
                </table>
            </div>

            <!-- Pie Chart Section -->
            <div class="card">
                <h2>EC Distribution</h2>
                {pie_chart_html}
            </div>
        </div>

        <footer>WILDkCAT</footer>
    """
    if shader:
        html += report_shader()
    else: 
        html += report_simple()
    html += """
    </body>
    </html>
    """

    # Save report
    os.makedirs(os.path.join(output_folder, "reports"), exist_ok=True)
    report_path = os.path.join(output_folder, "reports/extract_report.html")
    with open(report_path, "w", encoding="utf-8") as f:
        f.write(html)
    logging.info(f"HTML report saved to '{report_path}'")


def report_retrieval(df, output_folder, parameters, shader=False) -> None:
    """
    Generate a styled HTML report summarizing the kcat matching results,
    including kcat value distribution and matching score repartition.

    Parameters:
        df (pd.DataFrame): DataFrame containing data from the run_retrieval function.
        output_folder (str): Path to the output folder where the report will be saved.
        shader (bool, optional): If True, includes a shader canvas background in the report. Default is False.

    Returns:
        None: The function saves the generated HTML report to 'reports/retrieve_kcat_report.html'.
    """
    # Ensure numeric kcat values to avoid TypeError on comparisons
    kcat_values = pd.to_numeric(df['kcat'], errors='coerce').dropna()

    # Only use scores present in the data
    present_scores = sorted(df['penalty_score'].dropna().unique())
    score_counts = df['penalty_score'].value_counts().reindex(present_scores, fill_value=0)
    total = len(df)
    matched = len(kcat_values)
    match_percent = matched / total * 100 if total else 0
    score_percent = (score_counts / total * 100).round(2) if total else pd.Series(0, index=present_scores)

    # Gradient colors from green (best score) to red (worst score) # TODO: It could be better to create the scale dynamically based on present scores
    distinct_colors = [
        "#27ae60",
        "#43b76e",
        "#60c07c",
        "#7cc98a",
        "#98d298",
        "#b5dbb6",
        "#d1e4c4",
        "#e8e9b9",
        "#f1e9b6",
        "#f7d97c",
        "#f9c74f",
        "#f8961e",
        "#f3722c",
        "#e67e22",
        "#e74c3c",
        "#d35400",
        "#c0392b",
        "#a93226",
        "#7b241c"
    ]

    def score_color(score):
        idx = present_scores.index(score)
        return distinct_colors[idx % len(distinct_colors)]

    generated_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    
    # Parameters
    temp = parameters.get('Temperature')
    formatted_temp = f"{temp[0]} - {temp[1]}"

    ph = parameters.get('pH')
    formatted_ph = f"{ph[0]} - {ph[1]}"

    # Histogram with stacked bars for scores
    kcat_hist_base64 = ""
    if not kcat_values.empty:
        min_exp = int(np.floor(np.log10(max(1e-6, kcat_values.min()))))
        max_exp = int(np.ceil(np.log10(kcat_values.max())))
        bins = np.logspace(min_exp, max_exp, num=40)

        # Rm empty score groups (15 - 16)
        hist_data = []
        valid_scores = []
        for score in present_scores:
            vals = pd.to_numeric(df[df['penalty_score'] == score]['kcat'], errors='coerce')
            vals = vals[vals.notna()]
            if not vals.empty:
                hist_data.append(vals)
                valid_scores.append(score)

        fig, ax = plt.subplots(figsize=(12, 6))
        
        # Stacked histogram by score
        ax.hist(hist_data, bins=bins, stacked=True, 
                color=[score_color(s) for s in valid_scores],
                label=[f"{s}" for s in valid_scores],
                edgecolor='white')
        
        ax.set_xscale('log')
        ax.set_xlim([10**min_exp / 1.5, 10**max_exp * 1.5])
        ax.xaxis.set_major_formatter(LogFormatter(10))
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))

        ax.set_xlabel("kcat (s⁻¹)", fontsize=12)
        ax.set_ylabel("Count", fontsize=12)
        ax.set_title(f"", fontsize=13)
        
        ax.legend(
            title="Penalty Score", 
            fontsize=10, 
            title_fontsize=11,
            loc='center left', 
            bbox_to_anchor=(1, 0.5),
            frameon=False
        )
        
        # Style
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_color('#444444')
        ax.spines['bottom'].set_color('#444444')

        ax.grid(True, which='major', axis='y', linestyle='--', linewidth=0.6, alpha=0.4)
        ax.grid(False, which='major', axis='x') 

        plt.tight_layout(rect=[0, 0, 0.85, 1])
        buf = io.BytesIO()
        plt.savefig(buf, format='png', bbox_inches='tight')
        plt.close(fig)
        kcat_hist_base64 = base64.b64encode(buf.getvalue()).decode('utf-8')

    # HTML start
    html = f"""
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>Retrieve kcat Report</title>
        {report_style()}
    </head>
    <body>
        <header>
            <canvas id="shader-canvas"></canvas>
            <div class="overlay">
                <h1>Retrieve k<sub>cat</sub> Report</h1>
                <p>Generated on {generated_time}</p>
            </div>
        </header>

        <div class="container">
            <div class="card">
                <h2>Overview</h2>
                <div class="stats-grid">
                    <div class="stat-box">
                        <h3>{total}</h3>
                        <p>Total Entries</p>
                    </div>
                    <div class="stat-box">
                        <h3>{matched}</h3>
                        <p>Matched k<sub>cat</sub> ({match_percent:.2f}%)</p>
                    </div>
                </div>
            </div>

            <div class="card">
                <div class="card">
                <h2>Parameters</h2>
                <div class="stats-grid">
                    <div class="stat-box">
                        <h3>{parameters.get('Organism')}</h3>
                        <p>Organism name</p>
                    </div>
                    <div class="stat-box">
                        <h3>{formatted_ph}</h3>
                        <p>pH range</p>
                    </div>
                    <div class="stat-box">
                        <h3>{formatted_temp}</h3>
                        <p>Temperature range</p>
                    </div>
                     <div class="stat-box">
                        <h3>{parameters.get('database').capitalize()}</h3>
                        <p>Database(s)</p>
                    </div>
                </div>
            </div>
            </div>

            <div class="card">
                <h2>Penalty Score Distribution</h2>
                <div class="progress-stacked">
    """

    # Add progress bars only for present scores
    for score in present_scores:
        percent = score_percent.get(score, 0)
        if percent > 0:
            html += f'<div class="progress-bar" style="width:{percent}%;background:{score_color(score)};" title="Score {score}: {percent:.2f}%"></div>'

    html += """
            </div>
            <div class="legend">
    """

    # Add legend only for present scores
    for score in present_scores:
        html += f'<div class="legend-item"><div class="legend-color" style="background:{score_color(score)};"></div> Score {score}</div>'

    html += """
            </div>
            <table>
                <tr>
                    <th>Score</th>
                    <th>Count</th>
                    <th>Percent</th>
                </tr>
    """

    # Table rows only for present scores
    for score in present_scores:
        html += f'<tr><td>{score}</td><td>{score_counts[score]}</td><td>{score_percent[score]:.2f}%</td></tr>'

    html += """
            </table>
        </div>
    """

    # Histogram section (stacked by score)
    html += """
        <div class="card">
            <h2>Distribution of k<sub>cat</sub> values (Stacked by Penalty Score)</h2>
            <div class="img-section">
    """
    if kcat_hist_base64:
        html += f'<img src="data:image/png;base64,{kcat_hist_base64}" alt="k<sub>cat</sub> Distribution">'
    html += """
            </div>
        </div>
    """

    # Metadata section
    html += f"""
            <div class="card">
                <h2>Penalty Score</h2>
                <p>
                    The penalty score evaluates how well a candidate k<sub>cat</sub> entry fits the query enzyme and conditions. 
                    A lower score indicates a better match (0 = Best possible, 16 = No match).
                </p>
                <h3>Scoring process:</h3>
                <ul>
                    <li><b>Catalytic enzyme:</b> Check if the reported enzyme matches the expected catalytic enzyme(s).</li>
                    <li><b>Organism:</b> Penalize mismatches between the source organism and the target organism.</li>
                    <li><b>Enzyme variant:</b> Exclude or penalize mutant/engineered variants (wildtype preferred).</li>
                    <li><b>pH:</b> Check whether the reported pH is consistent with the desired experimental range.</li>
                    <li><b>Substrate:</b> Verify substrate compatibility with the catalytic reaction.</li>
                    <li><b>Temperature:</b> Penalize deviations from the target temperature; 
                        if possible, adjust kcat values using the Arrhenius equation.</li>
                </ul>

                <h3>Score breakdown:</h3>
                <table border="1" cellpadding="6" cellspacing="0" style="border-collapse: collapse; text-align: left;">
                    <tr>
                        <th>Criterion</th>
                        <th>Penalty</th>
                    </tr>
                    <tr>
                        <td>Substrate mismatch</td>
                        <td>+4</td>
                    </tr>
                    <tr>
                        <td>Catalytic enzyme mismatch</td>
                        <td>+3</td>
                    </tr>
                    <tr>
                        <td>Organism mismatch</td>
                        <td>+2</td>
                    </tr>
                    <tr>
                        <td>pH unknown</td>
                        <td>+1</td>
                    </tr>
                    <tr>
                        <td>pH out of range</td>
                        <td>+2</td>
                    </tr>
                    <tr>
                        <td>Temperature unknown</td>
                        <td>+1</td>
                    </tr>
                    <tr>
                        <td>Temperature out of range</td>
                        <td>+2</td>
                    </tr>
                    <tr>
                        <td>Enzyme variant unknown</td>
                        <td>+1</td>
                    </tr>
                </table>

                <p>
                    Candidates are then ranked by:
                    <ol>
                        <li>Lowest penalty-score</li>
                        <li>Highest sequence identity percentage to the target enzyme</li>
                        <li>Closest organism compared to the target organism.</li>
                        <li>Adjusted k<sub>cat</sub> value (favoring the highest value by default)</li>
                    </ol>

                <i>Please check the <a href="https://h-escoffier.github.io/WILDkCAT/explanation/explanation/#2-retrieve-experimental-kcat-values-from-brenda-andor-sabio-rk" target="_blank" rel="noopener noreferrer">documentation</a> for more details on the scoring system and the retrieval process.</i>
                </p>
            </div>
            
            <div class="card">
                <h2>Notes</h2>
                <p style="text-align: justify">
                    Please note that the number of rows may differ between the extraction and retrieval stages. 
                    Indeed, when a single reaction–enzyme combination is associated with multiple EC numbers, WILDkCAT 
                    automatically merges these rows after retrieval, keeping only the best entry according to the criteria 
                    described above. The EC number of the selected kcat is stored in the 'ec_code' column, while all EC 
                    numbers associated with the reaction are stored in the 'ec_codes' column.
                </p>
            </div>
        </div>

        <footer>WILDkCAT</footer>
    """
    if shader: 
        html += report_shader()
    else: 
        html += report_simple()
    html += """
    </body>
    </html>
    """

    # Save HTML
    os.makedirs(os.path.join(output_folder, "reports"), exist_ok=True)
    report_path = os.path.join(output_folder, "reports/retrieve_report.html")
    with open(report_path, "w", encoding="utf-8") as f:
        f.write(html)

    logging.info(f"HTML report saved to '{report_path}'")


def report_prediction_input(catapro_df, report_statistics, output_folder, shader=False) -> None: 
    """
    Generate a detailed HTML report summarizing the kcat prediction input statistics.

    Parameters:
        catapro_df (pd.DataFrame): DataFrame containing the CataPro input data.
        report_statistics (dict): Dictionary with statistics about the prediction input.
        output_folder (str): Path to the output folder where the report will be saved.
        shader (bool, optional): If True, includes a shader canvas background in the report. Default is False.

    Returns:
        None: The function saves the generated HTML report to 'reports/predict_kcat_report.html'.
    """
    # CataPro Statistics 
    total_catapro_entries = len(catapro_df) - 1

    # Report Statistics
    rxn_covered = report_statistics['rxn_covered']
    cofactors_covered = report_statistics['cofactor_identified']
    no_catalytic = report_statistics['no_catalytic']
    kegg_missing = report_statistics['kegg_no_matching']
    duplicates = report_statistics['duplicates_enzyme_substrates']
    missing_enzyme = report_statistics['missing_enzymes']
    predicted_transport =  report_statistics['predicted_transport']

    total_rxn = rxn_covered + no_catalytic + kegg_missing + missing_enzyme
    rxn_coverage = (rxn_covered / total_rxn * 100) if total_rxn > 0 else 0

    # Time
    generated_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    # Html report
    html = f"""
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>Predict kcat Report</title>
        {report_style()}
    </head>
    <body>
        <header>
            <canvas id="shader-canvas"></canvas>
            <div class="overlay">
                <h1>Predict k<sub>cat</sub> Report</h1>
                <p>Generated on {generated_time}</p>
            </div>
        </header>

        <div class="container">
            <!-- CataPro Overview -->
            <div class="card">
                <h2>Overview</h2>
                <div class="stats-grid">
                    <div class="stat-box">
                        <h3>{total_rxn}</h3>
                        <p>Total k<sub>cat</sub> values</p>
                    </div>
                    <div class="stat-box">
                        <h3>{rxn_covered}</h3>
                        <p>k<sub>cat</sub> to be predicted ({rxn_coverage:.2f}%)</p>
                    </div>
                </div>
            </div>

            <!-- Prediction kcat Table -->
            <div class="card">
                <h2>k<sub>cat</sub> Prediction Statistics</h2>
                <table>
                    <tr>
                        <th>Metric</th>
                        <th>Value</th>
                    </tr>
                    <tr>
                        <td>Total of entries in CataPro input file</td>
                        <td>{total_catapro_entries}</td>
                    </tr>
                    <tr>
                        <td>Number of cofactor identified</td>
                        <td>{cofactors_covered}</td>
                    </tr>
                </table>
            </div>

            <div class="card">
                <h2>Issues in k<sub>cat</sub> Predictions</h2>
                <table>
                    <tr>
                        <th>Metric</th>
                        <th>Value</th>
                    </tr>
                    <tr>
                        <td>Entries with no catalytic enzyme identified</td>
                        <td>{no_catalytic}</td>
                    </tr>
                    <tr>
                        <td>Entries with missing KEGG IDs</td>
                        <td>{kegg_missing}</td>
                    </tr>
                    <tr>
                        <td>Entries with missing enzyme information</td>
                        <td>{missing_enzyme}</td>
                    </tr>
                </table>
            </div>

            <div class="card">
                <h2>Duplicates</h2>
                <table>
                    <tr>
                        <th>Metric</th>
                        <th>Value</th>
                    </tr>
                    <tr>
                        <td>Number of duplicates</td>
                        <td>{duplicates}</td>
                    </tr>
                </table>
                <p style="text-align: justify">
                    Duplicates occur when multiple reactions share the same enzyme-substrate combination. 
                    A high number of duplicates may result from multiple enzyme complexes sharing the same catalytic enzyme.
                </p>
            </div>

            <div class="card">
                <h2>Warning</h2>
                <table>
                    <tr>
                        <th>Metric</th>
                        <th>Value</th>
                    </tr>
                    <tr>
                        <td>Number of transport reactions</td>
                        <td>{predicted_transport}</td>
                    </tr>
                </table>
                <p style="text-align: justify">
                    Transport reaction represents the transfer of one metabolite from one cellular compartment to another without altering their chemical structure. 
                    The enzymes catalysing these transport reactions are often not represented in experimental databases (BRENDA and SABIO-RK). 
                    We argue that kcat  predictions for such reactions tend to be less reliable due to the limited availability of experimental data, which limits the performance of deep learning models. 
                    Consequently, these values should be interpreted with caution. 
                    The full list of the predicted transport reactions is available in the log file 'prediction2_YYYY-MM-DD_HH-MM-SS.log'.
                </p>
            </div>

            <!-- Prediction Instructions -->
            <div class="card">
                <h2>Running k<sub>cat</sub> Predictions with CataPro</h2>
                <p style="text-align: justify">
                    This report provides the input needed to run the CataPro machine learning model 
                    (<a href="https://github.com/zchwang/CataPro" target="_blank">CataPro repository</a>). 
                    Follow the instructions in the repository to set up the environment and generate k<sub>cat</sub> predictions.
                </p>
            </div>

    <footer>WILDkCAT</footer>
    """
    if shader:
        html += report_shader()
    else: 
        html += report_simple()
    html += """
    </body>
    </html>
    """

    # Save report
    os.makedirs(os.path.join(output_folder, "reports"), exist_ok=True)
    report_path = os.path.join(output_folder, "reports/predict_report.html")
    with open(report_path, "w", encoding="utf-8") as f:
        f.write(html)
    logging.info(f"HTML report saved to '{report_path}'")


def report_final(model, final_df, output_folder, shader=False) -> None:
    """
    Generate a full HTML report summarizing retrieval results, including kcat distributions and coverage.

    Parameters:
        model (cobra.Model): The metabolic model object containing reactions, metabolites, and genes.
        final_df (pd.DataFrame): DataFrame containing the final kcat assignments from run_prediction_part2 function
        
    Returns: 
        None: The function saves the generated HTML report to 'reports/general_report.html'.
    """
    # Model information 
    nb_model_reactions = len(model.reactions)
    nb_model_metabolites = len(model.metabolites)
    nb_model_genes = len(model.genes)


    df = final_df.copy()
    df["db"] = df["db"].fillna("Unknown")
    generated_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    # Utility to convert matplotlib figures to base64 <img>
    def fig_to_base64(fig):
        buf = BytesIO()
        fig.savefig(buf, format="png", bbox_inches="tight")
        buf.seek(0)
        encoded = base64.b64encode(buf.read()).decode("utf-8")
        plt.close(fig)
        return f'<div class="plot-container"><img src="data:image/png;base64,{encoded}"></div>'

    # Distribution plots
    def plot_kcat_distribution_stacked(column_name, title, source):
        # Ensure numeric kcat
        df[column_name] = pd.to_numeric(df[column_name], errors='coerce')

        # Drop NaNs for both columns
        valid_df = df.dropna(subset=[column_name, source])
        kcat_values = valid_df[column_name]

        total = len(df)
        matched = len(kcat_values)
        match_percent = matched / total * 100 if total else 0

        if not kcat_values.empty:
            # Define log bins
            min_exp = int(np.floor(np.log10(max(1e-6, kcat_values.min()))))
            max_exp = int(np.ceil(np.log10(kcat_values.max())))
            bins = np.logspace(min_exp, max_exp, num=40)

            # Prepare data for stacked histogram
            sources = valid_df[source].unique()
            grouped_values = [valid_df.loc[valid_df[source] == src, column_name] for src in sources]

            # Fixed color mapping
            color_map = {
                "brenda": "#55bb55",   
                "sabio_rk": "#2277cc", 
                # "catapro": "#eedd00",  
                "CataPro": "#eedd00",  
                "Unknown": "#dddddd" 
            }

            label_map = {
                "brenda": "Brenda",
                "sabio_rk": "Sabio-RK",
                # "catapro": "CataPro",
                "CataPro": "CataPro",
                "Unknown": "Unknown"
            }

            sources = [src for src in sources if src in valid_df[source].unique()]
            colors = [color_map.get(src, "#999999") for src in sources]

            # Plot
            fig, ax = plt.subplots(figsize=(12, 6))
            ax.hist(grouped_values, bins=bins, stacked=True,
                    color=colors, label=[label_map[s] for s in sources],
                    edgecolor="white", linewidth=0.7)

            ax.set_xscale("log")
            ax.set_xlim([10**min_exp / 1.5, 10**max_exp * 1.5])
            ax.xaxis.set_major_formatter(LogFormatter(10))
            ax.yaxis.set_major_locator(MaxNLocator(integer=True))

            ax.set_xlabel("kcat (s⁻¹)", fontsize=12)
            ax.set_ylabel("Count", fontsize=12)
            ax.set_title(f"{title} (n={matched}, {match_percent:.1f}%)", fontsize=13)

            ax.legend(
                title="Source", 
                fontsize=10, 
                title_fontsize=11,
                loc='center left', 
                bbox_to_anchor=(1, 0.5),
                frameon=False
            )

            # Style
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['left'].set_color('#444444')
            ax.spines['bottom'].set_color('#444444')

            ax.grid(True, which='major', axis='y', linestyle='--', linewidth=0.6, alpha=0.4)
            ax.grid(False, which='major', axis='x') 
            
            plt.tight_layout(rect=[0, 0, 0.85, 1])

            return fig_to_base64(fig)

        return "<p>No valid values available for plotting.</p>"
    
    img_final = plot_kcat_distribution_stacked(
        'kcat', rf"{model.id} - $k_{{\mathrm{{cat}}}}$ Distribution", "db"
    )
    
    db_counts = df["db"].fillna("Unknown").value_counts()
    total_db = db_counts.sum()

    # Colors
    colors = {
        "brenda": "#55bb55",
        "sabio_rk": "#2277cc",
        # "catapro": "#eedd00",
        "CataPro": "#eedd00",
        "Unknown": "#ddd"
    }

    # Order
    # ordered_dbs = ["brenda", "sabio_rk", "catapro", "Unknown"]
    ordered_dbs = ["brenda", "sabio_rk", "CataPro", "Unknown"]

    progress_segments = ""
    legend_items = ""

    for db in ordered_dbs:
        count = db_counts.get(db, 0)
        if total_db > 0:
            percent = count / total_db * 100
        else:
            percent = 0

        color = colors.get(db, "#ddd")

        progress_segments += f"""
            <div class="progress-segment" style="width:{percent:.1f}%; background-color:{color};"
                title="{db.capitalize()}: {percent:.1f}%"></div>
        """

        legend_items += f"""
            <span style="display:flex; align-items:center; margin-right:15px; margin-bottom:5px;">
                <span style="display:flex; align-items:center; width:16px; height:16px; 
                            background:{color}; border:1px solid #000; margin-right:5px;"></span>
                {db.capitalize()} ({percent:.1f}%)
            </span>
        """

    progress_bar = f"""
        <div class="progress-multi" style="height: 18px; margin-bottom:18px; display:flex;">
            {progress_segments}
        </div>
        <div style="margin-top:10px; display:flex; justify-content:center; flex-wrap: wrap;">
            {legend_items}
        </div>
    """

    # Statistics 
    grouped = df.groupby("rxn")
    rxns_with_kcat = grouped["kcat"].apply(lambda x: x.notna().any())
    nb_reactions = df['rxn'].nunique()
    nb_rxn_with_kcat = rxns_with_kcat.sum()
    coverage = nb_rxn_with_kcat / nb_reactions
    coverage_total = nb_rxn_with_kcat / nb_model_reactions

    kcat_values = df["kcat"].dropna()
    total = len(df)
    matched = len(kcat_values)
    match_percent = matched / total

    # HTML
    html = f"""
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <title>WILDkCAT Report</title>
        {report_style()}
    </head>
    <body>
        <header>
            <canvas id="shader-canvas"></canvas>
            <div class="overlay">
                <h1>WILDkCAT Report</h1>
                <p>Generated on {generated_time}</p>
            </div>
        </header>

        <div class="container">
            <div class="card">
                <h2>Introduction</h2>
                <p style="text-align: justify;">
                    This report provides a summary of the performance of k<sub>cat</sub> value extraction, retrieval, and prediction for the specified metabolic model. 
                    It presents statistics on k<sub>cat</sub> values successfully retrieved, whether experimental or predicted.
                </p>
            </div>

            <div class="card">
                <h2>Model Overview</h2>
                <div class="stats-grid">
                    <div class="stat-box">
                        <h3>{model.id}</h3>
                        <p>Model ID</p>
                    </div>
                    <div class="stat-box">
                        <h3>{nb_model_reactions}</h3>
                        <p>Reactions</p>
                    </div>
                    <div class="stat-box">
                        <h3>{nb_model_metabolites}</h3>
                        <p>Metabolites</p>
                    </div>
                    <div class="stat-box">
                        <h3>{nb_model_genes}</h3>
                        <p>Genes</p>
                    </div>
                </div>
            </div>

            <div class="card" style="padding:20px; margin-bottom:20px;">
                <h2 style="margin-bottom:10px;">Coverage</h2>
                
                <!-- Explanation -->
                <p style="text-align: justify;">
                    The coverage section reports the number of k<sub>cat</sub> values retrieved for the model and the number of reactions that have at least one 
                    associated k<sub>cat</sub> value. This provides a measure of how extensively the model’s reactions are 
                    annotated with kinetic data.
                </p>
                <p style="text-align: justify;">
                    Higher coverage indicates that a larger fraction of reactions are constrained by k<sub>cat</sub> values, 
                    improving the accuracy and reliability of enzyme-constrained simulations.
                </p>

                <!-- Global coverage progress bar -->
                {progress_bar}        

                <!-- Detailed stats -->
                <table class="table" style="width:100%; border-spacing:0; border-collapse: collapse;">
                    <tbody>
                    <tr>
                            <td style="padding:8px 12px;">Eligible-reactions with at least one kcat value</td>
                            <td style="padding:8px 12px;">{nb_rxn_with_kcat} ({coverage:.1%})</td>
                            <td style="width:40%;">
                                <div class="progress" style="height:18px;">
                                    <div class="progress-bar-table" 
                                        style="width:{coverage:.1%}; background-color:#4caf50;">
                                    </div>
                                </div>
                            </td>
                        </tr>
                        <tr>
                            <td style="padding:8px 12px;">Model-wide reactions with at least one kcat value</td>
                            <td style="padding:8px 12px;">{nb_rxn_with_kcat} ({coverage_total:.1%})</td>
                            <td style="width:40%;">
                                <div class="progress" style="height:18px;">
                                    <div class="progress-bar-table" 
                                        style="width:{coverage_total:.1%}; background-color:#4caf50;">
                                    </div>
                                </div>
                            </td>
                        </tr>
                        <tr>
                            <td style="padding:8px 12px;">k<sub>cat</sub> values retrieved </td>
                            <td style="padding:8px 12px;">{matched} ({match_percent:.1%})</td>
                            <td style="width:40%;">
                                <div class="progress" style="height:18px;">
                                    <div class="progress-bar-table" 
                                        style="width:{match_percent:.1%}; background-color:#4caf50;">
                                    </div>
                                </div>
                            </td>
                        </tr>
                    </tbody>
                </table>
            </div>
            
            <div class="card">
                <h2>k<sub>cat</sub> Distribution</h2>
                <div class="img-section">
                    {img_final}
                </div>
            </div>
        </div>

        <footer>WILDkCAT</footer>
    """
    if shader:
        html += report_shader()
    else: 
        html += report_simple()
    html += """
    </body>
    </html>
    """

    os.makedirs(os.path.join(output_folder, "reports"), exist_ok=True)
    report_path = os.path.join(output_folder, "reports/general_report.html")
    with open(report_path, "w", encoding="utf-8") as f:
        f.write(html)

    logging.info(f"HTML report saved to '{report_path}'")
    return report_path


def report_style():
    """Return CSS script for report style."""
    return """
    <style>
        body {
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            font-size: 1rem;
            line-height: 1.7;
            padding: 2rem 1.5rem;
            background-color: #f4f6f9;
            margin: 0;
            padding: 0;
            color: #333;
        }
        header {
            position: relative;
            width: 100%;
            height: 150px;
            overflow: hidden;
            display: flex;
            align-items: center;
            justify-content: center;
            color: #fff;
            text-align: center;
        }
        header canvas {
            position: absolute;
            top: 0; left: 0;
            width: 100%;
            height: 100%;
            z-index: 0;
        }
        header::before {
            content: "";
            position: absolute;
            top: 0; left: 0; right: 0; bottom: 0;
            background: linear-gradient(
                rgba(0,0,0,0.5),
                rgba(0,0,0,0.3)
            );
            z-index: 1;
        }
        header .overlay {
            position: relative;
            z-index: 2;
            padding: 10px 20px;
            border-radius: 8px;
        }
        header h1 {
            margin: 0;
            font-size: 2.5rem;
            font-weight: bold;
            text-shadow: 0 2px 6px rgba(0,0,0,0.6);
        }
        header p {
            margin: 8px 0 0;
            font-size: 1.1rem;
            text-shadow: 0 1px 4px rgba(0,0,0,0.6);
        }
        p {
            margin-bottom: 1.2rem;
        }
        .container {
            max-width: 1100px;
            margin: 30px auto;
            padding: 20px;
        }
        .card {
            background: #fff;
            border-radius: 12px;
            padding: 20px;
            margin-bottom: 20px;
            box-shadow: 0 2px 8px rgba(0,0,0,0.05);
        }
        .card h2 {
            margin-top: 0;
            color: #2980b9;
            border-bottom: 2px solid #e6e6e6;
            padding-bottom: 10px;
            font-size: 1.5rem;
        }
        .stats-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 15px;
            margin-top: 15px;
        }
        .stat-box {
            background: #f9fafc;
            border-radius: 8px;
            padding: 15px;
            text-align: center;
            border: 1px solid #e2e2e2;
        }
        table {
            width: 100%;
            border-collapse: collapse;
            margin-top: 20px;
            font-size: 0.95rem;
        }
        table th, table td {
            border: 1px solid #ddd;
            padding: 10px;
            text-align: left;
        }
        table th {
            background-color: #2980b9;
            color: #fff;
        }
        table tr:nth-child(even) {
            background-color: #f2f2f2;
        }
        .progress {
            background-color: #ddd;
            border-radius: 10px;
            overflow: hidden;
            height: 18px;
            width: 100%;
            margin-top: 5px;
        }
        .progress-stacked {
            display: flex;
            height: 18px;
            border-radius: 10px;
            overflow: hidden;
            background-color: #ddd;
            font-size: 0.75rem;
            line-height: 18px;
            color: white;
            text-shadow: 0 1px 1px rgba(0,0,0,0.2);
            margin-bottom: 10px;
        }
        .progress-bar {
            display: flex;
            align-items: center;
            justify-content: center;
            height: 100%;
            white-space: nowrap;
            overflow: hidden;
        }
        .progress-bar-table {
            background-color: #27ae60;
            height: 100%;
            text-align: right;
            padding-right: 5px;
            color: white;
            font-size: 0.8rem;
            line-height: 18px;
        }
        .progress-multi {
            display: flex;
            width: 100%;
            height: 25px;
            border-radius: 12px;
            overflow: hidden;
            border: 1px solid #ccc;
        }
        .progress-segment {
            height: 100%;
        }
        .legend {
            display: flex;
            flex-wrap: wrap;
            gap: 10px;
            font-size: 0.85rem;
            margin-top: 5px;
        }
        .legend-item {
            display: flex;
            align-items: center;
            gap: 5px;
        }
        .legend-color {
            width: 14px;
            height: 14px;
            border-radius: 3px;
            border: 1px solid #aaa;
        }
        .img-section {
            display: flex;
            flex-wrap: wrap;
            gap: 30px;
            justify-content: center;
            align-items: flex-start;
            margin-top: 20px;
        }
        footer {
            text-align: center;
            font-size: 0.9rem;
            color: #777;
            padding: 15px;
            margin-top: 20px;
            border-top: 1px solid #ddd;
        }
    </style>
    """


def report_simple():
    """Return HTML code for report background."""
    return """
    <style>
        header {
            background-color: #2980b9; /* simple blue background */
            margin: 0;
            padding: 0;
        }
    </style>
    """


def report_shader(): 
    """Return HTML and GLSL shader code for report background. Adapted from localthunk (https://localthunk.com)"""
    return """
    <!-- Background adapted from original work by localthunk (https://localthunk.com) -->
    <script id="fragShader" type="x-shader/x-fragment">
    precision highp float;
    uniform vec2 iResolution;
    uniform float iTime;
    #define SPIN_ROTATION -1.0
    #define SPIN_SPEED 3.5
    #define OFFSET vec2(0.0)
    #define COLOUR_1 vec4(0.2, 0.4, 0.7, 1.0)
    #define COLOUR_2 vec4(0.6, 0.75, 0.9, 1.0)
    #define COLOUR_3 vec4(0.2, 0.2, 0.25, 1.0)
    #define CONTRAST 3.5
    #define LIGTHING 0.4
    #define SPIN_AMOUNT 0.25
    #define PIXEL_FILTER 745.0
    #define SPIN_EASE 1.0
    #define PI 3.14159265359
    #define IS_ROTATE false
    vec4 effect(vec2 screenSize, vec2 screen_coords) {
        float pixel_size = length(screenSize.xy) / PIXEL_FILTER;
        vec2 uv = (floor(screen_coords.xy*(1./pixel_size))*pixel_size - 0.5*screenSize.xy)/length(screenSize.xy) - OFFSET;
        float uv_len = length(uv);
        float speed = (SPIN_ROTATION*SPIN_EASE*0.2);
        if(IS_ROTATE) {
        speed = iTime * speed;
        }
        speed += 302.2;
        float new_pixel_angle = atan(uv.y, uv.x) + speed - SPIN_EASE*20.*(1.*SPIN_AMOUNT*uv_len + (1. - 1.*SPIN_AMOUNT));
        vec2 mid = (screenSize.xy/length(screenSize.xy))/2.;
        uv = (vec2((uv_len * cos(new_pixel_angle) + mid.x), (uv_len * sin(new_pixel_angle) + mid.y)) - mid);
        uv *= 30.;
        speed = iTime*(SPIN_SPEED);
        vec2 uv2 = vec2(uv.x+uv.y);
        for(int i=0; i < 5; i++) {
            uv2 += sin(max(uv.x, uv.y)) + uv;
            uv  += 0.5*vec2(cos(5.1123314 + 0.353*uv2.y + speed*0.131121),sin(uv2.x - 0.113*speed));
            uv  -= 1.0*cos(uv.x + uv.y) - 1.0*sin(uv.x*0.711 - uv.y);
        }
        float contrast_mod = (0.25*CONTRAST + 0.5*SPIN_AMOUNT + 1.2);
        float paint_res = min(2., max(0.,length(uv)*(0.035)*contrast_mod));
        float c1p = max(0.,1. - contrast_mod*abs(1.-paint_res));
        float c2p = max(0.,1. - contrast_mod*abs(paint_res));
        float c3p = 1. - min(1., c1p + c2p);
        float light = (LIGTHING - 0.2)*max(c1p*5. - 4., 0.) + LIGTHING*max(c2p*5. - 4., 0.);
        return (0.3/CONTRAST)*COLOUR_1 + (1. - 0.3/CONTRAST)*(COLOUR_1*c1p + COLOUR_2*c2p + vec4(c3p*COLOUR_3.rgb, c3p*COLOUR_1.a)) + light;
    }
    void mainImage(out vec4 fragColor, in vec2 fragCoord) {
        vec2 uv = fragCoord/iResolution.xy;
        fragColor = effect(iResolution.xy, uv * iResolution.xy);
    }
    void main() { mainImage(gl_FragColor, gl_FragCoord.xy); }
    </script>
    <script>
    const canvas = document.getElementById("shader-canvas");
    const gl = canvas.getContext("webgl");
    function resize() {
        canvas.width = canvas.clientWidth * window.devicePixelRatio;
        canvas.height = canvas.clientHeight * window.devicePixelRatio;
        gl.viewport(0, 0, canvas.width, canvas.height);
    }
    window.addEventListener("resize", resize);
    resize();
    const vertexSrc = `
    attribute vec2 position;
    void main() {
        gl_Position = vec4(position, 0.0, 1.0);
    }
    `;
    const fragSrc = document.getElementById("fragShader").text;
    function compileShader(src, type) {
    const shader = gl.createShader(type);
    gl.shaderSource(shader, src);
    gl.compileShader(shader);
    if (!gl.getShaderParameter(shader, gl.COMPILE_STATUS)) {
        console.error(gl.getShaderInfoLog(shader));
    }
    return shader;
    }
    const vertexShader = compileShader(vertexSrc, gl.VERTEX_SHADER);
    const fragmentShader = compileShader(fragSrc, gl.FRAGMENT_SHADER);
    const program = gl.createProgram();
    gl.attachShader(program, vertexShader);
    gl.attachShader(program, fragmentShader);
    gl.linkProgram(program);
    gl.useProgram(program);
    const positionBuffer = gl.createBuffer();
    gl.bindBuffer(gl.ARRAY_BUFFER, positionBuffer);
    gl.bufferData(gl.ARRAY_BUFFER, new Float32Array([
    -1, -1, 1, -1, -1, 1,
    -1, 1, 1, -1, 1, 1
    ]), gl.STATIC_DRAW);
    const positionLoc = gl.getAttribLocation(program, "position");
    gl.enableVertexAttribArray(positionLoc);
    gl.vertexAttribPointer(positionLoc, 2, gl.FLOAT, false, 0, 0);
    const iResolutionLoc = gl.getUniformLocation(program, "iResolution");
    const iTimeLoc = gl.getUniformLocation(program, "iTime");
    function render(time) {
    resize();
    gl.uniform2f(iResolutionLoc, canvas.width, canvas.height);
    gl.uniform1f(iTimeLoc, time * 0.001);
    gl.drawArrays(gl.TRIANGLES, 0, 6);
    requestAnimationFrame(render);
    }
    requestAnimationFrame(render);
    </script>
    """