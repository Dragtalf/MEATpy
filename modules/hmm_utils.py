import os
import re
import logging
import pandas as pd
from dataclasses import dataclass

@dataclass(repr=True)
class HMMThreshold:
    threshold: float
    score_type: str  # "full" or "domain"

    def __post_init__(self):
        if self.score_type not in ("full", "domain"):
            raise ValueError(f"Invalid score_type: {self.score_type} (must be 'full' or 'domain')")

class HMM_Manager:
    def __init__(self, paths, kofam_db_size):
        self.paths = paths
        self.kofam_db_size = kofam_db_size

    def load_hmm_templates(self):
        # Read hmm_table into DataFrame, skip comment lines but capture header
        with open(os.path.join(self.paths['templates'], 'hmm_table_template.txt'), 'r') as f:
            lines = f.readlines()
        header_line = next(line for line in lines if line.startswith('#'))
        hmm_table_head = header_line.strip().lstrip('#').split('\t')
        
        data_lines = [line for line in lines if not line.startswith('#')]
        hmm_df = pd.DataFrame([line.strip().split('\t') for line in data_lines], columns=hmm_table_head)

        # Dict[line_no] = original line string
        hmm_table_temp = {
            row[0]: "\t".join(row) for row in hmm_df.values.tolist()
        }

        # Build hmm2threshold for non-Kegg IDs
        hmm_thresholds = {}
        for _, row in hmm_df.iterrows():
            hmm_id = row.iloc[10]
            threshold = row.iloc[11]
            if hmm_id and not re.match(r"K\d{5}$", hmm_id):  # Not a Kegg ID
                hmm_thresholds[hmm_id] = threshold

        # Read the second table (template 2)
        hmm_df2 = pd.read_csv(
            os.path.join(self.paths['templates'], 'hmm_table_template_2.txt'),
            sep='\t', comment='#', header=None
        )
        hmm_table_temp_2 = {
            row[0]: "\t".join(map(str, row)) for row in hmm_df2.values.tolist()
        }

        return hmm_table_temp, hmm_table_head, hmm_thresholds, hmm_table_temp_2, hmm_df, hmm_df2

    def get_kofam_db_ko_threshold(self):
        # Determine the path to the prokaryote list based on the kofam_db_size
        logging.info(f" kofam db_size: {self.kofam_db_size}")        
        prok_list = os.path.join(self.paths['databases'], 'kofam_database', 'profiles', 'prokaryote.hal' if self.kofam_db_size == 'full' else 'All_Module_KO_ids.txt')
        kofam_thresholds_raw = {}

        # Read the ko_list file and populate the result dictionary
        try:
            with open(f"{self.paths['databases']}/kofam_database/ko_list", 'r') as file:
                for line in file:
                    if line.startswith('K'):
                        parts = line.strip().split('\t')
                        if len(parts) >= 3:
                            hmm_id = f"{parts[0]}.hmm"
                            if parts[1] == '-':
                                kofam_thresholds_raw[hmm_id] = (50.0, "full")
                            else:
                                try:
                                    threshold = float(parts[1])
                                    kofam_thresholds_raw[hmm_id] = (threshold, parts[2])
                                except ValueError:
                                    logging.warning(f"Skipping line with invalid number: {line.strip()}")
                        else:
                            logging.warning(f"Skipping malformed line: {line.strip()}")
        except FileNotFoundError:
            logging.error(f"File not found: {self.paths['databases']}/kofam_database/ko_list")
            raise

        # Read the prok_list file and populate the prok_list_entries set
        try:
            with open(prok_list, 'r') as file:
                prok_list_entries = {line.strip() for line in file}
        except FileNotFoundError:
            logging.error(f"File not found: {prok_list}")
            raise

        # Filter the result dictionary to only include entries present in prok_list_entries
        kofam_thresholds = {key: val for key, val in kofam_thresholds_raw.items() if key in prok_list_entries}

        return kofam_thresholds

    def merge_thresholds(self, hmm_thresholds, kofam_thresholds):
        total_hmm2threshold = {}

        # Handle thresholds
        for hmm_ids, value in hmm_thresholds.items():
            id_list = [h.strip() for h in re.split(r'[;,]', hmm_ids)]
            values = [v.strip() for v in value.split(',')]

            if len(id_list) != len(values):
                logging.warning(f"⚠️ Mismatched count for {hmm_ids}: HMMs={len(id_list)} vs thresholds={len(values)}")

            for hmm, val in zip(id_list, values):
                try:
                    threshold_str, score_type = val.split('|')
                    total_hmm2threshold[hmm] = HMMThreshold(float(threshold_str), score_type)
                except ValueError:
                    logging.error(f"⚠️ Skipping {hmm}: malformed threshold string '{val}'")
        
        # Handle Kofam thresholds
        for hmm_id, (threshold, score_type) in kofam_thresholds.items():
            if hmm_id in total_hmm2threshold:
                logging.debug(f"Overwriting curated threshold for {hmm_id} with KOfam entry")
            total_hmm2threshold[hmm_id] = HMMThreshold(threshold, score_type)

        return total_hmm2threshold