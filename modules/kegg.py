import os
import logging
import re
import json
from dataclasses import dataclass
from collections import defaultdict
from .functions import get_hmm_to_ko_hash

class KEGGAnnotator:
    def __init__(self, paths, hmm_table_temp, genome_ids, hmmscan_result, module_cutoff):
        self.paths = paths
        self.hmm_table_temp = hmm_table_temp
        self.genome_ids = genome_ids
        self.hmmscan_result = hmmscan_result
        self.module_cutoff = module_cutoff
        self.kegg_module = {}  # module_id -> list of KOs (from expression parsing)
        self.category_modules = self._load_module_categories()
        self.module_to_kos = {}  # module -> list of KOs (from expressions)
        self.kegg_module2name = {}  # module -> description
        self.module_presence_by_genome = {}
        self.step_presence_by_genome = {}  # step = module + KO, for worksheet 4
        self.hmm2ko = {}

    def run(self):
        logging.info("KEGG annotation pipeline...")
        self._parse_module_step_file()
        self.category_modules = self._load_module_categories()
        self._analyze_modules_by_genome()
        logging.info("KEGG annotation complete.")

        return (
            self.category_modules,
            self.kegg_module,
            self.kegg_module2name,
            self.module_presence_by_genome,
            self.step_presence_by_genome,
            self.hmm2ko
        )

    def _parse_module_step_file(self):
        logging.info("Parsing KEGG module step database...")
        step_db_path = os.path.join(self.paths['templates'], 'kegg_module_step_db.txt')

        category_modules = defaultdict(list)
        module_to_kos = defaultdict(list)
        kegg_module2name = {}
        kegg_module = {}

        with open(step_db_path, 'r') as file:
            next(file)  # Skip header
            for line in file:
                parts = line.strip().split('\t')
                if len(parts) != 3:
                    continue
                name, k_string, module_step = parts
                module, step_number = module_step.split('+')
                kegg_module[module_step] = [k_string, int(step_number)]
                kegg_module2name[module] = name
                category_modules[name].append(module)
                module_to_kos[module].append(k_string)

        self.kegg_module = kegg_module
        self.kegg_module2name = kegg_module2name
        self.category_modules = category_modules
        self.module_to_kos = module_to_kos
        logging.info("KEGG module step database parsed successfully.")

    def _load_module_categories(self):
        logging.info("Loading KEGG module categories from JSON...")
        category_modules = defaultdict(list)

        with open(os.path.join(self.paths['templates'], 'ko00002.json'), 'r') as f:
            data = json.load(f)

        def extract_modules(node, parent_categories=None):
            if parent_categories is None:
                parent_categories = []

            name = node.get("name", "")
            new_parents = parent_categories + [name] if 'children' in node else parent_categories

            if 'children' in node:
                for child in node['children']:
                    extract_modules(child, new_parents)
            elif name.startswith("M"):
                module_id = name.split()[0]
                if len(parent_categories) >= 3:
                    # Use 3nd-level category (e.g., "Central carbohydrate metabolism")
                    category = parent_categories[2]
                else:
                    # Fallback to higher category
                    category = parent_categories[1]
                category_modules[category].append(module_id)

        for top_level in data['children']:
            extract_modules(top_level)

        logging.info("KEGG module categories loaded successfully from JSON.")
        return category_modules

    def _analyze_modules_by_genome(self):
        logging.info("Analyzing modules by genome...")
        hmm2ko = get_hmm_to_ko_hash(self.hmm_table_temp)
        
        self.hmm2ko = hmm2ko

        step_presence_by_genome = {}
        module_to_steps = defaultdict(list)

        for m_step, (k_string, step_num) in self.kegg_module.items():
            module_id = m_step.split('+')[0]
            module_to_steps[module_id].append(step_num)
            for genome in self.genome_ids:
                genome_hmms = self.hmmscan_result.get(genome, {})
                ko_hits = [hmm2ko.get(hmm, hmm).replace('.hmm', '') for hmm in genome_hmms]
                val = self.determine_module_step(k_string, ko_hits)
                step_presence_by_genome.setdefault(m_step, {})[genome] = val

        module_presence_by_genome = {}
        for module, step_nums in module_to_steps.items():
            step_total = len(step_nums)
            for genome in self.genome_ids:
                present_count = sum(
                    int(step_presence_by_genome.get(f"{module}+{str(step).zfill(2)}", {}).get(genome, 0))
                    for step in step_nums
                )
                ratio = present_count / step_total
                module_presence_by_genome.setdefault(module, {})[genome] = (
                    "Present" if ratio >= self.module_cutoff else "Absent"
                )

        self.step_presence_by_genome = step_presence_by_genome
        self.module_presence_by_genome = module_presence_by_genome
        logging.info("Module analysis by genome completed.")

    @staticmethod
    def determine_module_step(k_string, ko_hits):
        ko_hash = set(k.replace('.hmm', '') for k in ko_hits)
        result = k_string

        for ko_number in re.findall(r'K\d{5}', k_string):
            result = result.replace(ko_number, '1' if ko_number in ko_hash else '0')

        safe_expr = (
            result.replace('&', ' and ')
                  .replace('|', ' or ')
                  .replace('!', ' not ')
                  .replace('KXXXXX', '0')
        )

        try:
            evaluation_result = eval(safe_expr, {"__builtins__": None}, {})
            return '1' if evaluation_result else '0'
        except Exception as e:
            if k_string == 'KXXXXX':
                return '0'
            else:
                logging.error(f"❌ Eval failed: {k_string} | Error: {e} | Result: {result}")
                return '0'
        
class KEGGProcessor:
    def __init__(self, paths, hmmscan_result, hmmscan_hits, hmm2ko, genome_ids, kofam_db_size):
        self.paths = paths
        self.hmmscan_result = hmmscan_result
        self.hmmscan_hits = hmmscan_hits
        self.hmm2ko = hmm2ko
        self.genome_ids = genome_ids
        self.kofam_db_size = kofam_db_size
        
    def run(self):
        self.process_kegg_hits()
        
    def process_kegg_hits(self):
        """
        Process HMM scan results and map to KO IDs, writing result and hit files with all KOs listed.

        Args:
            hmmscan_result (dict): {genome_id → {hmm → count}}
            hmmscan_hits (dict): {genome_id → {hmm → hits}}
            hmm2ko (dict): {hmm_id → ko_id.hmm}
            output_dir (str): Path to output directory
            genome_ids (list or set): List of genome IDs
            kofam_db_path (str): Path to Kofam database folder
            kofam_db_size (str): 'full' or 'small'
        """
        
        kegg_dir = os.path.join(self.paths['output_dir'], "KEGG_identifier_result")
        os.makedirs(kegg_dir, exist_ok=True)

        # Load full KO list from prokaryote.hal or All_Module_KO_ids.txt
        ko_list_file = os.path.join(self.paths['databases'], 'kofam_database', 'profiles', 'prokaryote.hal' if self.kofam_db_size == 'full' else 'All_Module_KO_ids.txt')
        try:
            with open(ko_list_file, 'r') as f:
                full_ko_set = sorted(line.strip().replace('.hmm', '') for line in f if line.strip().endswith('.hmm'))
        except Exception as e:
            raise RuntimeError(f"Failed to read KO list from {ko_list_file}: {e}")

        # Map each genome's HMM results to KO results
        ko_results_by_genome = {}
        ko_hits_by_genome = {}

        for genome in self.genome_ids:
            ko_results_by_genome[genome] = {}
            ko_hits_by_genome[genome] = {}

            for hmm in self.hmmscan_result.get(genome, {}):
                if hmm in self.hmm2ko:
                    ko = self.hmm2ko[hmm]
                elif re.match(r'^K\d{5}(\.hmm)?$', hmm):
                    ko = hmm if hmm.endswith('.hmm') else f"{hmm}.hmm"
                else:
                    continue
                ko_id = ko.replace('.hmm', '')
                ko_results_by_genome[genome][ko_id] = self.hmmscan_result[genome][hmm]
                ko_hits_by_genome[genome][ko_id] = self.hmmscan_hits.get(genome, {}).get(hmm, "")

        # Write out all KOs per genome, regardless of presence
        for genome in sorted(self.genome_ids):
            result_path = os.path.join(kegg_dir, f"{genome}.result.txt")
            hits_path = os.path.join(kegg_dir, f"{genome}.hits.txt")

            with open(result_path, 'w') as out_r, open(hits_path, 'w') as out_h:
                for ko_id in full_ko_set:
                    count = ko_results_by_genome[genome].get(ko_id, "")
                    hits = ko_hits_by_genome[genome].get(ko_id, "")
                    out_r.write(f"{ko_id}\t{count}\n")
                    out_h.write(f"{ko_id}\t{hits}\n")
