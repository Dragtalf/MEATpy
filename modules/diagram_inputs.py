import os
from collections import defaultdict
from statistics import mean
import re
import logging


class DiagramInputs:
    def __init__(self, paths, hmmscan_result, dbcan_results, threads, correction, user_tax, input_csv_check):
        self.paths = paths
        self.hmmscan_result = hmmscan_result
        self.dbcan_results = dbcan_results
        self.threads = threads
        self.correction = correction
        self.user_tax = user_tax
        self.input_csv_check = input_csv_check
        self.total_r_input = None
        self.genome_coverage = None
        self.bin2cat = None
        self.fig_input_dir = os.path.join(self.paths['output_dir'], 'MEATpy_Figures_Input')
        
    TAXONOMIC_LEVELS = {
    'phylum': 0, 'class': 1, 'order': 2,
    'family': 3, 'genus': 4, 'species': 5, 'bin': 6
    }
    
    def run(self):
        self.generate_nutrient_cycling_inputs()
        self.generate_coverage_inputs()
        self.generate_metabolic_handoff_summary_step1()
        self.generate_metabolic_handoff_summary_step2()
        self.generate_energy_flow_data()
        self.calculate_mw_score()
        self.finalize_mw_score()
    
    def generate_nutrient_cycling_inputs(self):
        os.makedirs(os.path.join(self.paths['output_dir'], "MEATpy_Figures_Input", "Nutrient_Cycling_Diagram_Input"), exist_ok=True)

        # Parse R pathways
        r_pathways = {}  # step -> hmm logic string
        r_hmm_ids = set()
        with open(os.path.join(self.paths['templates'], "R_pathways.txt")) as f:
            for line in f:
                step, hmm_logic = line.strip().split('\t')
                r_pathways[step] = hmm_logic
                if ';' in hmm_logic:
                    for clause in hmm_logic.split(';'):
                        for hmm in clause.split(','):
                            if 'NO' not in hmm:
                                r_hmm_ids.add(hmm)
                else:
                    for hmm in hmm_logic.split(','):
                        r_hmm_ids.add(hmm)

        # Analyze presence of pathways per genome
        total_r_input = defaultdict(dict)
        for genome, genome_hmms in self.hmmscan_result.items():
            r_input = {}
            for step, logic in r_pathways.items():
                value = 0
                if ';' not in logic:
                    for hmm in logic.split(','):
                        if hmm in genome_hmms:
                            value = 1
                            break
                else:
                    clause1, clause2 = logic.split(';')
                    clause1_present = any(h in genome_hmms for h in clause1.split(',') if 'NO' not in h)
                    clause2_present = any(h in genome_hmms for h in clause2.split(',') if 'NO' not in h)
                    clause2_negated = any(h in clause2 for h in genome_hmms if 'NO' in clause2)
                    if 'NO' in clause2:
                        if clause1_present and not clause2_negated:
                            value = 1
                    else:
                        if clause1_present and clause2_present:
                            value = 1
                r_input[step] = value
                total_r_input[step][genome] = value

            # Write genome-specific file
            genome_path = os.path.join(self.paths['output_dir'], "MEATpy_Figures_Input", "Nutrient_Cycling_Diagram_Input", f"{genome}.R_input.txt")
            with open(genome_path, 'w') as f_out:
                for step in sorted(r_input):
                    f_out.write(f"{step}\t{r_input[step]}\n")

        self.total_r_input = total_r_input
        return total_r_input

    def generate_coverage_inputs(self):
        logging.info("Parsing coverage values...")
        depth_file = os.path.join(self.paths['output_dir'], 'coverage', 'All_gene_collections_mapped.depth.txt')
        if not os.path.exists(depth_file):
            raise FileNotFoundError("Coverage output not found. Check if coverm ran successfully.")

        coverage_data = {}
        total_coverage = 0

        with open(depth_file) as f:
            header = f.readline().strip().split("\t")
            avg_idx = header.index("totalAvgDepth")
            for line in f:
                fields = line.strip().split("\t")
                contig_name = fields[0]
                genome = contig_name.split("~~")[0]
                depth = float(fields[avg_idx])
                if genome not in coverage_data:
                    coverage_data[genome] = []
                coverage_data[genome].append(depth)

        genome_coverage = {}
        for genome, depths in coverage_data.items():
            avg_cov = mean(depths)
            genome_coverage[genome] = avg_cov
            total_coverage += avg_cov

        if total_coverage == 0:
            raise ValueError("Total genome coverage is zero. Something went wrong in mapping.")

        genome_coverage_percent = {
            genome: cov / total_coverage
            for genome, cov in genome_coverage.items()
        }
        logging.info("Coverage values parsed.")

        total_txt = os.path.join(self.fig_input_dir, 'Nutrient_Cycling_Diagram_Input', 'Total.R_input.txt')
        os.makedirs(os.path.dirname(total_txt), exist_ok=True)

        with open(total_txt, "w") as out_f:
            for pth in sorted(self.total_r_input.keys()):
                gn_no = 0
                gn_cov_percent = 0
                for gn_id in sorted(self.hmmscan_result.keys()):
                    if self.total_r_input[pth].get(gn_id, 0):
                        gn_no += 1
                        gn_cov_percent += genome_coverage_percent.get(gn_id, 0)
                out_f.write(f"{pth}\t{gn_no}\t{gn_cov_percent:.6f}\n")

        self.genome_coverage = genome_coverage_percent
        return genome_coverage_percent

    def generate_metabolic_handoff_summary_step1(self):
        r_mh_01 = {}
        hmm_ids = set()
        letter_to_rxn = {}

        with open(os.path.join(self.paths['templates'], 'Sequential_transformations_01.txt')) as infile:
            for line in infile:
                line = line.strip()
                if ":" in line or "+" in line:
                    parts = line.split("\t")
                    if ":" in parts[0]:
                        step, reaction = parts[0].split(":", 1)
                        r_mh_01[step] = parts[2]
                        letter_to_rxn[step] = reaction
                        for hmm in parts[2].split(","):
                            hmm_ids.add(hmm)
                    else:
                        r_mh_01[parts[0]] = parts[2]

        r_mh_01_summary = {}
        for genome in sorted(self.hmmscan_result):
            for step in sorted(r_mh_01):
                r_mh_01_summary.setdefault(step, {})[genome] = 0
                step_hmms = r_mh_01[step]

                if "+" not in step:
                    for hmm in hmm_ids:
                        if hmm in step_hmms and self.hmmscan_result[genome].get(hmm):
                            r_mh_01_summary[step][genome] = 1
                            break
                else:
                    group_hits = 0
                    hmm_groups = step_hmms.split(";")
                    for group in hmm_groups:
                        found = False
                        for hmm in hmm_ids:
                            if hmm in group and self.hmmscan_result[genome].get(hmm):
                                found = True
                                break
                        if found:
                            group_hits += 1
                    if group_hits == len(hmm_groups):
                        r_mh_01_summary[step][genome] = 1

        total_r_hm_input_1 = {}
        for step in sorted(r_mh_01):
            genome_count = 0
            coverage_sum = 0.0
            for genome in sorted(self.hmmscan_result):
                if r_mh_01_summary[step].get(genome, 0):
                    genome_count += 1
                    coverage_sum += self.genome_coverage.get(genome, 0.0)
            total_r_hm_input_1[step] = (genome_count, coverage_sum)

        out_path = os.path.join(self.fig_input_dir, 'Sequential_Transformation_input_1.txt')
        os.makedirs(os.path.dirname(out_path), exist_ok=True)

        with open(out_path, "w") as out_f:
            for step, (count, cov_sum) in sorted(total_r_hm_input_1.items()):
                out_f.write(f"{step}\t{count}\t{cov_sum:.6f}\n")

        return r_mh_01_summary

    def generate_metabolic_handoff_summary_step2(self):
        cazy_map = {}
        with open(os.path.join(self.paths['templates'], 'CAZy_map.txt')) as f:
            for line in f:
                if not line.startswith("Family"):
                    parts = line.strip().split("\t")
                    gh, enzyme = parts[0], parts[1]
                    cazy_map[gh] = enzyme

        r_mh_02 = {}
        r_mh_02_enzyme_ids = set()
        letter_to_reaction = {}
        with open(os.path.join(self.paths['templates'], 'Sequential_transformations_02.txt')) as f:
            for line in f:
                if ":" in line or "+" in line:
                    parts = line.strip().split("\t")
                    if ":" in parts[0]:
                        step, reaction = parts[0].split(":")
                        r_mh_02[step] = parts[2]
                        letter_to_reaction[step] = reaction
                        for enz in parts[2].split("; "):
                            r_mh_02_enzyme_ids.add(enz)
                    else:
                        r_mh_02[parts[0]] = parts[2]

        r_mh_02_summary = defaultdict(lambda: defaultdict(int))

        for step in sorted(r_mh_02):
            for gn in sorted(self.dbcan_results):
                r_mh_02_summary[step][gn] = 0
                if "+" not in step:
                    enzymes = r_mh_02[step].split(";")
                    for enzyme in enzymes:
                        for gh, mapped_enzymes in cazy_map.items():
                            if enzyme in mapped_enzymes.split(";") and self.dbcan_results[gn].get(gh):
                                r_mh_02_summary[step][gn] = 1
                else:
                    hmm_blocks = r_mh_02[step].split("|")
                    total_required = len(hmm_blocks)
                    found_count = 0
                    for hmm_block in hmm_blocks:
                        logic_flag = False
                        for enzyme in hmm_block.split(";"):
                            for gh, mapped_enzymes in cazy_map.items():
                                if enzyme in mapped_enzymes.split(";") and self.dbcan_results[gn].get(gh):
                                    logic_flag = True
                        if logic_flag:
                            found_count += 1
                    if found_count == total_required:
                        r_mh_02_summary[step][gn] = 1

        total_r_hm_input_2 = {}
        for step in sorted(r_mh_02):
            gn_no = 0
            gn_cov_percent = 0
            for gn in sorted(self.hmmscan_result):
                if r_mh_02_summary[step][gn]:
                    gn_no += 1
                    gn_cov_percent += self.genome_coverage.get(gn, 0)
            total_r_hm_input_2[step] = (gn_no, gn_cov_percent)

        out_path = os.path.join(self.fig_input_dir, 'Sequential_Transformation_input_2.txt')
        os.makedirs(os.path.dirname(out_path), exist_ok=True)

        with open(out_path, "w") as f:
            for step in sorted(total_r_hm_input_2):
                gn_no, cov = total_r_hm_input_2[step]
                f.write(f"{step}\t{gn_no}\t{cov:.6f}\n")

        return r_mh_02_summary
        
    def generate_energy_flow_data(self):
        
        gtdbtk_dir = os.path.join(self.paths['output_dir'], 'intermediate_files', 'gtdbtk_Genome_files')        

        self.bin2cat = {}

        if self.correction == "NCBI":

            self._parse_gtdb_summary(os.path.join(gtdbtk_dir, "ncbi.bac120.correction.tsv"))
            self._parse_gtdb_summary(os.path.join(gtdbtk_dir, "ncbi.ar53.correction.tsv"))
            logging.info("Parsing GTDB-Tk classification completed with NCBI correction.")

        else:
            self._parse_gtdb_summary(os.path.join(gtdbtk_dir, "gtdbtk.bac120.summary.tsv"))
            self._parse_gtdb_summary(os.path.join(gtdbtk_dir, "gtdbtk.ar53.summary.tsv"))
            logging.info("Parsing GTDB-Tk classification completed.")

        logging.info("Compiling energy flow input tables...")

        if self.user_tax not in DiagramInputs.TAXONOMIC_LEVELS:
            raise ValueError(f"Invalid taxonomy level '{self.user_tax}'. Must be one of {list(DiagramInputs.TAXONOMIC_LEVELS.keys())}")

        hash_gn_path = {}
        total_r_community_cov = {}
        for pth in sorted(self.total_r_input):
            for gn in sorted(self.hmmscan_result):
                if self.genome_coverage.get(gn) and self.total_r_input[pth].get(gn):
                    cat = self.bin2cat.get(gn, ["Unknown"] * 7)[DiagramInputs.TAXONOMIC_LEVELS[self.user_tax]]
                    gn_path_key = f"{gn}\t{pth}"
                    hash_gn_path[gn_path_key] = 1
                    total_r_community_cov[gn_path_key] = f"{cat}\t{pth}\t{self.genome_coverage[gn]:.6f}"

        total_r_community_cov2 = {}
        for gn in sorted(self.hmmscan_result):
            path_keys = [k.split('\t')[1] for k in hash_gn_path if k.startswith(f"{gn}\t")]
            for i in range(len(path_keys)):
                for j in range(i + 1, len(path_keys)):
                    step1, step2 = path_keys[i], path_keys[j]
                    key = f"{gn}\t{step1}\t{step2}"
                    cov1 = float(total_r_community_cov.get(f"{gn}\t{step1}", "0\t0\t0").split('\t')[2])
                    cov2 = float(total_r_community_cov.get(f"{gn}\t{step2}", "0\t0\t0").split('\t')[2])
                    avg_cov = (cov1 + cov2) / 2
                    tax_label = self.bin2cat.get(gn, ["Unknown"] * 7)[DiagramInputs.TAXONOMIC_LEVELS[self.user_tax]]
                    total_r_community_cov2[key] = f"{tax_label}\t{avg_cov:.6f}"

        os.makedirs(os.path.join(self.fig_input_dir), exist_ok=True)

        with open(os.path.join(self.fig_input_dir, 'Metabolic_Sankey_diagram_input.txt'), "w") as f:
            for line in total_r_community_cov.values():
                f.write(f"{line}\n")

        with open(os.path.join(self.fig_input_dir, 'Functional_network_input.txt'), "w") as f:
            f.write("#Genome\tStep1\tStep2\tTaxonomic Group\tCoverage value(average)\n")
            for k, v in total_r_community_cov2.items():
                f.write(f"{k}\t{v}\n")
        
        
        return self.bin2cat
    
    def _parse_gtdb_summary(self, file_path):
        if not os.path.exists(file_path):
            return
        with open(file_path) as f:
            next(f)  # skip header
            for line in f:
                fields = line.strip().split('\t')
                if len(fields) < 2:
                    continue
                user_genome, taxonomy = fields[0], fields[1]
                cat_fields = {
                    'phylum': re.search(r"p__([\w\-]+)", taxonomy),
                    'class': re.search(r"c__([\w\-]+)", taxonomy),
                    'order': re.search(r"o__([\w\-]+)", taxonomy),
                    'family': re.search(r"f__([\w\-]+)", taxonomy),
                    'genus': re.search(r"g__([\w\-]+)", taxonomy),
                    'species': re.search(r"s__([\w\-]+)", taxonomy)
                }
                cat_list = [
                    (cat_fields[k].group(1) if cat_fields[k] else f"NK_{k}")
                    for k in ['phylum', 'class', 'order', 'family', 'genus', 'species']
                ]
                self.bin2cat[user_genome] = cat_list + [user_genome]
    
    def calculate_mw_score(self):

        if self.user_tax not in DiagramInputs.TAXONOMIC_LEVELS:
            raise ValueError("Your input taxonomy is wrong. It should be one of: phylum, class, order, family, genus, species, or bin.")

        logging.info("Calculating MW-score...")

        result_dir = os.path.join(self.paths['output_dir'], 'MW-score_result')
        os.makedirs(result_dir, exist_ok=True)

        mw_functions = {}
        mw_function_hmm_ids = set()

        with open(os.path.join(self.paths['templates'], 'MW-score_reaction_table.txt')) as f:
            for line in f:
                if not line.startswith("#"):
                    parts = line.strip().split("\t")
                    mw_functions[parts[0]] = parts[2]
                    if ";" not in parts[2]:
                        for hmm in parts[2].split(","):
                            mw_function_hmm_ids.add(hmm)
                    else:
                        for group in parts[2].split(";"):
                            for hmm in group.split(","):
                                if "NO" not in hmm:
                                    mw_function_hmm_ids.add(hmm)

        mw_score_hash = {}
        for gn in sorted(self.hmmscan_result):
            for key in sorted(mw_functions):
                mw_score_hash.setdefault(key, {})[gn] = 0
                hmms = mw_functions[key]
                if ";" not in hmms:
                    for hmm_id in mw_function_hmm_ids:
                        if hmm_id in hmms and self.hmmscan_result[gn].get(hmm_id):
                            mw_score_hash[key][gn] = 1
                else:
                    part1, part2 = hmms.split(";")
                    logic1 = logic2 = 0
                    if "NO" not in part2:
                        for hmm_id in mw_function_hmm_ids:
                            if hmm_id in part1 and self.hmmscan_result[gn].get(hmm_id):
                                logic1 = 1
                            if hmm_id in part2 and self.hmmscan_result[gn].get(hmm_id):
                                logic2 = 1
                        if logic1 and logic2:
                            mw_score_hash[key][gn] = 1
                    else:
                        logic1, logic2 = 0, 1
                        for hmm_id in mw_function_hmm_ids:
                            if hmm_id in part1 and self.hmmscan_result[gn].get(hmm_id):
                                logic1 = 1
                            if hmm_id in part2 and self.hmmscan_result[gn].get(hmm_id):
                                logic2 = 0
                        if logic1 and logic2:
                            mw_score_hash[key][gn] = 1

        mw_community_cov = {}
        if self.input_csv_check:
            for pth in sorted(mw_score_hash):
                for gn in sorted(self.hmmscan_result):
                    if self.genome_coverage.get(gn) and mw_score_hash[pth].get(gn):
                        cat = self.bin2cat.get(gn, ["Unknown"] * 7)[DiagramInputs.TAXONOMIC_LEVELS[self.user_tax]]
                        key = f"{gn}\t{pth}"
                        mw_community_cov[key] = f"{cat}\t{pth}\t{self.genome_coverage[gn]:.6f}"

        # Extended KO hit output with coverage and taxonomy
        ko_hits_cov_output = os.path.join(result_dir, "MW-score_KO_hits_detailed.txt")
        with open(ko_hits_cov_output, "w") as f:
            f.write("Genome\tFunction\tKO\tCoverage\tTaxonomic_Group\n")
            for func, genomes in mw_score_hash.items():
                hmm_string = mw_functions[func]
                ko_candidates = set()
                if ";" in hmm_string:
                    for group in hmm_string.split(";"):
                        ko_candidates.update(group.split(","))
                else:
                    ko_candidates.update(hmm_string.split(","))

                for gn in genomes:
                    if not mw_score_hash[func][gn]:  # skip genomes that didn't pass logic gate
                        continue
                    coverage = self.genome_coverage.get(gn)
                    if not coverage:
                        continue
                    taxon = self.bin2cat.get(gn, ["Unknown"] * 7)[DiagramInputs.TAXONOMIC_LEVELS[self.user_tax]]
                    for ko in ko_candidates:
                        if self.hmmscan_result[gn].get(ko):
                            f.write(f"{gn}\t{func}\t{ko}\t{coverage:.6f}\t{taxon}\n")

        mw_community_cov2 = {}
        for gn in sorted(self.hmmscan_result):
            paths = [k.split("\t")[1] for k in mw_community_cov if k.startswith(f"{gn}\t")]
            for i in range(len(paths)):
                for j in range(i + 1, len(paths)):
                    step1, step2 = paths[i], paths[j]
                    key = f"{gn}\t{step1}\t{step2}"
                    cov1 = float(mw_community_cov.get(f"{gn}\t{step1}", "0\t0\t0").split("\t")[2])
                    cov2 = float(mw_community_cov.get(f"{gn}\t{step2}", "0\t0\t0").split("\t")[2])
                    avg_cov = (cov1 + cov2) / 2
                    cat = self.bin2cat.get(gn, ["Unknown"] * 7)[DiagramInputs.TAXONOMIC_LEVELS[self.user_tax]]
                    mw_community_cov2[key] = f"{cat}\t{avg_cov:.6f}"

        result_table_path = os.path.join(result_dir, "MW-score_result_table_input.txt")

        with open(result_table_path, "w") as f:
            f.write("#Genome\tFunc1\tFunc2\tTaxonomic Group\tCoverage value(average)\n")
            for k, v in sorted(mw_community_cov2.items()):
                f.write(f"{k}\t{v}\n")

        return result_table_path

    def finalize_mw_score(self):
        input_path = os.path.join(self.paths['output_dir'], 'MW-score_result', 'MW-score_result_table_input.txt')
        output_path = os.path.join(self.paths['output_dir'], 'MW-score_result', 'MW-score_result.txt')

        input_data = {}
        with open(input_path) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                fields = line.strip().split("\t")
                if len(fields) < 5:
                    continue
                input_data[line.strip()] = fields

        output1 = defaultdict(lambda: defaultdict(float))  # func -> cat -> summed coverage
        output2 = defaultdict(float)  # func -> summed coverage
        cat_set = set()

        for values in input_data.values():
            genome, func1, func2, cat, coverage = values
            coverage = float(coverage)
            output1[func1][cat] += coverage
            output1[func2][cat] += coverage
            output2[func1] += coverage
            output2[func2] += coverage
            cat_set.add(cat)

        output3 = {}  # func -> percent of total
        total_cov = sum(output2.values())
        for func in output2:
            output3[func] = round((output2[func] / total_cov) * 100, 1) if total_cov else 0.0

        output4 = defaultdict(dict)  # func -> cat -> percent
        for func in output1:
            for cat in cat_set:
                if output2[func] > 0:
                    output4[func][cat] = round((output1[func].get(cat, 0) / output2[func]) * 100, 1)
                else:
                    output4[func][cat] = 0.0

        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        with open(output_path, "w") as out_f:
            categories = sorted(cat_set)
            out_f.write("Function\tMW-score for each function\t" + "\t".join(categories) + "\n")
            for func in sorted(output4):
                scores = [str(output4[func].get(cat, 0)) for cat in categories]
                out_f.write(f"{func}\t{output3[func]}\t" + "\t".join(scores) + "\n")
