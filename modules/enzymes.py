import os
from .functions import run_parallel
import logging
from collections import defaultdict
import re

class EnzymePipeline:
    def __init__(self, paths, threads, eval_cutoff=1e-15, coverage_cutoff=0.35):
        self.paths = paths
        self.threads = threads
        self.eval_cutoff = eval_cutoff
        self.coverage_cutoff = coverage_cutoff
    
    def run(self):
        # Run dbCAN2 only
        return self.run_dbcan()
    
    def run_dbcan(self):
        dbcan_output_dir = os.path.join(self.paths['output_dir'], "intermediate_files", "dbCAN2_Files")
        os.makedirs(dbcan_output_dir, exist_ok=True)

        pending_genomes = []

        # Iterate only over individual .faa files (skip total.faa)
        for faa_file in self.individual_faa_files:
            gn_id = os.path.splitext(faa_file)[0]
            dm_path = os.path.join(dbcan_output_dir, f"{gn_id}.dbCAN2.out.dm")

            if not os.path.exists(dm_path) or os.path.getsize(dm_path) == 0:
                faa_path = os.path.join(self.paths['faa_dir'], faa_file)
                pending_genomes.append((gn_id, faa_path))

        if pending_genomes:
            script_path = os.path.join(self.paths['output_dir'], "tmp_run_dbCAN2.sh")
            dbcan_hmm_path = os.path.join(self.paths['databases'], 'dbCAN2', 'dbCAN-fam-HMMs.txt')

            if not os.path.exists(dbcan_hmm_path):
                raise FileNotFoundError(f"Missing dbCAN2 HMM file: {dbcan_hmm_path}")

            with open(script_path, 'w') as script_file:
                for gn_id, faa_path in pending_genomes:
                    dm_out = os.path.join(dbcan_output_dir, f"{gn_id}.dbCAN2.out.dm")
                    text_out = os.path.join(dbcan_output_dir, f"{gn_id}.dbCAN2.out")
                    script_file.write(
                        f"hmmscan --domtblout {dm_out} --cpu 1 {dbcan_hmm_path} {faa_path} > {text_out};\n"
                    )

            run_parallel(script_path, self.paths['output_dir'], self.threads)
            os.remove(script_path)
        else:
            logging.info("✅ All dbCAN2 results already present — skipping hmmscan step.")

        # parse .dm files and run post-processing
        dbcan_results = {}
        dbcan_hits = {}

        for faa_file in sorted(os.listdir(self.paths['faa_dir'])):
            if self.individual_faa_files:
                gn_id = os.path.splitext(faa_file)[0]
                dm_path = os.path.join(dbcan_output_dir, f"{gn_id}.dbCAN2.out.dm")

                if not os.path.exists(dm_path):
                    logging.error(f"Missing dbCAN2 output file for {gn_id}: {dm_path}")
                    continue

                parsed_hits = self.parse_and_filter_dbcan_dm(dm_path)
                for hmm_id, seq_id in parsed_hits:
                    dbcan_results.setdefault(gn_id, {}).setdefault(hmm_id, 0)
                    dbcan_results[gn_id][hmm_id] += 1

                    dbcan_hits.setdefault(gn_id, {}).setdefault(hmm_id, [])
                    dbcan_hits[gn_id][hmm_id].append(seq_id)

        return dbcan_results, dbcan_hits

    def parse_and_filter_dbcan_dm(self, dm_path):
        """
        Parse and filter dbCAN2 domtblout (.dm) results.
        - Keeps only GH and PL families.
        - Filters by E-value and coverage.
        - Removes overlapping regions, keeping best e-value match.

        Returns:
            List of (formatted HMM ID, query name)
        """
        raw_hits = defaultdict(list)

        with open(dm_path) as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue

                parts = line.strip().split()
                if len(parts) < 23:
                    logging.warning(f"Skipping line with insufficient fields: {line.strip()}")
                    continue

                hmm_id_full = parts[0]        # GH123.hmm
                hmm_len = int(parts[2])       # HMM length
                query_id = parts[3]           # protein name
                evalue = float(parts[6])
                ali_from = int(parts[17])
                ali_to = int(parts[18])

                hmm_prefix = hmm_id_full.split('.')[0]
                if not hmm_prefix.startswith(("GH", "PL")):
                    continue  # Skip anything not GH or PL

                if hmm_len <= 0 or ali_to <= ali_from:
                    continue

                coverage = (ali_to - ali_from) / hmm_len
                if evalue > self.eval_cutoff or coverage < self.coverage_cutoff:
                    continue

                raw_hits[query_id].append({
                    "hmm_id": hmm_id_full,
                    "evalue": evalue,
                    "ali_from": ali_from,
                    "ali_to": ali_to,
                })

        # Collapse overlapping hits by query and region
        filtered = []
        for query_id, hits in raw_hits.items():
            hits.sort(key=lambda h: (h["ali_from"], h["ali_to"]))
            kept_hits = []

            for hit in hits:
                overlap = False
                for kept in kept_hits:
                    overlap_len = min(hit["ali_to"], kept["ali_to"]) - max(hit["ali_from"], kept["ali_from"])
                    len1 = hit["ali_to"] - hit["ali_from"]
                    len2 = kept["ali_to"] - kept["ali_from"]
                    if (
                        overlap_len > 0 and
                        (overlap_len / len1 > 0.5 or overlap_len / len2 > 0.5)
                    ):
                        if hit["evalue"] < kept["evalue"]:
                            kept_hits.remove(kept)
                            kept_hits.append(hit)
                        overlap = True
                        break
                if not overlap:
                    kept_hits.append(hit)

            for hit in kept_hits:
                hmm_id_base = hit["hmm_id"].split(".")[0]  # GH1 → GH
                match = re.match(r"([A-Z]+)(\d+)", hmm_id_base)
                if not match:
                    continue
                fam_prefix, fam_number = match.groups()
                hmm_id = f"{fam_prefix}{int(fam_number):03d}"  # GH001
                filtered.append((hmm_id, query_id))

        return filtered
    
    def run_merops(self): 
        """
        Run MEROPS peptidase analysis using DIAMOND against the pepunit database.
        Returns:
            - merops_out: genome => merops_id => count
            - merops_out_hits: genome => merops_id => semicolon-separated names
            - merops_ids: set of all merops IDs found
        """
        merops_dir_out = os.path.join(self.paths['output_dir'], 'intermediate_files', 'MEROPS_Files')
        os.makedirs(merops_dir_out, exist_ok=True)

        pending_genomes = []
        for faa_file in self.individual_faa_files:
            gn_id = os.path.splitext(faa_file)[0]
            out_path = os.path.join(merops_dir_out, f"{gn_id}.MEROPSout.m8")
            if not os.path.exists(out_path) or os.path.getsize(out_path) == 0:
                pending_genomes.append((gn_id, faa_file))

        if not pending_genomes:
            logging.info("✅ Skipping MEROPS: all results already exist. Proceeding to parse existing outputs.")
        else:
            logging.info(f"🧬 Running MEROPS on {len(pending_genomes)} genomes.")
        # You would likely trigger your DIAMOND alignment logic here
            # Generate DIAMOND blast shell script
            db_path = os.path.join(self.paths['databases'],'MEROPS', 'pepunit.dmnd')
            script_path = os.path.join(self.paths['output_dir'], 'tmp_run_MEROPS.sh')
            with open(script_path, 'w') as script_file:
                for faa_file in os.listdir(self.paths['faa_dir']):
                    if self.individual_faa_files:
                        gn_id = os.path.splitext(faa_file)[0]
                        faa_path = os.path.join(self.paths['faa_dir'], faa_file)
                        out_path = os.path.join(merops_dir_out, f"{gn_id}.MEROPSout.m8")
                        db_path = os.path.join(self.paths['databases'],'MEROPS', 'pepunit.dmnd')
                        script_file.write(
                            f"diamond blastp -d {db_path} -q {faa_path} -o {out_path} "
                            f"-k 1 -e 1e-10 --query-cover 80 --id 50 --quiet -p 1 2> /dev/null\n"
                        )

            run_parallel(script_path, self.paths['output_dir'], self.threads)
            os.remove(script_path)
            logging.info("MEROPS DIAMOND analysis completed.")
            
        # Parse MEROPS pepunit.lib
        logging.info("Parsing MEROPS pepunit.lib file.")
        merops_map = {}  # MER ID => full header line
        with open(os.path.join(self.paths['databases'],'MEROPS', 'pepunit.lib'), "r", encoding="latin-1") as lib_file:

            for line in lib_file:
                if line.startswith(">"):
                    line = line.strip().replace('\r', '')
                    mer_id = line.split(" ")[0][1:]
                    merops_map[mer_id] = line
        logging.info(f"Parsed {len(merops_map)} MEROPS entries.")

        # Parse DIAMOND output
        logging.info("Parsing MEROPS DIAMOND output.")
        merops_out = defaultdict(lambda: defaultdict(int))         # genome -> merops_id -> count
        merops_out_hits = defaultdict(lambda: defaultdict(list))   # genome -> merops_id -> list of hits
        merops_ids = set()

        for m8_file in os.listdir(merops_dir_out):
            if m8_file.endswith(".MEROPSout.m8"):
                gn_id = m8_file.split(".")[0]
                with open(os.path.join(merops_dir_out, m8_file)) as infile:
                    for line in infile:
                        parts = line.strip().split("\t")
                        if len(parts) < 2:
                            continue
                        query_name = parts[0]
                        subject_id = parts[1]
                        lib_id = subject_id.split("|")[1] if "|" in subject_id else subject_id
                        lib_entry = merops_map.get(lib_id)

                        if not lib_entry or "#" not in lib_entry:
                            logging.warning(f"Missing or malformed entry for {lib_id} in pepunit.lib.")
                            continue
                        match = re.search(r"#(.+?)#", lib_entry)
                        if not match:
                            continue

                        merops_id = match.group(1)
                        merops_ids.add(merops_id)

                        merops_out[gn_id][merops_id] += 1
                        merops_out_hits[gn_id][merops_id].append(query_name)

        # Convert hits to semicolon-separated strings
        formatted_hits = {
            gn: {mid: ";".join(names) for mid, names in merops_out_hits[gn].items()}
            for gn in merops_out_hits
        }
        return merops_out, formatted_hits, merops_ids

    @property
    def individual_faa_files(self):
        return [
            f for f in os.listdir(self.paths['faa_dir'])
            if f.endswith('.faa') and f != 'total.faa'
        ]
    