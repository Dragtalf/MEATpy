import os
import re
import logging
from Bio import SeqIO
from .functions import run_parallel, get_check_score
from .hmm_utils import HMMThreshold

class HMMScanner:
    def __init__(self, paths, total_hmm2threshold, motifs, motif_pairs, motif_filter, threads):
        self.paths = paths
        self.total_hmm2threshold = total_hmm2threshold
        self.motifs = motifs
        self.motif_pairs = motif_pairs
        self.motif_filter = motif_filter
        self.threads = threads
        self.hmm_out = os.path.join(self.paths['output_dir'], 'intermediate_files', 'Hmmsearch_Outputs')
        self.hmm_db_address = os.path.join(self.paths['databases'], 'METABOLIC_hmm_db')

    def extract_motif_sequence(self, sequence_id, faa_dir, motifs, motif_name):
        """
        Extracts a specific sequence from total.faa and compiles a motif regex.
        Returns the sequence and compiled motif if found, otherwise (None, None).
        """
        fasta_file = os.path.join(faa_dir, "total.faa")
        sequence = None
        seq_data = ''
        recording = False

        try:
            with open(fasta_file, 'r') as file:
                for line in file:
                    line = line.strip()
                    if line.startswith('>'):
                        current_id = line[1:].split()[0]
                        if recording:
                            # finished reading the sequence we care about
                            break
                        recording = (current_id == sequence_id)
                        logging.debug(f"🔍 Looking for sequence ID: {sequence_id} current ID: {current_id} (recording: {recording})")
                    elif recording:
                        seq_data += line
            if recording:
                sequence = seq_data
        except FileNotFoundError:
            logging.error(f"FASTA file not found: {fasta_file}")
            return None, None

        if sequence:
            motif_regex = motifs.get(motif_name, '').replace('X', '[ARNDCQEGHILKMFPSTWYV]')
            return sequence, motif_regex

        logging.warning(f"❌ Sequence ID '{sequence_id}' not found or no sequence collected.")
        return None, None
        
    def prepare_hmm_search_script(self):
        script_path = os.path.join(self.paths['output_dir'], 'tmp_run_hmmsearch.sh')
        total_faa_path = os.path.join(self.paths['faa_dir'], 'total.faa')

        try:
            if not os.path.exists(total_faa_path):
                logging.error(f"File {total_faa_path} does not exist.")
                return None

            os.makedirs(self.hmm_out, exist_ok=True)
            
            commands_written = 0

            with open(script_path, 'w') as script_file:
                for hmm, threshold_obj in self.total_hmm2threshold.items():
                    if not isinstance(threshold_obj, HMMThreshold):
                        logging.error(f"Invalid threshold object for {hmm}: {threshold_obj}")
                        continue

                    threshold = threshold_obj.threshold
                    score_type = threshold_obj.score_type

                    result_path = os.path.join(self.hmm_out, f"{hmm}.hmmsearch_result.txt")
                    if os.path.exists(result_path) and os.path.getsize(result_path) > 1024:
                        continue  # Already processed

                    # Decide which database folder to look in
                    hmm_db = (
                        os.path.join(self.paths['databases'], 'kofam_database', 'profiles')
                        if hmm.startswith('K') and hmm[1:6].isdigit()
                        else os.path.join(self.paths['databases'], 'METABOLIC_hmm_db')
                    )
                    
                    hmm_path = os.path.join(hmm_db, hmm)

                    logging.info(f"hmm_path: {hmm_path} | hmm: {hmm} | threshold: {threshold} | score_type: {score_type}")

                    if not os.path.exists(hmm_path):
                        logging.warning(f"⚠️ Skipping {hmm}: HMM file not found at {hmm_path}")
                        continue

                    cmd = f"hmmsearch -{'T' if score_type == 'full' else '-domT'} {threshold} --cpu 1 --tblout {result_path} {hmm_path} {total_faa_path}\n"
                    
                    logging.debug(f"Writing command: {cmd}")
                    
                    script_file.write(cmd)
                    commands_written += 1
                    
            if commands_written == 0:
                logging.info("✅ All HMM searches are complete. No commands written.")
                return None
            else:
                logging.info(f"📝 Prepared {commands_written} HMM search commands.")

        except Exception as e:
            logging.error(f"An error occurred: {e}")
        
        logging.info(f"Running HMM search script: {script_path}")    
        run_parallel(script_path, self.paths['output_dir'], self.threads)
        logging.info("HMM search script completed.")
        
        try:
            os.remove(script_path)
        except Exception as e:
            logging.error(f"Failed to clean up HMM search script: {str(e)}")
        
    def prepare_hmm_search_data(self):
        """
        Map sequence IDs to genome IDs and aggregate FAA and gene sequences.
        Returns:
        - genome_ids: set of genome IDs
        - total_faa_sequences: dict of all FAA sequences {seq_id: sequence}
        - total_gene_sequences: dict of all gene sequences {seq_id: sequence}
        - seqid_to_genomeid: dict mapping seq_id -> genome_id
        """
        genome_ids = set()
        seqid_to_genomeid = {}
        total_faa_seq = {}
        total_gene_seq = {}

        faa_files = [f for f in os.listdir(self.paths['faa_dir']) if f.endswith('.faa') and f != 'total.faa']

        for faa_file in faa_files:
            file_path = os.path.join(self.paths['faa_dir'], faa_file)
            genome_id = os.path.splitext(faa_file)[0]
            genome_ids.add(genome_id)

            for record in SeqIO.parse(file_path, "fasta"):
                total_faa_seq[record.id] = str(record.seq)
                seqid_to_genomeid[record.id] = genome_id

            if self.paths['gene_dir']:
                gene_path = os.path.join(self.paths['gene_dir'], f"{genome_id}.gene")
                if os.path.exists(gene_path):
                    for record in SeqIO.parse(gene_path, "fasta"):
                        total_gene_seq[record.id] = str(record.seq)

        return genome_ids, total_faa_seq, total_gene_seq, seqid_to_genomeid

    def summarize_hmmsearch_results(self, seqid_to_genomeid):
        Hmmscan_result = {}  # genome_name => hmm => count
        Hmmscan_hits = {}    # genome_name => hmm => hit details (e.g. sequences)
        Hmm_id = {}          # Unique list of HMMs used

        # Find all hmmsearch result files
        hmmsearch_files = [
            os.path.join(dp, f)
            for dp, _, filenames in os.walk(self.hmm_out)
            for f in filenames if f.endswith('.hmmsearch_result.txt')
        ]

        for file_name in hmmsearch_files:
            base_name = os.path.basename(file_name)
            match = re.match(r"(.+?)\.hmmsearch_result\.txt", base_name)
            if not match:
                logging.warning(f"⚠️ Could not parse HMM name from {base_name}")
                continue

            hmm = match.group(1)
            Hmm_id[hmm] = 1

            with open(file_name, 'r') as file:
                for line in file:
                    if line.startswith("#"):
                        continue

                    line = re.sub(r"\s+", "\t", line.strip())
                    tmp = line.split('\t')

                    seq_id = tmp[0]
                    gn_id = seqid_to_genomeid.get(seq_id)
                    
                    if not gn_id:
                        logging.warning(f"❌ Genome ID not found for sequence ID {seq_id}. Skipping...")
                        continue

                    # Get threshold and score type (full vs domain)
                    try:
                        threshold_obj = self.total_hmm2threshold.get(hmm)
                        if not threshold_obj:
                            logging.error(f"No threshold found for {hmm}")
                            continue
                        threshold = threshold_obj.threshold
                        score_type = threshold_obj.score_type
                    except Exception as e:
                        logging.error(f"❌ Failed to parse threshold for {hmm}: {e}")
                        continue

                    score = float(tmp[8]) if score_type == "domain" else float(tmp[5])

                    # Pass threshold → record hit
                    if score >= threshold:
                        self.process_hit(gn_id, hmm, tmp, Hmmscan_hits, Hmmscan_result)
                        logging.debug(f"✔️ Hit passed threshold ({score} >= {threshold}) for {hmm} in {gn_id}")
                    else:
                        logging.debug(f"❌ Score {score} below threshold {threshold} for {hmm} in {gn_id}")
        
        return Hmmscan_result, Hmmscan_hits, list(Hmm_id.keys())

    def process_hit(self, gn_id, hmm, tmp, Hmmscan_hits, Hmmscan_result):
        hmm_basename = hmm.split('.')[0]
        sequence_id = tmp[0].lstrip('>')
        logging.debug(f"Processing hit for {hmm_basename} in {gn_id} with sequence ID {sequence_id}")
        
        if hmm_basename in self.motifs:
            
            seq, motif_new = self.extract_motif_sequence(sequence_id, self.paths['faa_dir'], self.motifs, hmm_basename)
            if not seq:
                logging.warning(f"❌ No sequence found for {sequence_id}")
            elif not re.search(motif_new, seq):
                logging.debug(f"❌ Motif not matched for gn_id: {gn_id} {hmm_basename}: {motif_new} vs {seq[:50]}...")
            else:
                logging.debug(f"✅ Motif matched for {hmm_basename} in {gn_id}")
                self.update_results(gn_id, hmm, sequence_id, Hmmscan_hits, Hmmscan_result)
        elif hmm_basename in self.motif_pairs:
            logging.debug(f"🔍 Checking motif pairs for {hmm_basename} gn_id: {gn_id}")
            self.motif_filter.check_motif_pairs(gn_id, hmm, sequence_id, Hmmscan_hits, Hmmscan_result, self.motif_pairs, self.threads)
        else:
            self.update_results(gn_id, hmm, sequence_id, Hmmscan_hits, Hmmscan_result)

    def update_results(self, gn_id, hmm, sequence_id, Hmmscan_hits, Hmmscan_result):
        hmm = hmm.strip()
        if gn_id not in Hmmscan_hits:
            Hmmscan_hits[gn_id] = {}
        if hmm not in Hmmscan_hits[gn_id]:
            Hmmscan_hits[gn_id][hmm] = sequence_id
        else:
            Hmmscan_hits[gn_id][hmm] += "," + sequence_id

        if gn_id not in Hmmscan_result:
            Hmmscan_result[gn_id] = {}
        if hmm not in Hmmscan_result[gn_id]:
            Hmmscan_result[gn_id][hmm] = 0
        Hmmscan_result[gn_id][hmm] += 1
        
        logging.debug(f"✅ Updated HMM result: {gn_id} | {hmm} => Count: {Hmmscan_result[gn_id][hmm]} | Hits: {Hmmscan_hits[gn_id][hmm]}")

class Motif_Filter:
    def __init__(self, paths, threads):
        self.paths = paths
        self.threads = threads
        self.hmm_db_address = os.path.join(self.paths['databases'], 'METABOLIC_hmm_db')

    def check_motif_pairs(self, gn_id, hmm, sequence_id, Hmmscan_hits, Hmmscan_result, motif_pairs, threads):
        try:
            hmm_basename = hmm.split('.')[0]
    
            motif_hmm = os.path.join(self.hmm_db_address, f"{hmm_basename}.check.hmm")
            motif_anti_hmm = os.path.join(self.hmm_db_address, f"{motif_pairs[hmm_basename]}.check.hmm")

            input_faa = os.path.join(self.paths['faa_dir'], "total.faa")
            
            fasta_id = f">{sequence_id}"  # For matching FASTA lines
            
            input_seq_file = os.path.join(self.paths['output_dir'], f"{hmm_basename}.check.faa")
            result_file = os.path.join(self.paths['output_dir'], f"{hmm_basename}.check.hmmsearch_result.txt")
            result_file_anti = os.path.join(self.paths['output_dir'], f"{motif_pairs[hmm_basename]}.check.hmmsearch_result.txt")
            script_path = os.path.join(self.paths['output_dir'], "run_hmmsearch.sh")

            # Extract the sequence to a file
            with open(input_faa, 'r') as faa, open(input_seq_file, 'w') as out_faa:
                record = False
                for line in faa:
                    if line.startswith('>') and record:
                        break
                    if line.startswith(fasta_id):
                        record = True
                    if record:
                        out_faa.write(line)

            # Write commands to a script file
            with open(script_path, 'w') as script_file:
                script_file.write(f"hmmsearch --cpu 1 --tblout {result_file} {motif_hmm} {input_seq_file}\n")
                script_file.write(f"hmmsearch --cpu 1 --tblout {result_file_anti} {motif_anti_hmm} {input_seq_file}\n")

            # Run HMM searches in parallel
            run_parallel(script_path, self.paths['output_dir'], self.threads)

            # After running, calculate scores
            motif_check_score = get_check_score(result_file)
            motif_anti_check_score = get_check_score(result_file_anti)

            logging.debug(f"Motif check score for {hmm_basename}: {motif_check_score}")
            logging.debug(f"Motif anti-check score for {hmm_basename}: {motif_anti_check_score}")

            # Update hits and results based on scores
            if motif_check_score >= motif_anti_check_score and motif_check_score != 0:
                logging.debug(f"✅ Motif check > anti check score passed for {hmm_basename} in {gn_id}")
                # Ensure dictionaries exist
                if gn_id not in Hmmscan_hits:
                    Hmmscan_hits[gn_id] = {}
                if gn_id not in Hmmscan_result:
                    Hmmscan_result[gn_id] = {}

                # Now safely update
                if hmm not in Hmmscan_hits[gn_id]:
                    Hmmscan_hits[gn_id][hmm] = sequence_id
                    Hmmscan_result[gn_id][hmm] = 1
                    logging.debug(f"✅ New hit recorded for {hmm_basename} in {gn_id}: {sequence_id}")
                else:
                    Hmmscan_hits[gn_id][hmm] += "," + sequence_id
                    Hmmscan_result[gn_id][hmm] += 1
                    logging.debug(f"✅ Updated hit for {hmm_basename} in {gn_id}: {Hmmscan_hits[gn_id][hmm]}: {sequence_id}")
            else:logging.debug(f"❌ Motif check score failed for {hmm_basename} in {gn_id}")

        finally:
            # Cleanup
            for file in (input_seq_file, result_file, result_file_anti):
                try:
                    os.remove(file)
                except Exception as e:
                    logging.warning(f"Could not remove temporary file {file}: {e}")

    def cleanup_temp_files(self, basename):
        try:
            os.remove(os.path.join(self.paths['output_dir'], f"tmp.{basename}.check.faa"))
            os.remove(os.path.join(self.paths['output_dir'], f"{basename}.check.hmmsearch_result.txt"))
        except Exception as e:
            logging.error(f"Failed to clean up temporary files: {str(e)}")

class HMMCollector:
    def __init__(self, paths):
        self.paths = paths
    
    def generate_hmm_collections(self, hmm_ids, hmmscan_hits, total_faa_seq, total_gene_seq):

        output_dir = os.path.join(self.paths['output_dir'], "Each_HMM_Amino_Acid_Sequence")
        os.makedirs(output_dir, exist_ok=True)

        for hmm in sorted(hmm_ids):
            hmm_faa_seq = {}
            hmm_gene_seq = {}
            
            for gn_id, hits in hmmscan_hits.items():
                if hmm in hits:
                    for hit in hits[hmm].split(','):
                        hit = hit.strip().lstrip('>')  # strip > from beginning
                        if hit in total_faa_seq:
                            hmm_faa_seq[f">{gn_id}~~{hit}"] = total_faa_seq[hit]
                        else:
                            logging.warning(f"❌ {hit} not found in total_faa_seq")

                        if hit in total_gene_seq:
                            hmm_gene_seq[f">{gn_id}~~{hit}"] = total_gene_seq[hit]

            if hmm_faa_seq:
                with open(os.path.join(output_dir, f"{hmm}.collection.faa"), 'w') as f:
                    for key, seq in sorted(hmm_faa_seq.items()):
                        f.write(f"{key}\n{seq}\n")
            
            if hmm_gene_seq:
                with open(os.path.join(output_dir, f"{hmm}.collection.gene"), 'w') as f:
                    for key, seq in sorted(hmm_gene_seq.items()):
                        f.write(f"{key}\n{seq}\n")
