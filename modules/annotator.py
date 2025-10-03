import os
import logging
import csv
from Bio import SeqIO
from BCBio import GFF
from Bio.SeqRecord import SeqRecord
from .functions import run_parallel

class Annotator:
    def __init__(self, paths, p_method, threads):
        self.paths = paths
        self.prodigal_method = p_method
        self.threads = threads
        self.faa_dir = self.paths['faa_dir']
        self.gff_dir = self.paths['gff_dir']
        self.gene_dir = self.paths['gene_dir']

    def run_prodigal_annotation(self):
        """
        Annotate genomes using Prodigal and convert GFF output to .gene FASTA files.
        """
        
        annotation_script_path = os.path.join(self.paths['output_dir'], "run_annotate.sh")
        total_faa_file = os.path.join(self.faa_dir, "total.faa")

        # Clean existing total.faa
        if os.path.exists(total_faa_file):
            os.remove(total_faa_file)
            logging.info("Removed existing total.faa file.")

        rename_log_path = os.path.join(self.paths['output_dir'], "MAG_rename_log.csv")

        # If rename log exists, validate consistency
        if os.path.exists(rename_log_path) and os.path.exists(self.faa_dir):
            with open(rename_log_path) as f:
                reader = csv.DictReader(f)
                renamed = {row['new_filename'].replace('.fasta', '') for row in reader}

            # Get base names of existing faa/gff files
            existing_annotated = {
                os.path.splitext(f)[0]
                for f in os.listdir(self.faa_dir)
                if f.endswith(".faa")
            }

            # Check for mismatch
            missing = renamed - existing_annotated
            if missing:
                logging.error("Mismatch detected: renamed FASTA files exist, but some corresponding annotations are missing.")
                logging.error(f"Missing annotations for: {', '.join(sorted(missing))}")
                logging.error("If you've just renamed files, run with --oops to try and fix mismatches.")
                exit(1)

        # Discover input genomes
        fasta_files = sorted(f for f in os.listdir(self.paths['fasta_dir']) if f.endswith('.fasta'))
        if not fasta_files:
            logging.error("No FASTA files found in the input directory.")
            return
        
        os.makedirs(self.faa_dir, exist_ok=True)
        os.makedirs(self.gff_dir, exist_ok=True)

        commands = []
        for fasta_file in fasta_files:
            genome_id = os.path.splitext(fasta_file)[0]
            input_path = os.path.join(self.paths['fasta_dir'], fasta_file)
            faa_path = os.path.join(self.faa_dir, f"{genome_id}.faa")
            gff_path = os.path.join(self.gff_dir, f"{genome_id}.gff")

            if os.path.exists(faa_path) and os.path.exists(gff_path):
                continue

            cmd = f"prodigal -i {input_path} -a {faa_path} -o {gff_path} -f gff -p {self.prodigal_method} -q"
            commands.append(cmd)
            logging.debug(f"Scheduled Prodigal for {genome_id}")

        # Write and run script
        if commands:
            with open(annotation_script_path, 'w') as f:
                f.write("\n".join(commands) + "\n")
            logging.info("Running annotation script in parallel...")
            run_parallel(annotation_script_path, self.paths['output_dir'], self.threads)
        else:
            logging.info("All genomes annotated, skipping Prodigal run.")

        # Convert GFF to .gene FASTA
        logging.info("Converting GFF to .gene if missing...")
        os.makedirs(self.gene_dir, exist_ok=True)
        
        for fasta_file in fasta_files:
            genome_id = os.path.splitext(fasta_file)[0]
            fasta_path = os.path.join(self.paths['fasta_dir'], fasta_file)
            gff_path = os.path.join(self.gff_dir, f"{genome_id}.gff")
            gene_path = os.path.join(self.gene_dir, f"{genome_id}.gene")

            if not os.path.exists(gff_path) or os.path.getsize(gff_path) == 0:
                logging.warning(f"Missing or empty GFF for {genome_id}")
                continue
            elif os.path.exists(gene_path) and os.path.getsize(gene_path) > 0:
                continue
            
            # Load FASTA and preserve contig IDs
            fasta_dict = SeqIO.to_dict(SeqIO.parse(fasta_path, "fasta"))

            with open(gff_path) as gff_handle, open(gene_path, "w") as out_fasta:
                try:
                    for rec in GFF.parse(gff_handle, base_dict=fasta_dict):
                        contig_id = rec.id
                        feature_counter = 1

                        for feature in rec.features:
                            if 'ID' not in feature.qualifiers:
                                logging.debug(f"Skipping feature without ID in {genome_id}")
                                continue

                            # Give it a simple numbered suffix like _1, _2, etc.
                            feature_id = f"{contig_id}_{feature_counter}"
                            feature_counter += 1

                            feature_seq = feature.extract(rec.seq)
                            seq_record = SeqRecord(feature_seq, id=feature_id, description="")
                            SeqIO.write(seq_record, out_fasta, "fasta")
                except Exception as e:
                    logging.warning(f"⚠️ Failed to parse GFF for {genome_id}: {e}")

        if os.path.exists(annotation_script_path):
            os.remove(annotation_script_path)
            logging.info("Cleaned up annotation script.")

    def concatenate_faa_files(self):
        """
        Use Biopython to concatenate all .faa files, normalizing headers.
        Only includes `record.id` (no descriptions).
        """
        
        logging.info("Concatenating .faa files...")
        
        output_path = os.path.join(self.faa_dir, 'total.faa')
        faa_files = [
            os.path.join(self.faa_dir, f)
            for f in sorted(os.listdir(self.faa_dir))
            if f.endswith('.faa')
        ]

        if not faa_files:
            logging.warning("No .faa files found for concatenation.")
            return output_path

        if os.path.exists(output_path):
            output_mtime = os.path.getmtime(output_path)
            if all(os.path.getmtime(f) <= output_mtime for f in faa_files):
                return output_path

        all_records = []
        for filepath in faa_files:
            try:
                records = SeqIO.parse(filepath, "fasta")
                for rec in records:
                    # Normalize the header to be only the clean ID
                    clean_rec = SeqRecord(rec.seq, id=rec.id, description="")
                    all_records.append(clean_rec)
            except Exception as e:
                logging.warning(f"❌ Failed to parse {filepath}: {e}")

        SeqIO.write(all_records, output_path, "fasta")

        logging.info(f"✅ Concatenated {len(all_records)} records into {output_path}")

        return output_path
