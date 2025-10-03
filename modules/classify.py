import os
import re
import csv
import subprocess
import shutil
import urllib.request
import gzip
import logging
from .functions import show_download_progress

class Classifier:
    def __init__(self, paths, threads, correction, revert_names):
        self.paths = paths
        self.threads = threads
        self.correction = correction
        self.revert_names = revert_names
        self.gtdbtk_bac = None
        self.gtdbtk_ar = None
        self.ncbi_bac = None
        self.ncbi_ar = None

    def classify_genomes(self):
        gtdbtk_dir = os.path.join(self.paths['output_dir'], 'intermediate_files', 'gtdbtk_Genome_files')
        scratch_dir = os.path.join(self.paths['output_dir'], 'intermediate_files', 'gtdbtk_Genome_files', 'gtdbtk_scratch')
        gtdbtk_bac= os.path.join(gtdbtk_dir, "gtdbtk.bac120.summary.tsv")
        ncbi_bac= os.path.join(gtdbtk_dir, "ncbi.bac120.correction.tsv")
        gtdbtk_ar= os.path.join(gtdbtk_dir, "gtdbtk.ar53.summary.tsv")
        ncbi_ar= os.path.join(gtdbtk_dir, "ncbi.ar53.correction.tsv")
        
        self.gtdbtk_bac = gtdbtk_bac
        self.gtdbtk_ar = gtdbtk_ar
        self.ncbi_bac = ncbi_bac
        self.ncbi_ar = ncbi_ar
        
        # Skip GTDB-Tk if any known GTDB or corrected NCBI summary files exist
        gtdbtk_bac_summary = os.path.join(gtdbtk_dir, "gtdbtk.bac120.summary.tsv")
        ncbi_bac_correction = os.path.join(gtdbtk_dir, "ncbi.bac120.correction.tsv")

        if os.path.exists(gtdbtk_bac_summary) or os.path.exists(ncbi_bac_correction):
            logging.info("GTDB-Tk summary already exists. Skipping classification.")
        else:
            logging.info("Running GTDB-Tk classification...")
            os.makedirs(gtdbtk_dir, exist_ok=True)

            subprocess.run([
                "gtdbtk", "classify_wf",
                "--cpus", str(self.threads),
                "-x", "fasta",
                "--genome_dir", self.paths['fasta_dir'],
                "--skip_ani_screen",
                "--out_dir", gtdbtk_dir,
                "--pplacer_cpus=1",
                "--scratch_dir", scratch_dir
            ], check=True)

            logging.info("GTDB-Tk classification completed.")

        if self.correction == "NCBI" and not os.path.exists(ncbi_bac):
            logging.info("Correcting to NCBI taxonomy...")
            bac_metadata = os.path.join(gtdbtk_dir, "bac120_metadata.tsv")
            ar_metadata = os.path.join(gtdbtk_dir, "ar53_metadata.tsv")

            if not os.path.exists(bac_metadata):
                self.download_and_extract("https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/bac120_metadata.tsv.gz", bac_metadata)

            if not os.path.exists(ar_metadata):
                self.download_and_extract("https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/ar53_metadata.tsv.gz", ar_metadata)

            logging.info("Running GTDB to NCBI correction...")
            
            subprocess.run([
                "gtdb_to_ncbi_majority_vote.py",
                "--gtdbtk_output_dir", gtdbtk_dir,
                "--output_file", os.path.join(gtdbtk_dir, "ncbi.bac120.correction.tsv"),
                "--bac120_metadata_file", bac_metadata
            ], check=True)

            if os.path.exists(os.path.join(gtdbtk_dir, "gtdbtk.ar53.summary.tsv")):
                subprocess.run([
                    "gtdb_to_ncbi_majority_vote.py",
                    "--gtdbtk_output_dir", gtdbtk_dir,
                    "--output_file", os.path.join(gtdbtk_dir, "ncbi.ar53.correction.tsv"),
                    "--bac120_metadata_file", ar_metadata
                ], check=True)

            logging.info("NCBI correction completed. Merging with GTDB-Tk output...")
            
            if os.path.exists(os.path.join(gtdbtk_dir, "ncbi.bac120.correction.tsv")):
                self.merge_ncbi_with_gtdb(
                    gtdbtk_bac,
                    os.path.join(gtdbtk_dir, "ncbi.bac120.correction.tsv"),
                    ncbi_bac
            )

            if os.path.exists(os.path.join(gtdbtk_dir, "ncbi.ar53.correction.tsv")):
                self.merge_ncbi_with_gtdb(
                    gtdbtk_ar,
                    os.path.join(gtdbtk_dir, "ncbi.ar53.correction.tsv"),
                    ncbi_ar
                )
            logging.info("GTDB-Tk classification completed with NCBI correction.")
            
        else:
            logging.info("GTDB-Tk classification completed.")   
    
    def parse_gtdb_summary(self):

        if self.correction == "NCBI":
            files = [self.ncbi_bac, self.ncbi_ar]
        else:
            files = [self.gtdbtk_bac, self.gtdbtk_ar]

        summary_files = [f for f in files if os.path.exists(f)]
        if not summary_files:
            raise FileNotFoundError("No GTDB summary file found.")

        bin2cat = {}
        for file in summary_files:
            with open(file) as f:
                next(f)
                for line in f:

                    fields = line.strip().split('\t')
                    if len(fields) < 2:
                        continue
                    user_genome, taxonomy = fields[0], fields[1]
                    genus_match = re.search(r"g__([\w\-]+)", taxonomy)
                    species_match = re.search(r"s__([\w\-]+)", taxonomy)
                    genus = genus_match.group(1) if genus_match else "UnknownGenus"
                    species = species_match.group(1) if species_match else "UnknownSpecies"
                    bin2cat[user_genome] = f"{genus}_{species}"
        return bin2cat

    def patch_summary_file(self, summary_path, rename_log):
        if not os.path.exists(summary_path):
            return
        backup_path = summary_path + ".bak"
        shutil.copy2(summary_path, backup_path)
        logging.info(f"Backed up original summary to: {backup_path}")

        rename_dict = {
            orig.replace('.fasta', ''): new.replace('.fasta', '')
            for orig, new in rename_log
        }

        with open(summary_path, 'r') as infile:
            header = infile.readline()
            lines = infile.readlines()

        with open(summary_path, 'w') as outfile:
            outfile.write(header)
            for line in lines:
                parts = line.strip().split('\t')
                if parts[0] in rename_dict:
                    old = parts[0]
                    parts[0] = rename_dict[old]
                    logging.debug(f"Renamed {old} → {parts[0]} in summary")
                outfile.write('\t'.join(parts) + '\n')

        logging.info(f"✅ Patched summary file: {summary_path}")

    def rename_inputs(self, output_ext=".fasta"):
        import csv

        log_path = os.path.join(self.paths['output_dir'], "MAG_rename_log.csv")
        input_folder = self.paths['fasta_dir']
        output_folder = self.paths['fasta_dir']

        if self.revert_names:
            # === REVERT LOGIC ===
            if not os.path.exists(log_path):
                logging.error(f"Log file not found: {log_path}. Cannot revert names.")
                exit(1)
            elif os.path.getsize(log_path) < 1024:
                logging.error(f"Log file is too small: {log_path}. Possibly incomplete. Aborting.")
                exit(1)

            with open(log_path, 'r') as csvfile:
                reader = csv.DictReader(csvfile)
                rename_log = [(row['original_filename'], row['new_filename']) for row in reader]

            missing = []
            for original_name, new_name in rename_log:
                src_path = os.path.join(input_folder, new_name)
                dst_path = os.path.join(output_folder, original_name)

                if os.path.exists(src_path):
                    shutil.move(src_path, dst_path)
                    logging.info(f"Reverted {new_name} → {original_name}")
                else:
                    logging.warning(f"❌ File not found for reverting: {src_path}")
                    missing.append((new_name, original_name))

            if missing:
                logging.error(f"{len(missing)} file(s) failed to revert. Log file NOT deleted.")
                return rename_log, [new for _, new in rename_log]

            logging.info(f"✅ All names reverted.")
            return rename_log, [new for _, new in rename_log]

        # === FORWARD RENAME LOGIC ===
        if os.path.exists(log_path):
            logging.error(f"Rename Log file already exists: {log_path}. Skipping renaming. Run with --revert_names to revert.")
            logging.error(f"Check {log_path} and your input folder.")
            exit(1)

        bin2cat = self.parse_gtdb_summary()

        assembler_map = {'megahit': 'MH', 'spades': 'SP'}
        binner_map = {'concoct': 'CON', 'maxbin2': 'MB2', 'metabat2': 'MBT'}
        rename_log = []

        os.makedirs(os.path.join(self.paths['input_dir'], 'fasta_backup'), exist_ok=True)

        for root, _, files in os.walk(input_folder):
            for filename in files:
                if not (filename.endswith(".fa") or filename.endswith(".fasta")):
                    logging.warning(f"Skipping non-fasta file: {filename}")
                    continue

                src_path = os.path.join(root, filename)
                base = os.path.basename(filename)

                # Don't strip Refined- until after match
                match = re.match(
                    r"([^-]+)-([^-]+)Refined-(.+?)\.(\d+(?:_sub)?)\.(fa|fasta)$",
                    base
                )
                if not match:
                    logging.warning(f"Unrecognized MAG filename format: {filename}")
                    continue

                assembler, binner, sample_id, mag_num, _ = match.groups()

                assembler_abbr = assembler_map.get(assembler.lower(), assembler[:2].upper())
                binner_abbr = binner_map.get(binner.lower(), binner[:3].upper())
                base_name = f"{sample_id}_{assembler_abbr}-{binner_abbr}_{mag_num}"
                name_no_ext = os.path.splitext(base)[0]
                taxonomy = bin2cat.get(name_no_ext, "UnknownGenus_UnknownSpecies")
                new_name = f"{taxonomy}__{base_name}{output_ext}"
                dst_path = os.path.join(output_folder, new_name)

                shutil.copy2(src_path, os.path.join(self.paths['input_dir'], 'fasta_backup', base))
                shutil.move(src_path, dst_path)

                rename_log.append((filename, new_name))
                logging.info(f"Renamed {filename} → {new_name}")

        # Write log
        with open(log_path, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(['original_filename', 'new_filename'])
            writer.writerows(rename_log)

        self.patch_summary_file(self.gtdbtk_bac, rename_log)
        self.patch_summary_file(self.gtdbtk_ar, rename_log)
        self.patch_summary_file(self.ncbi_bac, rename_log)
        self.patch_summary_file(self.ncbi_ar, rename_log)

        logging.info(f"Processed {len(rename_log)} files. Log saved to: {log_path}")
        return rename_log, files

    def merge_ncbi_with_gtdb(self, summary_file, correction_file, output_file):
        # This is from GTDB-Tk
        corrections = {}

        # Load the correction file
        with open(correction_file) as f:
            next(f)
            for line in f:
                genome_id, gtdb_class, ncbi_class = line.strip().split('\t')[:3]

                # Fix NCBI spacing inconsistencies
                ncbi_class = ncbi_class.replace(' ', '_')

                # Fallback strategy: if NCBI is just 's__', use GTDB
                if ncbi_class == "s__":
                    ncbi_class = gtdb_class
                else:
                    # If genus or species missing in NCBI, add from GTDB
                    ncbi_has_genus = bool(re.search(r'g__\S+', ncbi_class))
                    ncbi_has_species = bool(re.search(r's__\S+', ncbi_class))
                    gtdb_genus = re.search(r'g__\S+', gtdb_class)
                    gtdb_species = re.search(r's__\S+', gtdb_class)

                    if not ncbi_has_genus and gtdb_genus:
                        ncbi_class += ";" + gtdb_genus.group()
                    if not ncbi_has_species and gtdb_species:
                        ncbi_class += ";" + gtdb_species.group()

                corrections[genome_id] = ncbi_class

        # Merge into the summary file
        with open(summary_file) as in_f, open(output_file, 'w') as out_f:
            header = in_f.readline()
            out_f.write(header)
            for line in in_f:
                parts = line.strip().split('\t')
                genome_id = parts[0]
                if genome_id in corrections:
                    parts[1] = corrections[genome_id]
                out_f.write("\t".join(parts) + "\n")

    @staticmethod       
    def download_and_extract(url, target_path):
        gz_path = target_path + ".gz"
        urllib.request.urlretrieve(url, gz_path, reporthook=show_download_progress)
        with gzip.open(gz_path, 'rb') as f_in, open(target_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
        os.remove(gz_path)

