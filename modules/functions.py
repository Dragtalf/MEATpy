import os
import sys
import subprocess
import datetime
import glob
import re
import logging
import math
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm

def get_hmm_to_ko_hash(Hmm_table_temp):
    """
    Convert HMM table (line number → line string) into a dict mapping:
    HMM ID → KO ID.hmm (e.g. TIGR02694 → K08355.hmm)
    """
    Hmm2ko = {}

    for line_no, temp in Hmm_table_temp.items():
        parts = temp.split('\t')
        hmm = parts[5]
        ko = parts[6]

        if hmm and ';' not in hmm:
            Hmm2ko[hmm] = f"{ko}.hmm"
        elif hmm:
            hmms = hmm.split('; ')
            kos = ko.split('; ')
            for h, k in zip(hmms, kos):
                Hmm2ko[h] = f"{k}.hmm"

    # Only include KO-like keys (start with 'K' followed by 5 digits)
    return {hmm: ko for hmm, ko in Hmm2ko.items() if re.match(r'^K\d{5}\.hmm$', ko)}

def get_motif(motif_file):
    """
    Load motifs from file. Replace 'X' with amino acid wildcard regex.
    Format: protein_id:motif_sequence
    Example line: dsrC:GPXKXXCXXXGXPXPXXCX
    """
    motifs = {}
    try:
        with open(motif_file, 'r') as file:
            for line in file:
                if ':' not in line:
                    logging.warning(f"Malformed motif line: {line.strip()}")
                    continue
                protein_id, motif_seq = line.strip().split(':', 1)
                # Replace 'X' with any amino acid (IUPAC one-letter codes)
                motifs[protein_id.strip()] = motif_seq.strip().replace('X', '[ARNDCQEGHILKMFPSTWYV]')
    except FileNotFoundError:
        logging.error(f"Motif file not found: {motif_file}")
    except Exception as e:
        logging.error(f"Error reading motif file: {e}")
    return motifs

def get_motif_pair(motif_pair_file):
    """
    Load motif pairs from file.
    Format: primary_protein_id:secondary_protein_id
    Example line: dsrC:tusE
    """
    motif_pairs = {}
    try:
        with open(motif_pair_file, 'r') as file:
            for line in file:
                if ':' not in line:
                    logging.warning(f"Malformed motif pair line: {line.strip()}")
                    continue
                primary, secondary = line.strip().split(':', 1)
                motif_pairs[primary.strip()] = secondary.strip()
    except FileNotFoundError:
        logging.error(f"Motif pair file not found: {motif_pair_file}")
    except Exception as e:
        logging.error(f"Error reading motif pair file: {e}")
    return motif_pairs

def get_genome_coverage_for_long_reads(reads, folder, output, test, omic_reads_type, sequencing_type, cpu_numbers):
    # Concatenate all the genes
    seq = {}
    for gene_file in glob.glob(f"{folder}/*.gene"):
        seq_name = os.path.basename(gene_file).split('.gene')[0]
        with open(gene_file, 'r') as file:
            for line in file:
                line = line.strip()
                if line.startswith('>'):
                    header = line.split()[0] if ' ' in line else line
                    seq[f">{seq_name}~~{header[1:]}"] = ""
                else:
                    seq[f">{seq_name}~~{header[1:]}"] += line

    # Write all concatenated gene sequences to a file
    all_genes_path = os.path.join(output, 'All_gene_collections.gene')
    with open(all_genes_path, 'w') as file:
        for header, sequence in seq.items():
            file.write(f"{header}\n{sequence}\n")

    # Prepare the alignment command script
    alignment_script_path = os.path.join(output, 'tmp_calculate_depth.sh')
    with open(alignment_script_path, 'w') as script:
        for index, read_path in enumerate(reads, start=1):
            ax_input = {
                'pacbio': 'map-pb',
                'nanopore': 'map-ont',
                'pacbio_hifi': 'map-hifi',
                'pacbio_asm20': 'asm20'
            }[sequencing_type]
            
            sam_path = os.path.join(output, f"All_gene_collections_mapped.{index}.sam")
            bam_path = sam_path.replace('.sam', '.bam')
            sorted_bam_path = bam_path.replace('.bam', '.sorted.bam')
            tmp_files_dir = os.path.join(output, f"sambamba_tmpfiles.{index}")

            script.write(
                f"minimap2 -ax {ax_input} {all_genes_path} {read_path} > {sam_path};\n"
                f"samtools view -bS {sam_path} > {bam_path} -@ {cpu_numbers};\n"
                f"mkdir -p {tmp_files_dir}; sambamba sort {bam_path} --tmpdir {tmp_files_dir} -o {sorted_bam_path};\n"
                f"samtools index {sorted_bam_path};\n"
                f"samtools flagstat {sorted_bam_path} > {sorted_bam_path.replace('.bam', '.stat')};\n"
                f"rm {sam_path} {bam_path}; rm -r {tmp_files_dir}\n"
            )

    # Execute the alignment script in parallel based on CPU number
    subprocess.run(["bash", alignment_script_path], check=True)
    os.remove(alignment_script_path)

    # Cleanup
    os.remove(all_genes_path)
    os.remove(f"{output}/*.bt2")
    for bam_file in glob.glob(f"{output}/*.bam"):
        os.remove(bam_file)
    for sorted_bam_file in glob.glob(f"{output}/*.sorted.bam"):
        os.remove(sorted_bam_file)

    return {'coverage_data': 'placeholder'}
 
def run_parallel(file, output_dir, cpu_numbers):
    
    logging.info(f"Running {file} in parallel with {cpu_numbers} threads...")
    
    log_file_path = os.path.join(output_dir, "run_parallel_log.txt")

    # Clear existing log
    if os.path.exists(log_file_path):
        os.remove(log_file_path)

    with open(file, 'r') as f:
        commands = [cmd for cmd in f.read().strip().split('\n') if cmd]

    def run_command(cmd):
        with open(log_file_path, 'a') as log_file:
            log_file.write(f"{datetime.datetime.now()} - Executing command: {cmd}\n")
            process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout, stderr = process.communicate()

            if stdout:
                log_file.write(f"stdout: {stdout.decode()}\n")
            if stderr:
                log_file.write(f"stderr: {stderr.decode()}\n")

            if process.returncode != 0:
                log_file.write(f"Error in command {cmd}\n")
            else:
                log_file.write("Command completed successfully.\n")
            log_file.write("\n")
        return cmd  # Return value for progress tracking

    # Show progress w/ tqdm
    with ThreadPoolExecutor(max_workers=cpu_numbers) as executor:
        futures = {executor.submit(run_command, cmd): cmd for cmd in commands}
        for _ in tqdm(as_completed(futures), total=len(futures), desc="Running in parallel"):
            pass

def show_download_progress(block_num, block_size, total_size):
    downloaded = block_num * block_size
    percent = min(100, downloaded * 100 / total_size) if total_size > 0 else 0
    progress_bar = f"[{'=' * int(percent // 2):<50}]"
    sys.stdout.write(f"\rDownloading: {progress_bar} {percent:.2f}%")
    sys.stdout.flush()

    if downloaded >= total_size:
        sys.stdout.write("\n")

def get_check_score(check_result_file):
    score = 0
    with open(check_result_file, 'r') as file:
        for line in file:
            if not line.startswith("#"):
                parts = line.split()
                score = max(score, float(parts[5]))  # Score is at index 5
    return score

def write_batched_scripts(stage_name, command_list, batch_dir, total_cpus, min_cpus_per_command=6):
    os.makedirs(batch_dir, exist_ok=True)
    num_cmds = len(command_list)
    cpus_per_cmd = max(min_cpus_per_command, total_cpus // max(num_cmds, 1))
    max_parallel_cmds = max(1, total_cpus // cpus_per_cmd)
    batches = math.ceil(num_cmds / max_parallel_cmds)

    script_paths = []
    for i in range(batches):
        batch_cmds = command_list[i * max_parallel_cmds:(i + 1) * max_parallel_cmds]
        script_path = os.path.join(batch_dir, f"{stage_name}_batch_{i+1}.sh")
        with open(script_path, "w") as f:
            f.write("\n".join(batch_cmds) + "\n")
        script_paths.append(script_path)
    return script_paths

def log_divider(emoji=None, title=None, position=None):
    if title and position == "top":
        logging.info("-" * 80)
        logging.info(f"{emoji} {title}")
    elif title and position == "bottom":
        logging.info(f"{emoji} {title}")
        logging.info("-" * 80)
    else:
        logging.info("-" * 80)

def repair_after_rename(output_dir, input_dir, threads):
    import os, csv, logging, shutil, pandas as pd, re

    rename_log_path = os.path.join(output_dir, "MAG_rename_log.csv")

    # Build rename_map from original base name (no extension) → new base name
    with open(rename_log_path) as f:
        reader = csv.DictReader(f)
        rename_map = {
            os.path.splitext(row['original_filename'])[0]: os.path.splitext(row['new_filename'])[0]
            for row in reader
        }

    def apply_rename(path, rename_map):
        if not os.path.isdir(path):
            logging.warning(f"Directory not found: {path}")
            return

        renamed = 0
        known_suffixes = [
            ".faa", ".gff", ".gene",
            ".result.txt", ".hits.txt",
            ".dbCAN2.out", ".dbCAN2.out.dm",
            ".MEROPSout.m8"
        ]

        for fname in os.listdir(path):
            base, ext = os.path.splitext(fname)
            match_found = False

            for suffix in known_suffixes:
                if fname.endswith(suffix):
                    core = fname.replace(suffix, "")
                    if core in rename_map:
                        new_core = rename_map[core]
                        new_fname = fname.replace(core, new_core, 1)
                        src = os.path.join(path, fname)
                        dst = os.path.join(path, new_fname)
                        if os.path.exists(src):
                            os.rename(src, dst)
                            logging.debug(f"Renamed: {fname} → {new_fname}")
                            renamed += 1
                            match_found = True
                            break  # Stop checking other suffixes

            if not match_found:
                logging.debug(f"No match for file: {fname}")

        logging.info(f"✅ Renamed {renamed} files in: {path}")

    def fix_headers_parallel(folder_path, output_dir, threads):
        script_path = os.path.join(output_dir, "tmp_fix_headers.sh")
        rename_log_path = os.path.join(output_dir, "MAG_rename_log.csv")

        # Regenerate rename_map for headers
        with open(rename_log_path) as f:
            reader = csv.DictReader(f)
            rows = list(reader)

        def build_fix_command(file_path, log_path):
            return f"python3 modules/fix_headers_helper.py '{file_path}' '{log_path}'"

        with open(script_path, "w") as sh:
            for fname in os.listdir(folder_path):
                if fname.endswith((".faa", ".gene")):
                    full_path = os.path.join(folder_path, fname)
                    cmd = build_fix_command(full_path, rename_log_path)
                    sh.write(cmd + "\n")

        run_parallel(script_path, output_dir, threads)
        os.remove(script_path)
        logging.info("✅ Parallel header fixing completed.")

    def fix_coverage_depth(filepath, rename_map):
        if not os.path.exists(filepath):
            logging.warning(f"Coverage file not found: {filepath}")
            return

        temp_path = filepath + ".tmp"
        made_changes = False

        pattern = re.compile(r'^(' + '|'.join(map(re.escape, rename_map)) + r')~~')

        with open(filepath, "r") as infile, open(temp_path, "w") as outfile:
            for line in infile:
                match = pattern.match(line)
                if match:
                    old_id = match.group(1)
                    new_id = rename_map.get(old_id, old_id)
                    new_line = line.replace(f"{old_id}~~", f"{new_id}~~", 1)
                    outfile.write(new_line)
                    made_changes = True
                else:
                    outfile.write(line)

        if made_changes:
            shutil.move(temp_path, filepath)
            logging.info(f"✅ Coverage file updated: {filepath}")
        else:
            os.remove(temp_path)
            logging.info(f"ℹ️ No changes made to coverage file: {filepath}")

    def fix_tsv(filepath):
        if os.path.exists(filepath):
            shutil.copy2(filepath, filepath + ".bak")
            df = pd.read_csv(filepath, sep="\t")
            if df.shape[1] > 0:
                df.iloc[:, 0] = df.iloc[:, 0].replace(rename_map)
                df.to_csv(filepath, sep="\t", index=False)
            else:
                logging.warning(f"No columns found in TSV: {filepath}")

    # Rename everything
    folders_main_input = ["faa", "gff", "gene"]
    folders_output_dir = [
        "KEGG_identifier_result",
        "intermediate_files/dbCAN2_Files",
        "intermediate_files/MEROPS_Files"
    ]

    for folder in folders_main_input:
        logging.info(f"Renaming files in {folder}...")
        apply_rename(os.path.join(input_dir, folder), rename_map)

    for folder in folders_output_dir:
        logging.info(f"Renaming files in {folder}...")
        apply_rename(os.path.join(output_dir, folder), rename_map)

    # Fix headers
    logging.info("Fixing headers in Each_HMM_Amino_Acid_Sequence...")
    fix_headers_parallel(os.path.join(output_dir, "Each_HMM_Amino_Acid_Sequence"), output_dir, threads)

    # Fix depth + metadata tables
    logging.info("Fixing headers in gene collection...")
    fix_coverage_depth(os.path.join(output_dir, "coverage", "All_gene_collections_mapped.depth.txt"), rename_map)

    logging.info("Fixing GTDBtk summary TSV...")
    fix_tsv(os.path.join(output_dir, "intermediate_files/gtdbtk_Genome_files/gtdbtk.bac120.summary.tsv"))

    logging.info("Fixing NCBI correction TSV...")
    fix_tsv(os.path.join(output_dir, "intermediate_files/gtdbtk_Genome_files/ncbi.bac120.correction.tsv"))

    logging.info("✅ Post-rename repair completed.")
    