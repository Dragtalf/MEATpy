import os
import glob
import subprocess
import logging
from Bio import SeqIO
from .functions import run_parallel, write_batched_scripts

class CoverageCalculator:
    def __init__(self, paths, reads_csv, threads, reads_type, sequencing_type):
        self.paths = paths
        self.coverage_out = os.path.join(paths['output_dir'], 'coverage')
        self.reads_csv = reads_csv
        self.threads = threads
        self.reads_type = reads_type
        self.sequencing_type=sequencing_type
        
    def get_genome_coverage(self):
        """
        Get genome coverage percentages using Illumina or long reads and create Total.R_input.txt.
        """
        depth_file = os.path.join(self.coverage_out, 'All_gene_collections_mapped.depth.txt')
        if (not os.path.exists(depth_file) or os.path.getsize(depth_file) == 0) and self.reads_csv:

            os.makedirs(self.coverage_out, exist_ok=True)
            gene_output_path = os.path.join(self.coverage_out, 'All_gene_collections.gene')

            # Concatenate and reformat all .gene FASTAs
            with open(gene_output_path, "w") as out_f:
                for gene_file in sorted(glob.glob(os.path.join(self.paths['gene_dir'], "*.gene"))):
                    genome_id = os.path.basename(gene_file).replace(".gene", "")
                    for record in SeqIO.parse(gene_file, "fasta"):
                        record.id = f"{genome_id}~~{record.id.split()[0]}"
                        record.description = ""
                        SeqIO.write(record, out_f, "fasta")

            batch_dir = os.path.join(self.paths['output_dir'], "tmp_scripts")

            # Short-read (Illumina) coverage using Bowtie2
            if self.reads_type == "metaG":
                subprocess.run([
                    "bowtie2-build", "--quiet",
                    gene_output_path,
                    os.path.join(self.coverage_out, 'All_gene_collections.gene.scaffold')
                ], check=True)

                read_pairs = []
                with open(self.reads_csv) as f:
                    for line in f:
                        line = line.strip()
                        if line and not line.startswith("#"):
                            parts = line.split(",")
                            if len(parts) == 2:
                                read_pairs.append((parts[0], parts[1]))

                commands = {"bowtie": [], "view": [], "sort": [], "index": [], "cleanup": []}
                for i, (r1, r2) in enumerate(read_pairs, 1):
                    prefix = f"{self.coverage_out}/All_gene_collections_mapped.{i}"
                    tmpdir = f"{self.coverage_out}/sambamba_tmpfiles.{i}"

                    commands["bowtie"].append(
                        f"bowtie2 -x {self.coverage_out}/All_gene_collections.gene.scaffold -1 {r1} -2 {r2} -S {prefix}.sam -p {self.threads} --quiet"
                    )
                    commands["view"].append(
                        f"samtools view -bS {prefix}.sam > {prefix}.bam -@ {self.threads} 2> /dev/null"
                    )
                    commands["sort"].append(
                        f"mkdir -p {tmpdir}; sambamba sort {prefix}.bam --tmpdir {tmpdir} -o {prefix}.sorted.bam 2> /dev/null"
                    )
                    commands["index"].append(
                        f"samtools index {prefix}.sorted.bam 2> /dev/null"
                    )
                    commands["cleanup"].append(
                        f"rm {prefix}.sam {prefix}.bam; rm -r {tmpdir}"
                    )

            # Long-read (Nanopore/PacBio) coverage using minimap2
            elif self.paths["sequencing_type"] in ["pacbio", "pacbio_hifi", "nanopore"]:
                read_files = []
                with open(self.reads_csv) as f:
                    for line in f:
                        line = line.strip()
                        if line and not line.startswith("#"):
                            parts = line.split(",")
                            read_files.append(parts[0])

                # Choose preset based on technology
                preset = "map-ont" if "nano" in self.paths["sequencing_type"] else "map-pb"
                commands = {"map": [], "sort": [], "index": [], "cleanup": []}
                for i, read_fp in enumerate(read_files, 1):
                    prefix = f"{self.coverage_out}/All_gene_collections_mapped.{i}"
                    tmpdir = f"{self.coverage_out}/sambamba_tmpfiles.{i}"

                    commands["map"].append(
                        f"minimap2 -ax {preset} {gene_output_path} {read_fp} | samtools view -bS - > {prefix}.bam -@ {self.threads}"
                    )
                    commands["sort"].append(
                        f"mkdir -p {tmpdir}; sambamba sort {prefix}.bam --tmpdir {tmpdir} -o {prefix}.sorted.bam"
                    )
                    commands["index"].append(
                        f"samtools index {prefix}.sorted.bam"
                    )
                    commands["cleanup"].append(
                        f"rm {prefix}.bam; rm -r {tmpdir}"
                    )
            else:
                raise ValueError("Unknown reads type or sequencing_type for coverage calculation.")

            # Write and execute batches
            for stage, cmd_list in commands.items():
                script_batches = write_batched_scripts(stage, cmd_list, batch_dir, self.threads)
                for script in script_batches:
                    run_parallel(script, self.paths['output_dir'], self.threads)
                    os.remove(script)

            # Calculate coverage
            logging.info("Calculating coverage...")
            subprocess.run(
                f"coverm contig --methods metabat --bam-files {self.coverage_out}/All_gene_collections_mapped.*.sorted.bam > {depth_file} 2> /dev/null",
                shell=True, check=True
            )

            logging.info("Coverage calculation completed. Cleaning up temporary files...")
            subprocess.run(
                f"rm {self.coverage_out}/*.bt2 {self.coverage_out}/All_gene_collections.gene {self.coverage_out}/*.sorted.stat {self.coverage_out}/*.bai",
                shell=True
            )
            logging.info("Temporary files removed.")