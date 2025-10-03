import os
import logging
import pickle
import datetime
import glob
from version import __version__
from .setup_check import SETUP_CHECK
from .classify import Classifier
from .annotator import Annotator
from .hmm_utils import HMM_Manager
from .hmm_scan import HMMScanner, HMMCollector, Motif_Filter
from .kegg import KEGGAnnotator, KEGGProcessor
from .enzymes import EnzymePipeline
from .worksheets import *
from .coverage import CoverageCalculator
from .diagram_inputs import DiagramInputs
from .visuals import DiagramGenerator
from .functions import log_divider, repair_after_rename

import csv

class MEATpyRunner:
    def __init__(self, args):
        self.args = args

        self.paths = {
            'input_dir': args.input_folder,
            'output_dir': args.output_folder,
            'fasta_dir': os.path.join(args.input_folder, 'fasta'),
            'faa_dir': os.path.join(args.input_folder, 'faa'),
            'gff_dir': os.path.join(args.input_folder, 'gff'),
            'gene_dir': os.path.join(args.input_folder, 'gene'),
            'templates': os.path.join('templates'),
            'databases': os.path.join('databases'),
            'r_scripts': os.path.join('Rscripts')
        }

        self.annotator = Annotator(self.paths, args.prodigal_method, args.threads)
        self.hmm_manager = HMM_Manager(self.paths, args.kofam_database_size)
        self.total_hmm2threshold = {}
        self._setup_logging()

    def _setup_logging(self):
        os.makedirs(self.paths['output_dir'], exist_ok=True)
        log_path = os.path.join(self.paths['output_dir'], 'MEATpy.log')

        # Reset handlers
        for handler in logging.root.handlers[:]:
            logging.root.removeHandler(handler)

        logging.basicConfig(
            filename=log_path,
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            datefmt='%d-%b-%y %H:%M:%S'
        )

        # Make sure root logger has the correct level
        logging.getLogger().setLevel(logging.INFO)

        # Console handler
        console = logging.StreamHandler()
        console.setLevel(logging.INFO)
        plain_formatter = logging.Formatter('%(message)s')
        console.setFormatter(plain_formatter)
        logging.getLogger('').addHandler(console)

        logging.info("Logging initialized.")

    def run(self):
        log_divider()
        logging.info(f"MEATpy: Microbial Ecosystem Annotation Toolkit (Python) {__version__} — © 2025 BM Flatgard")
        log_divider()
        
        self.check_setup()
        
        if self.args.oops:
            log_divider("🤦‍♂️🤦🤦‍♀️", "Fixing downstream naming to avoid running the pipeline again...", "top")
            repair_after_rename(self.paths['output_dir'], self.paths['input_dir'], self.args.threads)
        elif self.args.oops and self.args.rename or self.args.revert_names:
            log_divider("❌", "You can't use --oops and --rename/--revert_names at the same time. Exiting.", "top")
            exit(1)
        
        if self.args.skip_to_diagrams:
            log_divider("🖼️", "Running with --skip_to_diagrams. Skipping to diagram generation.", "top")
            required_files = ['hmmscan_result.pkl', 'dbcan_results.pkl']
            missing = [f for f in required_files if not os.path.exists(os.path.join(self.paths['output_dir'], f))]

            if missing:
                logging.error(f"❌ Missing required file(s) for --skip_to_diagrams: {', '.join(missing)}")
                raise FileNotFoundError("Cannot continue. Please run MEATpy without --skip_to_diagrams to generate required files.")

            # Load hmmscan result
            try:
                path = os.path.join(self.paths['output_dir'], 'hmmscan_result.pkl')
                with open(path, 'rb') as f:
                    data = pickle.load(f)
                    self.hmmscan_result = data['hmmscan_result']
                    self.hmmscan_hits = data['hmmscan_hits']
                    self.hmm_ids = data['hmm_ids']
                log_divider("✅", "Loaded hmmscan_result.pkl successfully.", "top")
            except Exception as e:
                logging.error(f"❌ Failed to load hmmscan_result.pkl: {e}")
                raise

            # Load dbCAN results
            try:
                path = os.path.join(self.paths['output_dir'], 'dbcan_results.pkl')
                with open(path, 'rb') as f:
                    data = pickle.load(f)
                    self.dbcan_results = data['dbcan_results']
                log_divider("✅", "Loaded dbcan_results.pkl successfully.")
            except Exception as e:
                logging.error(f"❌ Failed to load dbcan_results.pkl: {e}")
                raise
            #self.rename_patch()
            self.get_diagram_inputs()
            self.generate_diagrams()
            self.finish()
        
        self.classify_genomes()
        self.annotate_genomes()
        self.prepare_thresholds()
        self.run_hmmsearch()
        self.summarize_results()
        self.generate_worksheets_1_2()
        self.kegg_annotator()
        self.generate_worksheets_3_4()
        self.process_kegg_hits()
        self.run_dbcan()
        self.generate_worksheet_5()
        self.run_merops()
        self.generate_worksheet6()
        self.combine_worksheets()
        self.get_coverage()
        self.get_diagram_inputs()
        self.generate_diagrams()
        self.finish()

    def check_setup(self):
        log_divider("🔧", "Checking setup...", "top")
        check = SETUP_CHECK(
            self.paths,
            self.args.threads
        )
        check.check_setup()
        log_divider("✅", "Setup check complete", "bottom")

    def rename_patch(self):
        
        classifier = Classifier(
            self.paths,
            self.args.threads,
            self.args.correction,
            self.args.revert_names
        )
        
        gtdbtk_dir = os.path.join(self.paths['output_dir'], 'intermediate_files', 'gtdbtk_Genome_files')
        gtdbtk_bac= os.path.join(gtdbtk_dir, "gtdbtk.bac120.summary.tsv")
        ncbi_bac= os.path.join(gtdbtk_dir, "ncbi.bac120.correction.tsv")
        gtdbtk_ar= os.path.join(gtdbtk_dir, "gtdbtk.ar53.summary.tsv")
        ncbi_ar= os.path.join(gtdbtk_dir, "ncbi.ar53.correction.tsv")
        
        log_path = os.path.join(self.paths['output_dir'], 'MAG_rename_log.csv')
        with open(log_path, 'r') as csvfile:
            reader = csv.DictReader(csvfile)
            rename_log = [(row['original_filename'], row['new_filename']) for row in reader]

        classifier.patch_summary_file(gtdbtk_bac,rename_log)
        classifier.patch_summary_file(ncbi_bac,rename_log)
        classifier.patch_summary_file(gtdbtk_ar,rename_log)
        classifier.patch_summary_file(ncbi_ar,rename_log)
        
    def classify_genomes(self):
        log_divider("🧬", "Starting genome classification...", "top")
        classifier = Classifier(
            self.paths,
            self.args.threads,
            self.args.correction,
            self.args.revert_names
            )
        classifier.classify_genomes()
        log_divider("✅", "Genome classification complete.", "bottom")
        if self.args.rename or self.args.revert_names:
            log_divider("🪪", "Renaming input files...", "top")
            rename_log, files = classifier.rename_inputs()
            if len(rename_log) == 0:
                logging.error("⚠️ No files were renamed. Check the input folder.")
                exit(1)
            elif len(rename_log) < len(files):
                logging.error("⚠️ Some files were not renamed. Check the input folder.")
                exit(1)
            else:
                log_divider("✅", "Input files renamed.", "bottom")

    def annotate_genomes(self):
        log_divider("🧪", "Starting genome annotation...", "top")
        self.annotator.run_prodigal_annotation()
        self.total_faa_path = self.annotator.concatenate_faa_files()
        log_divider("✅", "Genome annotation complete.", "bottom")

    def prepare_thresholds(self):
        log_divider("📥", "Loading HMM templates and thresholds...", "top")
        (self.hmm_table_temp,
         self.hmm_table_head,
         hmm_thresholds,
         self.hmm_table_temp_2,
         self.hmm_df,
         self.hmm_df2) = self.hmm_manager.load_hmm_templates()
        kofam_thresholds = self.hmm_manager.get_kofam_db_ko_threshold()
        self.total_hmm2threshold = self.hmm_manager.merge_thresholds(hmm_thresholds, kofam_thresholds)
        
        from .functions import get_motif, get_motif_pair
        self.motifs = get_motif(os.path.join(self.paths['templates'], 'motif.txt'))
        self.motif_pairs = get_motif_pair(os.path.join(self.paths['templates'], 'motif.pair.txt'))
        
        self.motif_filter = Motif_Filter(self.paths, self.args.threads)
        self.hmm_scanner = HMMScanner(
            self.paths,
            self.total_hmm2threshold,
            self.motifs,
            self.motif_pairs,
            self.motif_filter,
            self.args.threads
    )
        log_divider("✅", "HMM templates and thresholds loaded successfully.", "bottom")

    def run_hmmsearch(self):
        log_divider("⚙️", "Preparing and running HMM search...", "top")
        self.hmm_scanner.prepare_hmm_search_script()
        log_divider("✅", "HMM search finished.", "bottom")

    def summarize_results(self):
        log_divider("📊", "Summarizing HMM results...", "top")
        (self.genome_ids,
         self.total_faa_seq,
         self.total_gene_seq,
         self.seqid_to_genomeid) = self.hmm_scanner.prepare_hmm_search_data()

        (self.hmmscan_result,
         self.hmmscan_hits,
         self.hmm_ids) = self.hmm_scanner.summarize_hmmsearch_results(self.seqid_to_genomeid)
        
        with open(os.path.join(self.paths['output_dir'], 'hmmscan_result.pkl'), 'wb') as f:
            pickle.dump({
                'hmmscan_result': self.hmmscan_result,
                'hmmscan_hits': self.hmmscan_hits,
                'hmm_ids': self.hmm_ids
            }, f) 
        log_divider("✅", "HMM results summarized.", "bottom")

        log_divider("🧱", "Generating HMM collections...", "top")
        collector = HMMCollector(self.paths)
        collector.generate_hmm_collections(
            self.hmm_ids,
            self.hmmscan_hits,
            self.total_faa_seq,
            self.total_gene_seq)
        log_divider("✅", "HMM collections generated.", "bottom")

    def generate_worksheets_1_2(self):
        log_divider("📄", "Generating Worksheets 1 and 2...", "top")
        logging.debug(f"results: {self.hmmscan_result}")
        logging.debug(f"hits: {self.hmmscan_hits}")
        generate_worksheet1(
            self.paths['output_dir'],
            self.genome_ids,
            self.hmm_table_head,
            self.hmm_table_temp,
            self.hmmscan_result,
            self.hmmscan_hits
        )
        generate_worksheet2(
            self.paths['output_dir'],
            self.genome_ids,
            self.hmm_table_head,
            self.hmm_table_temp_2,
            self.hmm_table_temp,
            self.hmmscan_result
        )
        log_divider("💾", "Worksheets 1 and 2 written.", "bottom")

    def kegg_annotator(self):
        log_divider("🧭", "Starting KEGG annotation...", "top")
        annotator = KEGGAnnotator(
            paths=self.paths,
            hmm_table_temp=self.hmm_table_temp,
            genome_ids=self.genome_ids,
            hmmscan_result=self.hmmscan_result,
            module_cutoff=self.args.module_cutoff
        )
        (
            self.category_modules,
            self.kegg_module,
            self.kegg_module2name,
            self.module_presence_by_genome,
            self.step_presence_by_genome,
            self.hmm2ko
        ) = annotator.run()
        log_divider("✅", "KEGG annotation complete.", "bottom")

    def generate_worksheets_3_4(self):
        log_divider("📄", "Generating Worksheets 3 and 4...", "top")
        generate_worksheet3(
            self.paths['output_dir'],
            self.genome_ids,
            self.kegg_module2name,
            self.category_modules,
            self.module_presence_by_genome
        )
        generate_worksheet4(
            self.paths['output_dir'],
            self.genome_ids,
            self.step_presence_by_genome,
            self.kegg_module,
            self.kegg_module2name,
            self.category_modules
        )
        log_divider("💾", "Worksheets 3 and 4 written.", "bottom")

    def process_kegg_hits(self):
        log_divider("🧬", "Processing KEGG hits...")
        kegg_processor = KEGGProcessor(
            self.paths,
            self.hmmscan_result,
            self.hmmscan_hits,
            self.hmm2ko,
            self.genome_ids,
            kofam_db_size=self.args.kofam_database_size
        )
        kegg_processor.run()
        log_divider("✅", "KEGG hits processed.", "bottom")

    def run_dbcan(self):
        log_divider("🍭", "Running dbCAN2...", "top")
        dbcan_runner = EnzymePipeline(
            self.paths,
            self.args.threads
        )
        self.dbcan_results, self.dbcan_hits = dbcan_runner.run()
        
        with open(os.path.join(self.paths['output_dir'], 'dbcan_results.pkl'), 'wb') as f:
            pickle.dump({
                'dbcan_results': self.dbcan_results
            }, f)
        log_divider("✅", "dbCAN2 searching completed.", "bottom")

    def generate_worksheet_5(self):
        log_divider("📄", "Generating worksheet 5...", "top")
        generate_worksheet5(
            self.paths['output_dir'],
            self.genome_ids,
            self.dbcan_results,
            self.dbcan_hits
        )
        log_divider("💾", "Worksheet 5 created.", "bottom")

    def run_merops(self):
        log_divider("🔪", "Running MEROPS...", "top")
        merops_runner = EnzymePipeline(
            self.paths,
            self.args.threads
        )
        self.merops_out, self.formatted_hits, self.merops_ids = merops_runner.run_merops()
        log_divider("✅", "MEROPS searching completed.", "bottom")

    def generate_worksheet6(self):
        log_divider("📄", "Generating worksheet 6...", "top")
        generate_worksheet6(
            self.paths['output_dir'],
            self.genome_ids,
            self.merops_out,
            self.formatted_hits,
            self.merops_ids
        )
        log_divider("💾", "Worksheet 6 created.", "bottom")
        
    def combine_worksheets(self):
        log_divider("📚", "Combining worksheets into Excel...", "top")
        combine_worksheets_to_excel(
            self.paths['output_dir']
            )
        log_divider("✅", "Worksheets combined into Excel.", "bottom")

    def get_coverage(self):
        log_divider("📏", "Calculating genome coverage...", "top")
        coverage_calculator = CoverageCalculator(
            self.paths,
            self.args.input_csv,
            self.args.threads,
            self.args.reads_type,
            self.args.sequencing_type
        )
        coverage_calculator.get_genome_coverage()
        log_divider("✅", "Genome coverage calculated.", "bottom")

    def get_diagram_inputs(self):
        log_divider("🖼️", "Generating diagram inputs...", "top")
        diagram_inputs = DiagramInputs(
            self.paths, 
            self.hmmscan_result,
            self.dbcan_results,
            self.args.threads,
            self.args.correction,
            self.args.taxonomy,
            self.args.input_csv
        )
        diagram_inputs.run()
        log_divider("✅", "All diagram inputs generated.", "bottom")

    def generate_diagrams(self):
        log_divider("📈", "Generating all diagrams...", "top")
        diagrams= DiagramGenerator(
            self.paths
        )
        diagrams.run()
        log_divider("✅", "All diagrams saved.", "bottom")
    
    def finish(self):
        log_divider("🧹", "Finishing up and removing temporary files...", "top")

        # Remove temporary shell scripts
        sh_files = glob.glob(os.path.join(self.paths['output_dir'], '*.sh'))
        for f in sh_files:
            try:
                os.remove(f)
                logging.debug(f"Removed temp file: {f}")
            except Exception as e:
                logging.warning(f"Could not remove file {f}: {e}")

        # Remove parallel log file if exists
        parallel_log = os.path.join(self.paths['output_dir'], 'run_parallel_log.txt')
        if os.path.exists(parallel_log):
            try:
                os.remove(parallel_log)
            except Exception as e:
                logging.warning(f"Could not remove {parallel_log}: {e}")

        # Write version record
        version_path = os.path.join(self.paths['output_dir'], 'MEATpy.version')
        with open(version_path, 'w') as vf:
            vf.write(f"MEATpy version: {__version__}\n")
            vf.write(f"Run completed on: {datetime.datetime.now()}")

        log_divider("🏁🏁🏁", "Pipeline complete!", "bottom")
        exit(0)
