import shutil
import subprocess
import os
import logging
import tarfile
import urllib.request
import gzip
import importlib
from .functions import run_parallel, show_download_progress

class SETUP_CHECK:

    def __init__(self, paths, threads):
        self.paths = paths
        self.threads = threads

    def check_setup(self):
        self.check_libs()
        self.check_dependencies()
        self.check_r_packages()
        self.check_template_files()
        self.check_database_files()

    def check_libs(self):
        required_packages = {
            "Bio": "biopython",
            "pandas": "pandas",
            "BCBio.GFF": "bcbio-gff",
            "tqdm": "tqdm",
            "openpyxl": "openpyxl",
            "boolean.py": "boolean.py"
        }

        for module_name, pip_name in required_packages.items():
            try:
                module = importlib.import_module(module_name)
                version = getattr(module, "__version__", "unknown")
                logging.info(f"{module_name} ({pip_name}) is installed, version: {version}")
            except ImportError:
                logging.warning(f"{module_name} is not installed. Install it using: pip install {pip_name}")
                
    def check_dependencies(self):
        tool_versions = {
            "hmmscan": ["hmmscan", "-h"],
            "prodigal": ["prodigal", "-v"],
            "sambamba": ["sambamba", "--version"],
            "coverm": ["coverm", "--version"],
            "Rscript": ["Rscript", "--version"],
            "diamond": ["diamond", "--version"],
            "samtools": ["samtools", "--version"],
            "bowtie2": ["bowtie2", "--version"],
            "minimap2": ["minimap2", "--version"], #only for long reads
            "gtdbtk": ["gtdbtk", "--version"],
        }

        logging.info("Checking required tools and their versions...")
        missing = []

        for tool, cmd in tool_versions.items():
            path = shutil.which(cmd[0])
            if path is None:
                missing.append(tool)
                continue

            try:
                result = subprocess.run(
                    cmd,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True
                )
                version_output = result.stdout.strip() or result.stderr.strip()
                logging.info(f"--- {tool} found: {path}")
                logging.info(f"   └─ Version info: {version_output.splitlines()[0]}")
            except Exception as e:
                logging.warning(f"⚠️ Could not determine version for {tool}: {e}")

        if missing:
            logging.error(f"❌ Missing required tools: {', '.join(missing)}")
            raise EnvironmentError(f"Missing tools: {', '.join(missing)}")

        # Check GTDBTK_DATA_PATH
        gtdbtk_data_path = os.environ.get("GTDBTK_DATA_PATH")
        if not gtdbtk_data_path:
            logging.error("❌ GTDBTK_DATA_PATH environment variable is not set.")
            raise EnvironmentError("GTDBTK_DATA_PATH is not set.")
        elif not os.path.isdir(gtdbtk_data_path):
            logging.error(f"❌ GTDBTK_DATA_PATH points to invalid path: {gtdbtk_data_path}")
            logging.error("Please set GTDBTK_DATA_PATH to a valid path.")
            logging.error("For Conda: conda env config vars set GTDBTK_DATA_PATH=/path/to/unarchived/gtdbtk/data")
            logging.error("For Pip: export GTDBTK_DATA_PATH=/path/to/unarchived/gtdbtk/data")
            logging.error("For more details, visit: https://ecogenomics.github.io/GTDBTk/installing/index.html")
            raise EnvironmentError(f"GTDBTK_DATA_PATH invalid: {gtdbtk_data_path}")
        else:
            logging.info(f"GTDBTK_DATA_PATH is set: {gtdbtk_data_path}")
            
    def check_r_packages(self):
        r_packages = [
            "ggplot2"
        ]

        for package in r_packages:
            try:
                subprocess.run(["Rscript", "-e", f"library({package})"], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                logging.info(f"{package} is installed.")
            except subprocess.CalledProcessError:
                logging.warning(f"{package} is not installed. Please install it using R.")
                
    def check_template_files(self):
        template_files = [
            "ko00002.json",
            "R_pathways.txt",
            "Sequential_transformations.tsv",
            "order_of_input_01.txt",
            "order_of_input_02.txt",
            "CAZy_map.txt",
            "MW-score_reaction_table.txt",
            "motif.txt",
            "motif.pair.txt",
            "hmm_table_template.txt",
            "hmm_table_template_2.txt",
            "kegg_module_step_db.txt"
        ]

        missing_files = []
        for fname in template_files:
            path = os.path.join(self.paths['templates'], fname)
            if not os.path.exists(path):
                missing_files.append(fname)

                if fname == "ko00002.json":
                    logging.warning(f"⚠️ {fname} is missing. This file must be downloaded manually from KEGG.")
                    logging.warning("➡️  Visit: https://www.kegg.jp/kegg-bin/get_htext?ko00002.keg")
                    logging.warning("👉  Click 'Download json', then save it as 'ko00002.json' in your templates directory.")

                elif fname in ["hmm_table_template.txt", "hmm_table_template_2.txt", "kegg_module_step_db.txt"]:
                    logging.warning(f"⚠️ {fname} is missing. Generating with setup_generate_templates.py...")
                    try:
                        subprocess.run(
                            ["python3", "modules/setup_generate_templates.py"],
                            check=True
                        )
                    except subprocess.CalledProcessError as e:
                        logging.error(f"❌ Template generation failed: {e}")
                        continue  # keep marked as missing

                    # re-check after generation
                    if os.path.exists(path):
                        logging.info(f"✅ {fname} successfully generated.")
                        if fname in missing_files:
                            missing_files.remove(fname)
                    else:
                        logging.error(f"❌ Failed to generate {fname}. Please check for errors.")

                else:
                    logging.warning(f"⚠️ {fname} is missing. Please ensure it is present in {self.paths['templates']}")

        if not missing_files:
            logging.info("All template files found.")
        else:
            logging.warning(f"⚠️ {len(missing_files)} template files missing. Check logs for details.")
            raise EnvironmentError(f"⚠️ {len(missing_files)} template files missing. Check logs for details.")
                
    def check_database_files(self):
        # === KEGG kofam_database ===
        kofam_dir = os.path.join(self.paths['databases'], 'kofam_database')
        profiles_dir = os.path.join(kofam_dir, 'profiles')
        done_marker_ko = os.path.join(kofam_dir, 'kofam_done')
        os.makedirs(kofam_dir, exist_ok=True)

        ko_list_path = os.path.join(kofam_dir, 'ko_list')
        profiles_tar_path = os.path.join(kofam_dir, 'profiles.tar.gz')

        if os.path.isdir(profiles_dir) and os.path.exists(ko_list_path) and os.path.exists(done_marker_ko):
            logging.info("kofam_database exists.")
        else:
            logging.warning("⚠️ kofam_database not found or incomplete.")
            logging.info("Downloading and setting up kofam_database...")

            # Download ko_list.gz
            ko_list_gz_url = "ftp://ftp.genome.jp/pub/db/kofam/ko_list.gz"
            ko_list_gz_path = ko_list_path + ".gz"
            urllib.request.urlretrieve(ko_list_gz_url, ko_list_gz_path, reporthook=show_download_progress)

            logging.info("Extracting ko_list.gz.")
            with gzip.open(ko_list_gz_path, 'rb') as f_in, open(ko_list_path, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
            os.remove(ko_list_gz_path)

            # Download and extract profiles
            logging.info("Downloading profiles.tar.gz...")
            profiles_tar_url = "ftp://ftp.genome.jp/pub/db/kofam/profiles.tar.gz"
            urllib.request.urlretrieve(profiles_tar_url, profiles_tar_path, reporthook=show_download_progress)

            logging.info("Extracting profiles.tar.gz.")
            with tarfile.open(profiles_tar_path, "r:gz") as tar:
                tar.extractall(kofam_dir)
            os.remove(profiles_tar_path)

            # Move All_Module_KO_ids.txt if it exists
            logging.info("Moving All_Module_KO_ids.txt to profiles directory...")
            alt_ko_path = os.path.join(self.paths['databases'], 'All_Module_KO_ids.txt')
            if os.path.exists(alt_ko_path):
                shutil.copy(alt_ko_path, profiles_dir)

            # Run hmmpress on all .hmm files in kofam profiles
            logging.info("Running hmmpress on kofam HMM profiles...")
            # Collect all hmmpress commands
            hmmpress_commands = []
            for hmm_file in os.listdir(profiles_dir):
                if hmm_file.endswith('.hmm'):
                    hmm_path = os.path.join(profiles_dir, hmm_file)
                    hmmpress_commands.append(f"hmmpress {hmm_path}")

            # Write commands to a temporary shell script
            script_path = os.path.join(self.paths['output_dir'], 'tmp_run_hmmpress.sh')
            with open(script_path, 'w') as script_file:
                script_file.write("\n".join(hmmpress_commands) + "\n")

            # Run commands in parallel
            run_parallel(script_path, self.paths['output_dir'], self.threads)
            os.remove(script_path)

            # Create marker file to signal successful setup
            with open(done_marker_ko, 'w') as f:
                f.write("kofam_database setup complete.\n")

            logging.info("kofam_database setup complete.")

        # === dbCAN2 HMM database ===
        dbcan_dir = os.path.join(self.paths['databases'], 'dbCAN2')
        done_marker_dbcan = os.path.join(dbcan_dir, 'dbCAN2_done')
        os.makedirs(dbcan_dir, exist_ok=True)


        dbcan_hmm = os.path.join(dbcan_dir, 'dbCAN-fam-HMMs.txt')
        if os.path.exists(dbcan_hmm) and os.path.exists(dbcan_hmm + ".h3f") and os.path.exists(done_marker_dbcan):
            logging.info("dbCAN2 HMM database already present and pressed.")
        else:
            logging.warning("⚠️ dbCAN2 HMM database not found or incomplete.")
            logging.info("Downloading dbCAN2 HMM file...")


            dbcan_url = "https://bcb.unl.edu/dbCAN2/download/Databases/dbCAN-old@UGA/dbCAN-fam-HMMs.txt.v14"
            tmp_file = dbcan_hmm + ".v14"
            urllib.request.urlretrieve(dbcan_url, tmp_file, reporthook=show_download_progress)


            # Rename to strip off .v14 and keep only .txt
            if os.path.exists(tmp_file):
                os.replace(tmp_file, dbcan_hmm)


            logging.info("Running hmmpress on dbCAN2 HMM file...")
            subprocess.run(['hmmpress', dbcan_hmm], check=True)


            # Create marker file to signal successful setup
            with open(done_marker_dbcan, 'w') as f:
                f.write("dbCAN2 setup complete.\n")


            logging.info("dbCAN2 HMM database setup complete.")
        
        # === MEROPS database ===
        MEROPS_db_base = os.path.join(self.paths['databases'], 'MEROPS', 'pepunit')
        MEROPS_db = MEROPS_db_base + ".dmnd"
        if not os.path.exists(MEROPS_db):
            logging.error(f"Missing MEROPS database: {MEROPS_db}, checking for pepunit.lib")

            lib_path = os.path.join(self.paths['databases'],'MEROPS', 'pepunit.lib')
            if not os.path.exists(lib_path):
                logging.error("pepunit.lib not found. Dowloading...")
                
                os.makedirs(os.path.join(self.paths['databases'],'MEROPS'), exist_ok=True)
                lib_url = "ftp://ftp.ebi.ac.uk/pub/databases/merops/current_release/pepunit.lib"
                urllib.request.urlretrieve(lib_url, lib_path, reporthook=show_download_progress)
                logging.info("pepunit.lib downloaded successfully. Generating pepunit.dmnd...")
            else:    
                logging.info(f"Found pepunit.lib, generating pepunit.dmnd...")

            tmp_seq = os.path.join(self.paths['databases'], 'MEROPS', 'tmp.seq')

            # Clean pepunit.lib and reformat
            seq_data = {}
            with open(lib_path, 'r', encoding='latin-1') as infile:
                for line in infile:
                    line = line.strip().replace('\r', '')
                    if line.startswith('>'):
                        header = line
                        seq_data[header] = ''
                    else:
                        seq_data[header] += line.replace(' ', '')

            with open(tmp_seq, 'w') as outfile:
                for h in sorted(seq_data):
                    outfile.write(f"{h}\n{seq_data[h]}\n")

            subprocess.run([
                'diamond', 'makedb', '--in', tmp_seq, '-d', MEROPS_db_base, '-p', str(self.threads),
            ], check=True)

            os.remove(tmp_seq)
            logging.info("pepunit.dmnd successfully created.")
        else:
            logging.info("pepunit.dmnd exists.")
        
        # === Currated HMM database ===
        # This is obtained from https://github.com/AnantharamanLab/METABOLIC/
        
        METABOLIC_hmm_db = os.path.join(self.paths['databases'], 'METABOLIC_hmm_db')
        if not os.path.exists(METABOLIC_hmm_db):
            logging.error(f"Missing METABOLIC_hmm_db: {METABOLIC_hmm_db}, checking for METABOLIC_hmm_db.tgz")
            
            METABOLIC_hmm_db_tgz = os.path.join(self.paths['databases'], 'METABOLIC_hmm_db.tgz')
            
            if not os.path.exists(METABOLIC_hmm_db_tgz):
                logging.error("METABOLIC_hmm_db.tgz not found. Downloading...")
                
                METABOLIC_hmm_db_url = "https://github.com/AnantharamanLab/METABOLIC/raw/refs/heads/master/METABOLIC_hmm_db.tgz"
                urllib.request.urlretrieve(METABOLIC_hmm_db_url, METABOLIC_hmm_db_tgz, reporthook=show_download_progress)
                logging.info("METABOLIC_hmm_db.tgz downloaded successfully. Decompressing...")
                with tarfile.open(METABOLIC_hmm_db_tgz, "r:gz") as tar:
                    tar.extractall(os.path.join(self.paths['databases']))
                os.remove(METABOLIC_hmm_db_tgz)
        else:
            logging.info("METABOLIC_hmm_db exists.")
