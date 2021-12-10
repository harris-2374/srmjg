"""
Author: Andrew Harris
Email: ajharris@cvm.tamu.edu
Version: 0.0.5
Last update: 10/21/2021
"""
import argparse
import configparser
import logging
import os
import sys
from pathlib import Path
## Dependencies
import pandas as pd

########################### Logging Init ###########################
def set_logger(filename, OUTPUT, LOG_LEVEL):
    logger = logging.getLogger(__name__)
    logger.setLevel(LOG_LEVEL)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s: %(message)s', "%Y-%m-%d %H:%M:%S")
    log_output_dir = OUTPUT / "logs"
    log_output_dir.mkdir(parents=True, exist_ok=True)
    log_file_name = log_output_dir / filename
    file_handler = logging.FileHandler(log_file_name)
    stream_handler = logging.StreamHandler()
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)
    return logger
 
########################## Error Classes ##########################
class Error(Exception):
    """Base class for other exceptions"""
    pass
class InvalidHeader(Error):
    """Raise when not input is provided"""
    pass
class NoModulesFound(Error):
    """Raise when no modules are given for SLURM module load"""
    pass
class SlurmOutDirInvalid(Error):
    "Raise if no project output provided"
    pass
class SlurmTMPSizeInvalidType(Error):
    "Raise if TMPDIR type invalid"
    pass
class SlurmTMPSizeInvalidSize(Error):
    "Raise if TMPDIR size invalid"
    pass
########################## Helper Functions ##########################
# --- General Functions ---
def convert_log_level(l):
    valid_types = {'NOTSET': 0, 'DEBUG': 10, 'INFO': 20, 'WARNING': 30, 'ERROR': 40, 'CRITICAL': 50}
    if l in valid_types.keys():
        return valid_types[l]
    else:
        return False


def get_input_df(INPUT_PATH, logger):
    """Read input file into pandas dataframe, 
    ensure headers are valid, and return dataframe

    :param INPUT_PATH: Path to input file
    :type INPUT_PATH: Object
    :raises InvalidHeader: Header provided does not match required header values
    :return: dataframe with input file information
    :rtype: Pandas DataFrame
    """
    # Read file into DF
    if 'xls' in INPUT_PATH.name:
        try:
            import openpyxl
        except ImportError or ModuleNotFoundError:
            # The blabla module does not exist, display proper error message and exit
            logger.info(f"Module openpyxl not installed. Install and rerun")
            exit(1)
        df = pd.read_excel(INPUT_PATH, engine='openpyxl')
    elif 'csv' in INPUT_PATH.name:
        df = pd.read_csv(INPUT_PATH)
    elif 'tsv' in INPUT_PATH.name:
        df = pd.read_csv(INPUT_PATH, sep='\t')
    elif 'txt' in INPUT_PATH.name:
        df = pd.read_csv(INPUT_PATH, sep='\t')
    # Check headers are valid
    try:
        header_list = ['QueryID', 'QueryLibID', 'QueryRun', 'QueryR1', 'QueryR2', 'ReferenceID', 'ReferencePath']
        if df.columns.to_list() != header_list:
            raise InvalidHeader
    except InvalidHeader:
        print("Input file has invalid column headers - please verify and rerun")
    return df


def project_dir_setup_script(PDIR, OUTPUT_PATH, logger):
    """Writes file that creates requried directories for job files

    :param OUTPUT_PATH: Path to output directory
    :type OUTPUT_PATH: pathlib Path object
    """
    bam_dir = PDIR / "BAM"
    command = f"mkdir {bam_dir}"
    filename = OUTPUT_PATH / "project_dir_setup.sh"    
    with open(filename, 'w') as oh:
        oh.write(f"{command}")
    logger.info("-------------")
    logger.info(f"Project directory setup file written: {filename}")
    logger.info("-------------")
    return

# --- Config Template Functions ---
def slurm_config_template(outFile):
    with open(outFile, 'w') as oh:
        oh.write("""[SLURM INPUT]\nproject_directory = project/dir/path\njob_type = scsn\ntime = 7:00:00\nqueue = medium\nnodes = 1\ntasks_per_node = 1\nntasks = 1\ncpus_per_task = 48\nmemory = 10G\naccount = \nemail = \ntmp = \nmodules = """)
        return


def bash_config_template(outFile):
    with open(outFile, 'w') as oh:
        oh.write("""[BASH INPUT]\nproject_directory = project/dir/path\nthreads = 1""")
        return


# --- SLURM Functions ---
def validate_slurm_input(
    logger,
    SLURM_OUTDIR,  # Done
    SLURM_JOBTYPE,
    SLURM_TIME,
    SLURM_NODES,
    SLURM_TASKS_PER_NODE,
    SLURM_NTASKS,
    SLURM_CPUS,
    SLURM_MEMORY,
    SLURM_ACCOUNT,
    SLURM_EMAIL,
    SLURM_TMPDIR,
    SLURM_MODULES,
):
    """Ensure SLURM input values are valid"""
    try:
        if not SLURM_OUTDIR:
            raise SlurmOutDirInvalid
        else:
            pass
    except SlurmOutDirInvalid:
        logger.info("No project output directory provided. Ensure one is provided and rerun")

    return True


def generate_SLRUM_JobFiles(
    row,
    OUTPUT_PATH,
    PIPELINE,
    SLURM_OUTDIR,
    SLURM_BAM_OUTDIR,
    SLURM_JOBTYPE,
    SLURM_TIME,
    SLURM_NODES,
    SLURM_TASKS_PER_NODE,
    SLURM_NTASKS,
    SLURM_CPUS,
    SLURM_MEMORY,
    SLURM_ACCOUNT,
    SLURM_EMAIL,
    SLURM_TMPDIR,
    SLURM_MODULES,
    logger,
):
    def slurm_header(
        sampleID,
        refID,
        SLURM_JOBTYPE,
        SLURM_TIME,
        SLURM_NODES,
        SLURM_TASKS_PER_NODE,
        SLURM_NTASKS,
        SLURM_MEMORY,
        SLURM_ACCOUNT,
        SLURM_EMAIL,
        SLURM_MODULES,
        SLURM_TMPDIR,
    ):
        header = []
        if SLURM_JOBTYPE == 'scsn':
            """
            #!/bin/bash
            ##NECESSARY JOB SPECIFICATIONS
            #SBATCH --job-name=JobExample1       #Set the job name to "JobExample1"
            #SBATCH --time=01:30:00              #Set the wall clock limit to 1hr and 30min
            #SBATCH --ntasks=1                   #Request 1 task
            #SBATCH --mem=2560M                  #Request 2560MB (2.5GB) per node
            #SBATCH --output=Example1Out.%j      #Send stdout/err to "Example1Out.[jobID]"

            ## OPTIONAL JOB SPECIFICATIONS
            ##SBATCH --account=123456             #Set billing account to 123456
            ##SBATCH --mail-type=ALL              #Send email on all job events
            ##SBATCH --mail-user=email_address    #Send all emails to email_address
            """
            # -- Required headers --
            header.append("#!/bin/bash")
            header.append(f"#SBATCH --job-name={sampleID}_to_{refID}")
            header.append(f"#SBATCH --time={SLURM_TIME}")
            header.append(f"#SBATCH --ntasks={SLURM_NTASKS}")
            header.append(f"#SBATCH --mem={SLURM_MEMORY}")
            header.append(f"#SBATCH --output={sampleID}_to_{refID}.%j")
            header.append(f"#SBATCH --error={sampleID}_to_{refID}.%j")
            header.append(f"\n## OPTIONAL JOB SPECIFICATIONS")
            if len(SLURM_TIME.split("-")) > 1:
                split = SLURM_TIME.split("-")
                if int(split[0]) > 7:
                    header.append("#SBATCH --partition xlong")
                else:
                    pass
            if SLURM_ACCOUNT:
                header.append(f"#SBATCH --account={SLURM_ACCOUNT}")
            if SLURM_EMAIL:
                header.append("#SBATCH --mail-type=ALL")
                header.append(f"#SBATCH --mail-user={SLURM_EMAIL}")
            if SLURM_TMPDIR:
                try:
                    if type(SLURM_TMPDIR) != int:
                        raise SlurmTMPSizeInvalidType
                    elif SLURM_TMPDIR < 10240:
                        raise SlurmTMPSizeInvalidSize
                except SlurmTMPSizeInvalidType:
                    logger.info("SLURM_TMPDIR input type is invalid")
                    exit(1)
                except SlurmTMPSizeInvalidSize:
                    logger.info("SLURM_TMPDIR size must be larger than 10240 Mb")
                    exit(1)
                header.append("#SBATCH --tmp={SLURM_TMPDIR}")
            # -- Add Module Loading --
            header.append('')
            SLURM_MODULES = "\n".join(SLURM_MODULES.split(";"))
            header.append(f"{SLURM_MODULES}")
            return header
        elif SLURM_JOBTYPE == 'mcsn':
            """
            #!/bin/bash
            ##NECESSARY JOB SPECIFICATIONS
            #SBATCH --job-name=JobExample2       #Set the job name to "JobExample2"
            #SBATCH --time=6:30:00               #Set the wall clock limit to 6hr and 30min
            #SBATCH --nodes=1                    #Request 1 node
            #SBATCH --ntasks-per-node=8          #Request 8 tasks/cores per node
            #SBATCH --mem=8G                     #Request 8GB per node 
            #SBATCH --output=Example2Out.%j      #Send stdout/err to "Example2Out.[jobID]" 

            ## OPTIONAL JOB SPECIFICATIONS
            ##SBATCH --account=123456             #Set billing account to 123456
            ##SBATCH --mail-type=ALL              #Send email on all job events
            ##SBATCH --mail-user=email_address    #Send all emails to email_address 
            """
            # -- Required headers --
            header.append("#!/bin/bash")
            header.append(f"#SBATCH --job-name={sampleID}_to_{refID}")
            header.append(f"#SBATCH --time={SLURM_TIME}")
            header.append(f"#SBATCH --nodes=1")
            header.append(f"#SBATCH --ntasks-per-node={SLURM_TASKS_PER_NODE}")
            header.append(f"#SBATCH --mem={SLURM_MEMORY}")
            header.append(f"#SBATCH --output={sampleID}_to_{refID}.%j")
            header.append(f"#SBATCH --error={sampleID}_to_{refID}.%j")
            header.append(f"\n## OPTIONAL JOB SPECIFICATIONS")
            if len(SLURM_TIME.split("-")) > 1:
                split = SLURM_TIME.split("-")
                if int(split[0]) > 7:
                    header.append("#SBATCH --partition xlong")
                else:
                    pass
            if SLURM_ACCOUNT:
                header.append(f"#SBATCH --account={SLURM_ACCOUNT}")
            if SLURM_EMAIL:
                header.append("#SBATCH --mail-type=ALL")
                header.append(f"#SBATCH --mail-user={SLURM_EMAIL}")
            if SLURM_TMPDIR:
                try:
                    if type(SLURM_TMPDIR) != int:
                        raise SlurmTMPSizeInvalidType
                    elif SLURM_TMPDIR < 10240:
                        raise SlurmTMPSizeInvalidSize
                except SlurmTMPSizeInvalidType:
                    logger.info("SLURM_TMPDIR input type is invalid")
                    exit(1)
                except SlurmTMPSizeInvalidSize:
                    logger.info("SLURM_TMPDIR size must be larger than 10240 Mb")
                    exit(1)
                header.append("#SBATCH --tmp={SLURM_TMPDIR}")
            # -- Add Module Loading --
            header.append('\n')
            SLURM_MODULES = "\n".join(SLURM_MODULES.split(";"))
            header.append(f"{SLURM_MODULES}")
            return header
        elif SLURM_JOBTYPE == 'mcmn':
            """
            #!/bin/bash
            ##NECESSARY JOB SPECIFICATIONS
            #SBATCH --job-name=JobExample3       #Set the job name to "JobExample3"
            #SBATCH --time=1-12:00:00            #Set the wall clock limit to 1 Day and 12hr
            #SBATCH --ntasks=8                   #Request 8 tasks
            #SBATCH --ntasks-per-node=2          #Request 2 tasks/cores per node
            #SBATCH --mem=4096M                  #Request 4096MB (4GB) per node 
            #SBATCH --output=Example3Out.%j      #Send stdout/err to "Example3Out.[jobID]"
            
            ## OPTIONAL JOB SPECIFICATIONS
            ##SBATCH --account=123456             #Set billing account to 123456
            ##SBATCH --mail-type=ALL              #Send email on all job events
            ##SBATCH --mail-user=email_address    #Send all emails to email_address 
            """
            # -- Required headers --
            header.append("#!/bin/bash")
            header.append(f"#SBATCH --job-name={sampleID}_to_{refID}")
            header.append(f"#SBATCH --time={SLURM_TIME}")
            header.append(f"#SBATCH --ntasks={SLURM_NTASKS}")
            header.append(f"#SBATCH --ntasks-per-node={SLURM_TASKS_PER_NODE}")
            header.append(f"#SBATCH --mem={SLURM_MEMORY}")
            header.append(f"#SBATCH --output={sampleID}_to_{refID}.%j")
            header.append(f"#SBATCH --error={sampleID}_to_{refID}.%j")
            header.append(f"\n## OPTIONAL JOB SPECIFICATIONS")
            if len(SLURM_TIME.split("-")) > 1:
                split = SLURM_TIME.split("-")
                if int(split[0]) > 7:
                    header.append("#SBATCH --partition xlong")
                else:
                    pass
            if SLURM_ACCOUNT:
                header.append(f"#SBATCH --account={SLURM_ACCOUNT}")
            if SLURM_EMAIL:
                header.append("#SBATCH --mail-type=ALL")
                header.append(f"#SBATCH --mail-user={SLURM_EMAIL}")
            if SLURM_TMPDIR:
                try:
                    if type(SLURM_TMPDIR) != int:
                        raise SlurmTMPSizeInvalidType
                    elif SLURM_TMPDIR < 10240:
                        raise SlurmTMPSizeInvalidSize
                except SlurmTMPSizeInvalidType:
                    logger.info("SLURM_TMPDIR input type is invalid")
                    exit(1)
                except SlurmTMPSizeInvalidSize:
                    logger.info("SLURM_TMPDIR size must be larger than 10240 Mb")
                    exit(1)
                header.append("#SBATCH --tmp={SLURM_TMPDIR}")
            # -- Add Module Loading --
            header.append('\n')
            SLURM_MODULES = "\n".join(SLURM_MODULES.split(";"))
            header.append(f"{SLURM_MODULES}")
            return header
        elif SLURM_JOBTYPE == 'custom':
            """
            #!/bin/bash
            ##NECESSARY JOB SPECIFICATIONS
            #SBATCH --job-name=JobExample3       #Set the job name to "JobExample3"
            #SBATCH --time=1-12:00:00            #Set the wall clock limit to 1 Day and 12hr
            #SBATCH --ntasks=8                   #Request 8 tasks
            #SBATCH --ntasks-per-node=2          #Request 2 tasks/cores per node
            #SBATCH --mem=4096M                  #Request 4096MB (4GB) per node 
            #SBATCH --output=Example3Out.%j      #Send stdout/err to "Example3Out.[jobID]"
            
            ## OPTIONAL JOB SPECIFICATIONS
            ##SBATCH --account=123456             #Set billing account to 123456
            ##SBATCH --mail-type=ALL              #Send email on all job events
            ##SBATCH --mail-user=email_address    #Send all emails to email_address 
            """
            # -- Required headers --
            header.append("#!/bin/bash")
            header.append(f"#SBATCH --job-name={sampleID}_to_{refID}")
            header.append(f"#SBATCH --time={SLURM_TIME}")
            header.append(f"#SBATCH --ntasks={SLURM_NTASKS}")
            header.append(f"#SBATCH --ntasks-per-node={SLURM_TASKS_PER_NODE}")
            header.append(f"#SBATCH --mem={SLURM_MEMORY}")
            header.append(f"#SBATCH --output={sampleID}_to_{refID}.%j")
            header.append(f"#SBATCH --error={sampleID}_to_{refID}.%j")
            header.append(f"\n## OPTIONAL JOB SPECIFICATIONS")
            if len(SLURM_TIME.split("-")) > 1:
                split = SLURM_TIME.split("-")
                if int(split[0]) > 7:
                    header.append("#SBATCH --partition xlong")
                else:
                    pass
            if SLURM_ACCOUNT:
                header.append(f"#SBATCH --account={SLURM_ACCOUNT}")
            if SLURM_EMAIL:
                header.append("#SBATCH --mail-type=ALL")
                header.append(f"#SBATCH --mail-user={SLURM_EMAIL}")
            if SLURM_TMPDIR:
                try:
                    if type(SLURM_TMPDIR) != int:
                        raise SlurmTMPSizeInvalidType
                    elif SLURM_TMPDIR < 10240:
                        raise SlurmTMPSizeInvalidSize
                except SlurmTMPSizeInvalidType:
                    logger.info("SLURM_TMPDIR input type is invalid")
                    exit(1)
                except SlurmTMPSizeInvalidSize:
                    logger.info("SLURM_TMPDIR size must be larger than 10240 Mb")
                    exit(1)
                header.append("#SBATCH --tmp={SLURM_TMPDIR}")
            # -- Add Module Loading --
            header.append('\n')
            SLURM_MODULES = "\n".join(SLURM_MODULES.split(";"))
            header.append(f"{SLURM_MODULES}")
            return header

    sampleID, sampleLibID, sampleRun, sampleR1, sampleR2, refID, refPath = row
    output_filename = OUTPUT_PATH / f"{sampleID}_to_{refID}.sh"
    # -- Generate header --
    header = slurm_header(
        sampleID,
        refID,
        SLURM_JOBTYPE,
        SLURM_TIME,
        SLURM_NODES,
        SLURM_TASKS_PER_NODE,
        SLURM_NTASKS,
        SLURM_MEMORY,
        SLURM_ACCOUNT,
        SLURM_EMAIL,
        SLURM_MODULES,
        SLURM_TMPDIR,
    )
    header_formatted = "\n".join(header)
    # -- Generate command --
    if PIPELINE == 'bwa-mem2':
        logger.debug(f"Running {PIPELINE} Pipeline")
        sample_unsorted_sam = SLURM_BAM_OUTDIR / f"{sampleID}_to_{refID}.sam"
        sample_unsorted_bam = SLURM_BAM_OUTDIR / f"{sampleID}_to_{refID}.bam"
        sample_markdups_bam = SLURM_BAM_OUTDIR / f"{sampleID}_to_{refID}.sort.md.bam"
        sample_markdups_metrics = SLURM_BAM_OUTDIR / f"{sampleID}_to_{refID}.sort.md.metrics.txt"
        sample_completion_file = SLURM_OUTDIR / f"{sampleID}_to_{refID}.complete"
        if SLURM_TMPDIR:
            command = Rf"""bwa-mem2 mem -t {SLURM_CPUS} -Y -R "@RG\tID:{sampleID}\tPL:ILLUMINA\tLB:{sampleLibID}\tDS:{sampleID}_{sampleRun}\tPU:{sampleID}_{sampleRun}\tSM:{sampleID}" -o $TMPDIR/{sample_unsorted_sam.name} {refPath} {sampleR1} {sampleR2} && samtools view -bS $TMPDIR/{sample_unsorted_sam.name} > $TMPDIR/{sample_unsorted_bam.name} && gatk MarkDuplicatesSpark -I $TMPDIR/{sample_unsorted_bam.name} -O {sample_markdups_bam.as_posix()} -M {sample_markdups_metrics.as_posix()} --conf 'spark.executor.cores={SLURM_CPUS}' && touch {sample_completion_file.as_posix()}"""
        else:
            command = Rf"""bwa-mem2 mem -t {SLURM_CPUS} -Y -R "@RG\tID:{sampleID}\tPL:ILLUMINA\tLB:{sampleLibID}\tDS:{sampleID}_{sampleRun}\tPU:{sampleID}_{sampleRun}\tSM:{sampleID}" -o {sample_unsorted_sam.as_posix()} {refPath} {sampleR1} {sampleR2} && samtools view -bS {sample_unsorted_sam.as_posix()} > {sample_unsorted_bam.as_posix()} && gatk MarkDuplicatesSpark -I {sample_unsorted_bam.as_posix()} -O {sample_markdups_bam.as_posix()} -M {sample_markdups_metrics.as_posix()} --conf 'spark.executor.cores={SLURM_CPUS}' && touch {sample_completion_file.as_posix()}"""
        logger.debug(f"{command}")
    return output_filename, header_formatted, command


# --- BASH Functions ---
def generate_BASH_JobFiles(
    row,
    OUTPUT_PATH,
    BASH_OUTDIR,
    BASH_BAM_OUTDIR,
    BASH_CPU,
    PIPELINE,
    logger,
):
    sampleID, sampleLibID, sampleRun, sampleR1, sampleR2, refID, refPath = row
    output_filename = OUTPUT_PATH / f"{sampleID}_to_{refID}.sh"
    # -- Generate command --
    if PIPELINE == 'bwa-mem2':
        logger.debug(f"Running {PIPELINE} Pipeline")
        sample_unsorted_sam = BASH_BAM_OUTDIR / f"{sampleID}_to_{refID}.sam"
        sample_unsorted_bam = BASH_BAM_OUTDIR / f"{sampleID}_to_{refID}.bam"
        sample_markdups_bam = BASH_BAM_OUTDIR / f"{sampleID}_to_{refID}.sort.md.bam"
        sample_markdups_metrics = BASH_BAM_OUTDIR / f"{sampleID}_to_{refID}.sort.md.metrics.txt"
        sample_completion_file = BASH_OUTDIR / f"{sampleID}_to_{refID}.complete"
        command = Rf"""bwa-mem2 mem -t {BASH_CPU} -Y -R "@RG\tID:{sampleID}\tPL:ILLUMINA\tLB:{sampleLibID}\tDS:{sampleID}_{sampleRun}\tPU:{sampleID}_to_{refID}\tSM:{sampleID}" -o {sample_unsorted_sam.as_posix()} {refPath} {sampleR1} {sampleR2} && samtools view -bS {sample_unsorted_sam.as_posix()} > {sample_unsorted_bam.as_posix()} && gatk MarkDuplicatesSpark -I {sample_unsorted_bam.as_posix()} -O {sample_markdups_bam.as_posix()} -M {sample_markdups_metrics.as_posix()} --conf 'spark.executor.cores={BASH_CPU}' && touch {sample_completion_file.as_posix()}"""
        logger.debug(f"{command}")
    return output_filename, command


########################### Main Function ###########################
def main():
    # General input variable that are passed to both SLURM and BASH options
    parser = argparse.ArgumentParser(description='')
    parser.add_argument(
        '-i',
        '--input',
        type=str,
        action='store',
        default=None,
        help=".xlsx, .tsv, or .csv file - Required column headers = [sampleID | sampleRun | sampleR1 | sampleR2 | ReferenceID | ReferencePath]",
        metavar='',
    )
    parser.add_argument(
        '-o',
        '--output',
        type=str,
        action='store',
        default=os.getcwd(),
        help='Output location - default: cwd',
        metavar='',
    )
    parser.add_argument(
        '--pipeline',
        type=str,
        choices=['bwa-mem2'],
        default='bwa-mem2',
        help='Mapping pipeline',
    )
    parser.add_argument(
        '--scheduler',
        type=str,
        choices=['SLURM', 'BASH'],
        help='Type of job file to create',
    )
    parser.add_argument(
        '--log',
        type=str,
        action='store',
        default='ShortReadMappingPipeline.log',
        help='Name of log file',
        metavar='',
    )
    parser.add_argument(
        '--log_level',
        type=str,
        action='store',
        choices=['NOTSET', 'DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
        default='INFO',
        help='Log level',
    )
    parser.add_argument(
        '--slurm_config_template',
        action='store_true',
        default=False,
        help='Output SLURM parameter input config file',
    )
    parser.add_argument(
        '--bash_config_template',
        action='store_true',
        default=False,
        help='Output SLURM parameter input config file',
    )
    parser.add_argument(
        '--pdir',
        type=str,
        action='store',
        default=None,
        help='Output directory for job scripts - Path where intermediate and final job files will be written to.',
        metavar='',
    )
    # SLURM specific input parameters - (copies of options in config file)
    slurm = parser.add_argument_group('SLURM arguments')
    slurm.add_argument(
        '--slurm_config',
        type=str,
        action='store',
        default=None,
        help='config.ini file with SLURM arguments',
        metavar='',
    )
    slurm.add_argument(
        '--slurm_jobtype',
        type=str,
        choices=['scsn', 'mcsn', 'mcmn', 'custom'],
        help='Type of SLURM job - scsn=SingleCore-SingleNode, mcsn=MultiCore-SingleNode, mcmn=MultiCore-MultiNode, custom=FillsWithAllProvidedInfo',
    )
    slurm.add_argument(
        '--time',
        type=str,
        action='store',
        default='7:00:00',
        help='Per-job runtime - Input format: day-hr:min:sec or 0-0:00:00',
        metavar='',
    )
    slurm.add_argument(
        '--modules',
        type=str,
        action='store',
        default=None,
        help='List of modules to load (i.e. module load gcc;module load OpenMPI OR module load Anaconda;source activate env)',
        metavar='',
    )
    slurm.add_argument(
        '--nodes',
        type=str,
        action='store',
        default=None,
        help='Number of nodes per-job.',
        metavar='',
    )
    slurm.add_argument(
        '--ntasks',
        type=str,
        action='store',
        default=None,
        help='Number of tasks to run per-job.',
        metavar='',
    )
    slurm.add_argument(
        '--cpus_per_task',
        type=str,
        action='store',
        default=None,
        help='Number of cpus to dedicate to each task.',
        metavar='',
    )
    slurm.add_argument(
        '--memory',
        type=str,
        action='store',
        default=None,
        help='Amount of memory per-job - Provided as: value[K|M|G|T]',
        metavar='',
    )
    slurm.add_argument(
        '--tasks_per_node',
        type=str,
        action='store',
        default=None,
        help='Number of tasks to run per-node',
        metavar='',
    )
    # -- Optional SLURM inputs --
    slurm.add_argument(
        '--account',
        type=str,
        action='store',
        default=None,
        help='Account number for SU payments',
        metavar='',
    )
    slurm.add_argument(
        '--email',
        type=str,
        action='store',
        default=None,
        help='Email to send job notifications to. Automatically sets --mail-type=ALL',
        metavar='',
    )
    slurm.add_argument(
        '--tmp',
        type=int,
        action='store',
        default=None,
        help='Size of $TMPDIR to create - will output intermediates to $TMPDIR rather than output dir [minSize=10240]',
        metavar='',
    )
    # -- SLURM Header direct input --
    slurm.add_argument(
        '--custom_header',
        type=int,
        action='store',
        default=None,
        help='Proivde file with SLURM job header ',
        metavar='',
    )
    # BASH specific input parameters - (copies of options in config file)
    bash = parser.add_argument_group('BASH arguments')
    bash.add_argument(
        '--bash_config',
        type=str,
        action='store',
        default=None,
        help='config.ini file with BASH arguments',
        metavar='',
    )
    bash.add_argument(
        '--threads',
        type=str,
        action='store',
        default=None,
        help='Number of threads/cpus per-job',
        metavar='',
    )
    args = parser.parse_args()
    # --- Raw General Input Variables ---
    INPUT_RAW = args.input
    OUTPUT_RAW = args.output
    PIPELINE = args.pipeline
    SCHEDULER = args.scheduler
    MAKE_SLURM_CONFIG = args.slurm_config_template
    MAKE_BASH_CONFIG = args.bash_config_template
    # --- print help if no arguments are provided ---
    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    else:
        pass
    # --- Ensure input is a valid path ---
    try:
        OUTPUT_PATH = Path(OUTPUT_RAW)
    except TypeError:
        print("ERROR: Invalid output location")
        exit()
    # --- Output SLURM/BASH input parameter config file ---
    if MAKE_SLURM_CONFIG:
        output_file = OUTPUT_PATH / 'SLURM_config.ini'
        slurm_config_template(output_file)
        exit()
    if MAKE_BASH_CONFIG:
        output_file = OUTPUT_PATH / 'BASH_config.ini'
        bash_config_template(output_file)
        exit()
    # --- Ensure valid input file ---
    try:
        INPUT_PATH = Path(INPUT_RAW)
    except TypeError:
        print("ERROR: Input not valid")
        exit()
    # --- Initiate log file ---
    LOG_LEVEL = convert_log_level(args.log_level)
    logger = set_logger(args.log, OUTPUT_PATH, LOG_LEVEL)
    # --- Load input file into dataframe ---
    df = get_input_df(INPUT_PATH, logger)
    # --- Log basic info about input ---
    numRefs = len(df["ReferenceID"].unique())
    numSamples = len(df["QueryID"].unique())
    logger.info("===============================")
    logger.info("*** Input Info ***")
    logger.info(f"Scheduler: {SCHEDULER}")
    logger.info(f"Pipeline: {PIPELINE}")
    logger.info(f"Total Runs: {len(df)}")
    logger.info(f"Total Samples: {numSamples}")
    logger.info(f"Total References:{numRefs}")
    logger.info("===============================")
    # --- Run by JobType ---
    if SCHEDULER == 'SLURM':
        # --- Load SLURM input variables ---
        if args.slurm_config:
            # Load config file
            config = configparser.ConfigParser(allow_no_value=True)
            config.read(args.slurm_config)
            SLURM_OUTDIR = config['SLURM INPUT']['project_directory']
            SLURM_TIME = config['SLURM INPUT']['time']
            SLURM_NODES = config['SLURM INPUT']['nodes']
            SLURM_TASKS_PER_NODE = config['SLURM INPUT']['tasks_per_node']
            SLURM_NTASKS = config['SLURM INPUT']['ntasks']
            SLURM_CPUS = config['SLURM INPUT']['cpus_per_task']
            SLURM_MEMORY = config['SLURM INPUT']['memory']
            SLURM_ACCOUNT = config['SLURM INPUT']['account']
            SLURM_EMAIL = config['SLURM INPUT']['email']
            SLURM_MODULES = config['SLURM INPUT']['modules']
            SLURM_JOBTYPE = config['SLURM INPUT']['job_type']
            try:
                SLURM_TMPDIR = int(config['SLURM INPUT']['tmp'])
            except ValueError:
                if not config['SLURM INPUT']['tmp']:
                    SLURM_TMPDIR = None
                    pass
                else:
                    logger.info(f"Invalid value for tmp - must be integer larger than 10240")
                    exit(1)
            try:
                if not SLURM_MODULES:
                    raise NoModulesFound
                else:
                    pass
            except NoModulesFound:
                logger.error(f"No modules provided, please provide and rerun")
                exit()
        else:
            # Parse SLURM input variables
            SLURM_JOBTYPE = args.slurm_jobtype
            SLURM_TIME = args.time
            SLURM_NODES = args.nodes
            SLURM_TASKS_PER_NODE = args.tasks_per_node
            SLURM_NTASKS = args.ntasks
            SLURM_CPUS = args.cpus_per_task
            SLURM_MEMORY = args.memory
            SLURM_ACCOUNT = args.account
            SLURM_EMAIL = args.email
            SLURM_TMPDIR = args.tmp
            SLURM_OUTDIR = args.pdir
            try:
                if not args.modules:
                    raise NoModulesFound
                else:
                    SLURM_MODULES = args.modules
            except NoModulesFound:
                logger.error(f"No modules provided, please provide and rerun")
                exit()
        # --- Ensure inputs are valid ---
        result = validate_slurm_input(
            logger,
            SLURM_OUTDIR,
            SLURM_JOBTYPE,
            SLURM_TIME,
            SLURM_NODES,
            SLURM_TASKS_PER_NODE,
            SLURM_NTASKS,
            SLURM_CPUS,
            SLURM_MEMORY,
            SLURM_ACCOUNT,
            SLURM_EMAIL,
            SLURM_TMPDIR,
            SLURM_MODULES,
        )
        if not result:  ## Print message and exit if input invalid
            print("INPUT INVALID")
            exit()
        # --- Generate project directory dir script ---
        project_dir_setup_script(SLURM_OUTDIR, OUTPUT_PATH, logger)
        # --- Project Output Directories ---
        SLURM_OUTDIR = Path(SLURM_OUTDIR)
        SLURM_BAM_OUTDIR = SLURM_OUTDIR / "BAM"
        SLURM_ALL_IN_ONE_RUN = OUTPUT_PATH / "run_all_jobs.sh"
        # --- Generate SLURM JobScripts ---
        logger.debug(f"Start SLURM job file creation")
        all_sample_output_files = list()
        for row in df.itertuples(index=False):
            output_filename, header_formatted, command = generate_SLRUM_JobFiles(
                row,
                OUTPUT_PATH,
                PIPELINE,
                SLURM_OUTDIR,
                SLURM_BAM_OUTDIR,
                SLURM_JOBTYPE,
                SLURM_TIME,
                SLURM_NODES,
                SLURM_TASKS_PER_NODE,
                SLURM_NTASKS,
                SLURM_CPUS,
                SLURM_MEMORY,
                SLURM_ACCOUNT,
                SLURM_EMAIL,
                SLURM_TMPDIR,
                SLURM_MODULES,
                logger,
            )
            # -- Log filename --
            logger.info(f"File written: {output_filename}")
            # -- Write output file --
            with open(output_filename, 'w') as oh:
                oh.write(header_formatted)
                oh.write("\n")
                oh.write("\n")
                oh.write(command)
            all_sample_output_files.append(output_filename)
        # --- Generate File To Run All Jobs ---
        logger.info("-------------")
        logger.info(f"Batch job script written: {SLURM_ALL_IN_ONE_RUN}")
        logger.info("-------------")
        with open(SLURM_ALL_IN_ONE_RUN, 'w') as oh:
            oh.write("dir='provide/path/to/jobfiles'\n")
            for f in all_sample_output_files:
                oh.write(f"sbatch $dir/{f.name} &\n")
                continue
    # --- Generate bash JobScripts ---
    elif SCHEDULER == 'BASH':
        logger.debug(f"Start BASH job file creation")
        # --- Load SLURM input variables ---
        if args.bash_config:
            # Load config file
            config = configparser.ConfigParser(allow_no_value=True)
            config.read(args.bash_config)
            BASH_OUTDIR = config['BASH INPUT']['project_directory']
            BASH_CPU = config['BASH INPUT']['threads']
        else:
            # Parse SLURM input variables
            BASH_OUTDIR = args.pdir
            BASH_CPU = args.threads
        # --- Generate project directory dir script ---
        project_dir_setup_script(BASH_OUTDIR, OUTPUT_PATH, logger)
        # --- Project Output Directories ---
        BASH_OUTDIR = Path(BASH_OUTDIR)
        BASH_BAM_OUTDIR = BASH_OUTDIR / "BAM"
        BASH_ALL_IN_ONE_RUN = OUTPUT_PATH / "run_all_jobs.sh"
        # --- Generate BASH JobScripts ---
        all_sample_output_files = list()
        for row in df.itertuples(index=False):
            output_filename, command = generate_BASH_JobFiles(
                row,
                OUTPUT_PATH,
                BASH_OUTDIR,
                BASH_BAM_OUTDIR,
                BASH_CPU,
                PIPELINE,
                logger,
            )
            # -- Write output file --
            with open(output_filename, 'w') as oh:
                oh.write(command)
            all_sample_output_files.append(output_filename)
            # -- Log filename --
            logger.info(f"File written: {output_filename}")
        # --- Generate File To Run All Jobs ---
        logger.info("-------------")
        logger.info(f"Batch job script written: {BASH_ALL_IN_ONE_RUN}")
        logger.info("-------------")
        with open(BASH_ALL_IN_ONE_RUN, 'w') as oh:
            oh.write("dir='provide/path/to/jobfiles'\n")
            for f in all_sample_output_files:
                oh.write(f"./$dir/{f.name} &\n")
                continue
        pass
    return


if __name__ == "__main__":
    main()
