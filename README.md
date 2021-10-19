# Short Read Mapping Job Generator (srmjg)
Short Read Mapping Job Generator creates job files for short read mapping pipelines and aims to standardize the commands used.

# Installation
To install srmjg, clone the git respository locally and run the pip installation

    $ git clone https://github.com/harris-2374/srmjg.git
    $ cd srmjg
    $ pip install .


# Dependencies
Pandas is the only required dependency, but if you plan to use .xlsx files you will also need openpyxl. Below are commands for conda enviornments and python virtual environments.

    $ conda install python=3.9.5 pandas openpyxl 
    - or - 
    $ pip install -r requirements.txt

# Input file format
srmjg takes in a tab or comma delimited file with five required column headers. Each row represents an individual mapping event, so this can be used to make mass job files with a combination of reference-query pairs. 

    SampleID | SampleSRR | SampleR1    | SampleR2	 | ReferenceID | ReferencePath
    -------- | --------- | ----------  | ----------  | ----------- | ----------------
    SampleA  | SRR123456 | sA_R1.fasta | sA_R2.fasta | Reference1  | Reference1.fasta
    SampleB  | SRR789101 | sB_R1.fasta | sB_R2.fasta | Reference1  | Reference1.fasta
    SampleA  | SRR123456 | sA_R1.fasta | sA_R2.fasta | Reference2  | Reference2.fasta
    SampleB  | SRR789101 | sB_R1.fasta | sB_R2.fasta | Reference2  | Reference2.fasta

# Usage
srmjg currently supports two different job script types, SLURM job files and bash files for local runs. The SLURM job scripts are based on Texas A&M High Performance Research Computing's SLURM scheduler on their Grace cluster. Visit their [wiki](https://hprc.tamu.edu/wiki/Grace:Batch) for more information on the scheduler and examples of job script types. 

There are two ways to create jobs, through a config file or by command line arguments. Config files are reccomended as they are easier to set up and are better for tracking and reproducibility purposes.
    
    usage: srmjg [-h] [-i] [-o] [--pipeline {bwa-mem2}] [--scheduler {SLURM,BASH}] [--log] [--log_level {NOTSET,DEBUG,INFO,WARNING,ERROR,CRITICAL}] [--slurm_config_template] [--bash_config_template] [--pdir] [--slurm_jobtype {scsn,mcsn,mcmn,custom}] [--slurm_config] [--time] [--modules] [--nodes] [--ntasks] [--cpus_per_task] [--memory] [--tasks_per_node] [--account] [--email] [--tmp] [--bash_config] [--threads]

    optional arguments:
    -h, --help            show this help message and exit
    -i , --input          .xlsx, .tsv, or .csv file - Required column headers = [SampleID | SampleSRR | SampleR1 | SampleR2 | ReferenceID | ReferencePath]
    -o , --output         Output location - default: cwd
    --pipeline {bwa-mem2}
                            Type of job file to create
    --scheduler {SLURM,BASH}
                            Type of job file to create
    --log                 Name of log file
    --log_level {NOTSET,DEBUG,INFO,WARNING,ERROR,CRITICAL}
                            Log level
    --slurm_config_template
                            Output SLURM parameter input config file
    --bash_config_template
                            Output SLURM parameter input config file
    --pdir                Output directory for job scripts - Path where intermediate and final job files will be written to.

    SLURM arguments:
    --slurm_config        config.ini file with SLURM arguments
    --slurm_jobtype {scsn,mcsn,mcmn,custom}
                            Type of SLURM job - scsn=SingleCore-SingleNode, mcsn=MultiCore-SingleNode, mcmn=MultiCore-MultiNode, custom=FillsWithAllProvidedInfo
    --time                Per-job runtime - Input format: day-hr:min:sec or 0-0:00:00
    --modules             List of modules to load (i.e. module load gcc;module load OpenMPI -OR- module load Anaconda;source activate env)
    --nodes               Number of nodes per-job.
    --ntasks              Number of tasks to run per-job.
    --cpus_per_task       Number of cpus to dedicate to each task.
    --memory              Amount of memory per-job - Provided as: value[K|M|G|T]
    --tasks_per_node      Number of tasks to run per-node
    --account             Account number for SU payments
    --email               Email to send job notifications to. Automatically sets --mail-type=ALL
    --tmp                 Size of $TMPDIR to create - will also output intermediates to $TMPDIR rather than output dir

    BASH arguments:
    --bash_config         config.ini file with BASH arguments
    --threads             Number of threads/cpus per-job

# Example Commands
    SLURM command with configuration file:
        $ srmjg -i ./tests/input.xlsx -o ./tests/output --scheduler BASH --bash_config ./tests/BASH_config.ini
    
    SLURM command without configuration file:
        $ srmjg -i test_txt.txt -o path/to/output --pdir project/dir/path --scheduler SLURM --slurm_jobtype scsn --slurm_config SLURM_config.ini --tmp 10240
    
    BASH command with configuration file:
        $ srmjg -i ./tests/input.xlsx -o ./tests/output --scheduler BASH --bash_config ./tests/BASH_config.ini

    BASH command without configuration file:
        $ srmjg -i ./tests/input.xlsx -o ./tests/output --scheduler BASH --pdir project/dir/path --threads 5


# Supported SLURM arguments
    1. #SBATCH --time=<time>
        - Input style = day-hr:min:sec

    2. #SBATCH --partition=<queue>
        - TAMU HPRC options: short, medium, long, xlong

    3. #SBATCH --ntasks=<ntasks>

    4. #SBATCH --cpus-per-task=<cpus-per-task>

    5. #SBATCH --mem=<memory>
        - Provided as: value[K|M|G|T]

    6. #SBATCH --account=<accountnumber>

    7. #SBATCH --ntasks-per-node=<tasks-per-node>

    8. ##SBATCH --mail-user=<email>
        - Auto set this value when email is provided - ##SBATCH --mail-type=ALL

    9. #SBATCH --nodes=<[min[-max]]>

    10. #SBATCH --tmp=<tmp>
        - Provide in MB
        - If provided, writes all but last output file to $TMPDIR

# Example SLURM config file (i.e. slurm_config.ini)
    [SLURM INPUT]
    project_directory = project/dir/path
    job_type = snsc
    time = 7:00:00
    queue = medium
    nodes = 1
    tasks_per_node = 1
    ntasks = 1
    cpus_per_task = 1
    memory = 10G
    account = 
    email = 
    tmp = 
    modules = 

# Example Bash config file (i.e. bash_config.ini)
    [BASH INPUT]
    project_directory = project/dir/path
    threads = 1


# SLURM Limitations
Currently, srmjg does not support GPU jobs and only provides a subset of commonly used SLURM arguments. Open an issue if you would like other features to be added in future versions.




