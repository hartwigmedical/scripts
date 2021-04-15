# UMI-undup README

UMI-undup is a script for unmarking reads as duplicate based on their UMI after sambamba markdup has marked duplicates 
and UMI-tools group has grouped reads by UMI.
 
## Installation
If you want to run the code, please generate a local Python 3 venv and install the requirements:

```bash
$ python3 -m venv [path/to/new/virtual/environment, for example: ./umi]
$ source [path/to/new/venv, for example: ./umi/bin/activate]
(ensembl) $ pip install -r requirements.txt
```

Or run the scripts create_umi_venv and run_umi_undup

## Usage
Remember to source the virtualenv before running `main.py`.

####General usage
```(umi) $ python main.py input.txt```

####Arguments
* `input_bam`: (Required) Path to input bam file.
* `output_bam`: (Required) Path to output bam file.