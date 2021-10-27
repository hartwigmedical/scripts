# SearchEnsembl README

SearchEnsembl is a script for searching Ensembl Gene Ids and canonical names for human genes.
 
## Installation
(There are no non-default requirements yet, so creating a venv is not yet necessary.)

If you want to run the code, please generate a local Python 3 venv and install the requirements:

```bash
$ python3 -m venv [path/to/new/virtual/environment, for example: ./ensembl]
$ source [path/to/new/venv, for example: ./ensembl/bin/activate]
(ensembl) $ pip install -r requirements.txt
```

## Usage
(There are no non-default requirements yet, so creating a venv is not yet necessary.)
Remember to source the virtualenv before running `main.py`.

####General usage
```(ensembl) $ python3 main.py input.txt```

####Arguments
* `input`: (Required) Path to file with gene names.