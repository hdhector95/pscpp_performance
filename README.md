# Protein Side Chain Packing Algorithms Performance Comparison
This repository contains the code for the performance comparison of protein side chain packing algorithms. The algorithms used are:
* [FASPR](https://zhanggroup.org/FASPR/)
* [SCWRL4](http://dunbrack.fccc.edu/scwrl4/license/index.html)
* [OSCAR-star](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3187653/)
* [OPUS-Rota4](https://github.com/thuxugang/opus_rota4)

A specific parameter was employed to gauge the accuracy of structural predictions. This assessment involved measuring the percentage of correctly predicted side-chain dihedral angles (ranging from x1 to x1+2) while adhering to a tolerance criterion of 20 degrees.
```
pscp_performance/
    datasets/
        # Contains the datasets used for the comparison
        # Native structures for all datasets
        # Predicted structures for all datasets
        listas/
            # Contains the lists of pdb files used for the comparison
    lib/
        # Contains the libraries used for the comparison
    results/
        # Contains the results of the comparison agrupated by dataset
    config.ini # Contains the parameters for the algorithms
    main.py # Main script for the comparison
```
## Usage

### Dependency
```
Python 3.7
requeriment.txt
```

### RUN Performance Comparison
1. The `config.ini` file contains the parameters for the algorithms. The `pdb_origin_file` and `pdb_processed_files` parameters are the lists of pdb files to be processed and the list of pdb files already processed, respectively. The lists are located in the `datasets/listas` folder.
```
config.ini
    [configuracion]
    correct_angle = 20 # Angle threshold for correct prediction
    print_pdb_console = false # Print results in console
    print_pdb_excel = true # Print results in excel
    pdb_origin_file = pdb_list_origin_XRAY # List of pdb files to be processed
    pdb_processed_files = pdb_list_processed_XRAY # List of pdb files processed
```
 
2. Use `main.py` to run the performance comparison. The script will run the comparison algorithm and save the results in the `results` folder. 

## Results
#### The performance of different side-chain modeling methods on native backbone test sets measured by all residues
#### NMR(21)
|   Method   | χ1(%)  | χ1+2(%)  |
|:----------:|:------:|:------:|
|   FASPR    | 	56.92 | 	34.83 |
|   SCWRL4   | 	54.70 | 	32.88 |
| OSCAR-star | 	58.16 | 	34.60 |
| OPUS-Rota4 | 	62.04 | 	35.47 |
#### XRAY(21)
|   Method   | χ1(%)  |  χ1+2(%) |
|:----------:|:------:|:------:|
|   FASPR    | 	77.30 | 	61.92 |
|   SCWRL4   | 	77.44 | 	61.30 |
| OSCAR-star | 	78.83 | 	61.86 |
| OPUS-Rota4 | 	83.22 | 	68.13 |



## Availability 
This project is freely available for academic usage only.
