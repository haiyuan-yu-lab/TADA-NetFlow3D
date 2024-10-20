
# SupTable5_script.py

This is a Python script to generate a coding count analysis table.

## Usage

To run `SupTable5_script.py`, your command should follow this format (the contents in `[]` are optional):

```bash
python SupTable5_script.py -I <JOB_NAME> [-s SIB_FILE] [-R RAW_FILE_PATH] [-r ROLES_FILE] [-g SAMPLE_INFORMATION] [-m MUTABILITY_FILE] [-o OUTPUT_PATH]
```

### Required Argument
- `-I <job_name>`: Replace `<job_name>` with a name of your choice for the current job. The output files will be stored in a folder with this name.

### Optional Arguments
- `-s , --sib_file`: Path to the sibling information file. Used to exclude siblings with ASD risk. Default path is: `"/fs/cbsuhy01/storage/NetFlow_ASD/Week8_429/individuals_registration.csv"`.

- `-R , --raw_file_path`: Path to UCSF raw mutation files folder. Default path is: `"/fs/cbsuhy01/storage/NetFlow_ASD/Week8_429/UCSF_raw/"`.

- `-r , --roles_file`: Path to roles file. Used to exclude multiplex families and siblings with ASD risk. Default path is: `"/fs/cbsuhy01/storage/NetFlow_ASD/Week8/SSC+SPARK+KOR-STAR_parsed.txt"`. 

- `-g , --sample_information`: Path to sample information file. Contains the information on sex and phenotype (whether case or sibling). The default path is `"/fs/cbsuhy01/storage/NetFlow_ASD/Week8/SSC+SPARK+KOR-STAR_parsed.txt"`.

- `-m , --Mutability_file`: Path to the Mutability file. The default path is `"/fs/cbsuhy01/storage/NetFlow_ASD/Week8/New_Analysis/coding_analysis/41588_2022_1104_MOESM3_ESM.xlsx"`.

- `-o , --output_path`: Path to the job folder. If not specified, the job folder will be stored in the default output folder `./Test_Output/`.

### Example Command
```bash
python SupTable5_script.py -I test_script
```

If you run the above command, the output should be found in your specified `output_path + /job_name/` or (`./Test_Output + /job_name/` by default). The output will include the following 4 folders: `/DNV/`, `/Transmitted/`, `/Untransmitted/`, `/Summary/`.

Note: If using default output_path, the folder will be created in your current directory.

To get an idea of what the output files should look like, please see the example output files in `/fs/cbsuhy01/storage/NetFlow_ASD/Week8/Final_ASD_github/test_902`. The main output file (coding count table) will be stored in the path `./Summary/Total_CodingAnalysis_Table.txt`.

## Output Files

#### `{job_name}/{Summary}/`
This folder contains the summary files. There are three files in this folder. Two of them are count tables for De Novo and Transmitted variants, respectively. The overall coding count table is:

#### `Total_CodingAnalysis_Table.txt`
This is a tab-separated file containing the count information for each coding variant type. The first line is a header, and the table includes 21 columns:

1. `gene`: Gene name
2. `gene_gencodeV33`: Gene name in gencodeV33 annotation
3. `gene_id`: ENSEMBL gene ID
4. `chromosome`: Chromosome of the gene
5. `pLI`: Probability of being loss-of-function intolerant (constraint metric)
6. `LOEUF`: Loss-of-function observed/expected upper bound fraction (constraint metric)
7. `mut.ptv`: De Novo protein truncating (LoF) mutation rate
8. `mut.misb`: De Novo Missense B (misB; MPC ≥ 2) mutation rate
9. `mut.misa`: De Novo Missense A (misA; 1 ≤ MPC < 2) mutation rate
10. `DNV.ptv`: Number of De Novo PTVs in affected offspring
11. `DNV.ptv.sib`: Number of De Novo PTVs in unaffected offspring
12. `DNV.misb`: Number of De Novo misBs in affected offspring
13. `DNV.misb.sib`: Number of De Novo misBs in unaffected offspring
14. `DNV.misa`: Number of De Novo misAs in affected offspring
15. `DNV.misa.sib`: Number of De Novo misAs in unaffected offspring
16. `Transmitted.ptv`: Number of transmitted PTVs to affected offspring
17. `Transmitted.ptv.sib`: Number of transmitted PTVs to unaffected offspring
18. `Transmitted.misb`: Number of transmitted misBs to affected offspring
19. `Transmitted.misb.sib`: Number of transmitted misBs to unaffected offspring
20. `Transmitted.misa`: Number of transmitted misAs to affected offspring
21. `Transmitted.misa.sib`: Number of transmitted misAs to unaffected offspring

#### `/{job_name}/{DNV}/`
This folder contains all the intermediate files produced for De Novo variant analysis.

#### `/{job_name}/{Transmitted}/`
This folder contains all the intermediate files produced for Transmitted variant analysis.

#### `/{job_name}/{Untransmitted}/`
This folder contains all the intermediate files produced for Untransmitted variant analysis.
