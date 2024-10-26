
# Parse_3_df.py

This is a Python script to parse raw ASD data files.

## Usage

To run `Parse_3_df.py`, your command should follow this format (the contents in `[]` are optional):

```bash
python Parse_3_df.py [-o OUTPUT_PATH]
```

### Optional Arguments
- `-o , --output_path`: Path to the job folder. If not specified, the job folder will be stored in the default output folder `./Test_Output/`.

### Example Command
```bash
python SupTable5_script.py -o /output_path/
```

If you run the above command, the output should be found in your specified `output_path` or (`./Test_Output` by default). The output will include the following folder: `/output_path/Run_NetFlow/`.

To get an idea of what the output files should look like, please see the example output files in `/fs/cbsuhy01/storage/NetFlow_ASD/Final_ASD_github/Example_Output/Run_NetFlow`. The main output file (coding count table) will be stored in the path `./Run_NetFlow/new_final_combined_RunNF.txt`.

## Output Files

#### `{output_path}/Run_NetFlow/`
This folder contains all the intermediate files generated. The main, parsed data file is:

#### `new_final_combined_RunNF.txt`
This is a tab-separated file containing the count information for each coding variant type. The first line is a header. This file contains the raw mutation data, annotated with VEP and parsed to merge all three datasets together. The datasets have been parsed and filtered with CADD_PHRED >= 25, MPC >= 2. REVEL scores are annotated.
