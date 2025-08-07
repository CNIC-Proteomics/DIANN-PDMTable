## DIANN-PDMTable
Variaton of the [PDM Table Maker Module](https://github.com/CNIC-Proteomics/PTM-compass/blob/main/MODULES.md#1-pdmtablemaker) to create the PDM Table from DIA-NN output. Its purpose is to adapt the DIA-NN output to the [iSanXoT](https://github.com/CNIC-Proteomics/iSanXoT) analysis and then the [Report Analysis](https://github.com/CNIC-Proteomics/nf-PTM-Analyzer).

## Usage
The _PDMTableMaker_DIA.py_ script does all the work, but it is advisable to review the _params_PDMTableMaker_DIA.ini_ configuration file to avoid errors. The main parameters to configure are:

 * `Sequence_column_name`: name of the column in the input data table that contains the (modified) sequence of the peptides. Typically, in DIA-NN reports, this column corresponds to _Modified.Sequence_.
 * `MasterProtein_column_name`: name of the column that contains the protein identifier. It may be _Protein.Group_ if the DIA-NN report is used directly, but it's advisable to execute [Protein Assigner](https://github.com/CNIC-Proteomics/PTM-compass/blob/main/src/tools/ProteinAssigner.py) module over the DIA-NN report and before this module, therefore each peptidoform is assigned to a single protein.
 * `Outfile_suffix`: suffix to be added to the output PDMTable file.

**Input**:
 - Tabular file containing at least one column with the peptidoforms sequences, one column with the protein identifiers and the intensity (or abundance) of each peptidoform per sample. It could be the `report.pr_matrix.tsv` file obtained directly from DIA-NN, or the output of Protein Assigner module. Should be passed as parameter to the script using the `-i` or `--infile` flag.
 - Fasta file with the protein sequences of the organism analyzed. Decoy sequences should **not** be included. Should be passed as parameter to the script using the `-f` or `--fastafile` flag.
 - Config file with the parameters described above. Should be passed as parameter to the script using the `-c` or `--config` flag.

**Output**:
 - `*_plus.tsv`: same tabular file as the one used as input, but with three additional columns: _d_ (deltamass of the peptidoform), _g_ (group of the peptidoform) and _pgm_ (identifier of the peptidoform with the group and the position of the modification).
 - `*_PDMTable_DIAversion.txt`: PDM Table with all the necessary information per _pgm_.

**Execution**:
> ⚠️ _Adapt the paths to the script and the files to match your case!_
```
python PDMTableMaker_DIA.py -i DIA-NN_report.tsv -f organism_fasta.fasta -c params_PDMTableMaker_DIA.ini
```

<pre style="white-space: pre-wrap;">
 python PDMTableMaker_DIA.py -i DIA-NN_report.tsv -f organism_fasta.fasta -c params_PDMTableMaker_DIA.ini
</pre>
