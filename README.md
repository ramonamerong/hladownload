# hladownload<br>*A tool for downloading HLA alleles, epitopes and frequencies* 
As HLA alleles, epitopes/eplets and frequencies often need to be manually downloaded and parsed, I developed several scripts which could automatically download and parse these data for my Bachelor's thesis project. In order to make this available to the public I combined these scripts and developed this tool.

With this tool you will be able to automatically download HLA allele sequencies/alignments from the [IMGT/HLA database's github](https://github.com/ANHIG/IMGTHLA), epitopes/eplets from [the eplet registry](https://www.epregistry.com.br) and HLA frequencies from [the HLA frequency database](http://www.allelefrequencies.net). 

This can either be done through the [command line program](#command-line-program) or in python through [the separate modules](#python-modules) (`hla`, `eplet` and `frequency`). When doing the latter additional functionalties will be avaiable. The `hla` module, for example, will give you access to an `HLA` object which allows you to programatically retrieve, parse and slice HLA alleles by their codon positions.

## Installation
There are several ways to use and install this program.
* ```python -m hladownload``` When this or an OS equivelant command is run in a folder with the files of [the github repository of this project](https://github.com/ramonamerong/hladownload), the command line program can be used without any further installation. However, it is required that the python dependencies listed in `setup.py` are installed.
* ```python setup.py install``` This command must also be run in a folder with the repository files. The tool will then be installed under the current used python installation and will be available on the command line or in python through the name '`hladownload`'.
* ```pip install hladownload``` This command will download the tool through [pip](https://pypi.org/project/hladownload).

For uninstallation you can do:
```pip uninstall hladownload```
You could also get a prompt asking if 'possible manually installed files' should be deleted. These cache files are created by the program in the folder `hladownload/alignmentCache`  and can safely be deleted.

## Command line program
The command line program has several arguments which will be explained and highlighted in an example here.

### Output flags
These flags will enable you to specify the desired output: HLA allele sequences, eplets or frequencies. The output can further be restricted by specifying alleles for which only the sequences, eplets and/or frequencies should be reported (see `-a` or `-A` under [Other arguments](#other-arguments)). Not all HLA loci are available for every type of output:
* `-H, --HLA` Flag to indicate that HLA alleles should be saved as multiple sequence alignments in fasta format per locus. These allele frequencies are retrieved from [the IMGT/HLA database github](https://github.com/ANHIG/IMGTHLA).<br>*Supported loci: A, B, C, DRA, DRB1, DRB3, DRB4, DRB5, DRB345, DQA1, DQA2, DQB1, DPA1, DPB1*
* `-E, --eplet` Flag to indicate that eplets should retrieved from [the eplet registry database](https://www.epregistry.com.br) and outputted them into three csv files containing a summary (`summary_table-eplets.csv`), allele-eplet associations (`allele_table-eplets.csv`) and residues per eplet (`residue_table-eplets.csv`).<br>*Supported loci: A, B, C, DRB1, DRB3, DRB4, DRB5, DRB345, DQA1, DQB1, DPA1, DPB1*
* `-F, --frequency` Flag to indicate that HLA allele frequencies should be retrieved from [the allele frequency database](http://www.allelefrequencies.net) and outputted into two table files: One with the average allele frequency per region and the other globally (averaged over all (specified) regions). The total sample size is also calculated per row.<br>*Supported loci: A, B, C, DRB1, DQA1, DQB1, DPA1, DPB1*

### `-H` additional arguments
  * `-t, --slice_type {no_gaps,gaps,amino}` Can only be used if `-H` is given as a flag. Indicate whether the sliced codons should include gaps, no gaps or should be translated to amino acids. These gaps are only referring to those that are introduced into the reference sequence due to the alignment. Gaps are included by default.
  * `-s, --slice_range SLICE_LOCI SLICE RANGES` Can only be used if `-H` is given as a flag. Indicate the codons of the HLA alleles per locus that should be preserved in the outputted alignment files. The arguments should consist of a comma-separated list of loci, followed by a combination of ranges and comma separated numbers for the colons to preserve: e.g. `-s A,B,C 1-5,10` means codons 1 up to (but not including) 5 and 10 of loci A, B and C should be preserved. You can also specify different codons for different loci by repeating the argument, e.g. `-s A,B 1-5 -s DRB1,DRB345 5-8 etc.`. By default all codons of the aligned alleles of a locus will be preserved if that locus is not specified.

### `-E` additional arguments
* `-l, --all_eplets` Flag to indicate that all eplets should be included, as any eplets which are not associated with any (specified) alleles are removed from the output by default.

### `-F` additional arguments
* `-r, --region REGION1 REGION2 etc...` Indicate the regions separated by spaces for which the HLA frequencies should be retrieved. A region consisting of multiple words should be enclosed in quotes. If no regions are specified all of the regions will be reported. If it is specified the global frequencies will only be based on the average of those regions. <br>*Usable regions are: Australia, Europe, North Africa, North America, North-East Asia, Oceania, South and Central America, South Asia, South-East Asia, Sub-Saharan Africa, Western Asia*
* `-hr, --higher_resolution` Flag to indicate that frequencies of higher resolution alleles of the specified alleles (see '-a' or '-A') should also be reported.
* `-lr, --lower_resolution` Flag to indicate that frequencies of lower resolution alleles of the specified alleles (see '-a' or '-A') should also be reported.

### Other arguments
* `-h, --help` Show an help message and exit.
* `-o, --output_directory OUTPUT_DIRECTORY` Directory to save the output files to. For every output flag (`-H`, `-E` or `-F`) that is enabled, a separate directory (`/alleles`, `/eplets` or `/frequencies`) will be created in the output folder with the output files. The output folder defaults to the current folder.
* `-d, --database_version DATABASE_VERSION` The IMGT/HLA database version github (https://github.com/ANHIG/IMGTHLA) to use to retrieve HLA alleles from. Defaults to the latest branch.
* `-a, --alleles ALLELES` Indicate the HLA alleles for which an alignment, frequencies and/or eplets should be outputted. The alleles can be a combination of ranges and comma separated allele names without spaces: 'A\*01:01:01:01-A\*01:01:01:09,B\*01:01:01:01'. If an allele cannot be found, it will be automatically replaced by its higher/lower resolution variants in the IMGT/HLA database. These alleles will then be reported in the eventual alignment and eplet files. However, for allele frequencies higher/lower resolution alleles will only serve as a replacement if both their allele frequencies can be found and if the flags `->` and/or `-<` have been used. The default is that all alleles of the given IMGT/HLA database version are included (for HLA alignment and eplets) or all alleles for which frequencies can be found (for frequencies).
* `-A, --alleles_file ALLELES_FILE` Indicate in a file the HLA alleles for which an alignment, frequencies and/or eplets should be outputted. The allele names should be separated by new lines in the file. This argument overrides the alleles specified with `-a`.
* `-n, --no_null` Flag to indicate whether HLA null alleles should be excluded.
* `-g, --group_DRB345` Flag to indicate whether HLA alleles of the DRB3, -4 and -5 loci should be grouped under a single locus 'DRB345'.

### Examples
`hladownload -a A*31:01 -o output -H -s A 1-10 -E -F -r Europe Australia -> -n`
Using this command, no alleles corresponding to 'A\*31:01' will be found in the IMGT/HLA alleles. Therefore, higher resolution alleles are returned instead for which only codons 1 up to 9 are preserved in the eventual sequence and only eplets occuring in these alleles will be reported. In addition, allele frequencies of A*\31:01 or of higher resolution alleles will be reported from the region Europe and Australia. Null alleles, however, will all be excluded from the results (otherwise allele A\*31:01:02:03N could have been included). The output files of this example can be found under `example_output` on [the github repository of this project](https://github.com/ramonamerong/hladownload).

## Python modules
When this program has been installed, different submodules can be found under `hladownload`:
```
from hladownload import hla, eplet, frequency
```
See the previous section about [output flags](#output-flags) to see which HLA loci each module supports.

### HLA submodule
The `hla` submodule uses IMGT/HLA per locus allele alignments to allow easy retrieval and slicing of HLA allele sequences/alignments.

#### Retrieving HLA allele alignment files
The `hla.downloadAlignments(outFolder, molecule, loci)` function can be used to download the alignments per locus from [the database's GitHub](https://github.com/ANHIG/IMGTHLA/tree/Latest/msf) for all the loci that are present as strings in a list as part of the third parameter. The second parameter should specify whether to download nucleotide (`'nt', 'nuc', 'nucleotide'`) or amino acid (`'aa', 'prot', 'amino acid'`) alignments. The last parameter specifies which from which database version the files should be downloaded.

Example code:
```
downloadAlignments('nuc_alignments', 'nt', databaseVersion = 'Latest');
```

#### Creating the main HLA alignment object
The class `hla.<HLA>(alignmentInput, format, codonStarts, databaseVersion, noNull)` is mainly used for retrieval and slicing of sequences. In order to create an object from this class either a folder name with alignment files (see [the previous section](#retrieving-hla-allele-alignments), a list of file names or nothing at all must be provided for `alignmentInput`. In the latter case, the class uses the `hla.downloadalignments()` to download the files from the internet and you can therefore again specify the database version as the fourth parameter (by default it used the `Latest` github branch for the alignment files). The second parameter must specify the file format (which is `.msf` by default). Since this class uses Biopython's AlignIO module for multisequence alignment objects, the supported formats can be found [in their documentation](http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec80). The third parameter, `codonStarts`, uses a (default) dictionary with locus names as keys and codon counting ofsets as values. The last argument is a boolean which specifies whether null alleles should be filtered out or not (by default they are left in).

Example code:
```
#Load alignment files from folder
HLAObjectFolder = HLA('nuc_alignments', format = 'msf')

#Download alignment files from the github with database version 3.43.0
HLAObjectOnline = HLA(databaseVersion = '3.43.0')

#Load alignment from a pickle object (see the section further on)
HLAObjectPickle = loadAlignments('pickles')
```

#### Retrieving alleles using the main HLA alignment object
The alignments are subsequently stored as [Biopython MultipleSeqAlignments](http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec79) in the property `<HLA>.alignments` in a dictionary under their respective locus key. Specifically, the dictionary values are custom `<LocusAlignment>` objects which function as a wrapper for the Biopython alignment. The Biopython alignment is located in the `<LocusAlignment>.alignment` property. Therefore, in order to access a Biopython alignment you must do `<HLA>.alignments[locus].alignment` or use the shorthand version `<HLA>.alignments.align(locus).alignment`.

You can also retrieve specific alleles in the following ways:
1. `<HLA>.names(optionalLocus, returnDict)`: Retrieve either all (`optionalLocus = None`) or locus specific allele names as strings (`optionalLocus = locusName`) in either a list (`returnDict = False`) or dictionary (`returnDict = True`).
2. `<HLA>[Name]`: Retrieve a single allele as a [Biopython SeqReq object](http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec17) by providing the name between slicing operators. When the allele is not present a list of alleles with either preferably higher or otherwise lower resolution alleles will be returned instead.
2. `<HLA>[Locus]`: Retrieve the reference allele of a locus as a Biopython SeqReq object. The reference sequence is derived from  the first sequence of a locus alignment when creating the class, which is subsequently stored in `<HLA>.alignments[locus].ref`.
3. `<HLA>[Name1:Name2, Name3, Name4]`: Create a new child HLA object (specifically from the `<HLASlice>` class) with only the alleles with the names separated by commas. These can either be ranges of alleles (`Name1:Name2`) or specific alleles (`Name3, Name4`). The alleles are automatically assigned to their appropiate locus alignment. When providing only a single range, you do not have to use commas (`<HLA>[Name1:Name2]`). Alternatively, you can also provide a specific locus name instead of a specific allele name in order to include the reference allele in the Biopython alignment. However, irrespective of whether you include the locus name or not in your selection, the reference allele can always be found in `<HLA>.alignments[locus].ref`. You cannot use locus names when specifying a range of alleles.

Example code
```
print('Total number of alleles:', len(HLAObject.names()), '\n')              #1.
print('Total number of locus A alleles:', len(HLAObject.names('A')), '\n')   #1.
print(HLAObject['A*01:01:01:01'], '\n')                                      #2.
print(HLAObject['A'], '\n')                                                  #3.
print(HLAObject['A', 'C*01:02:01:14':'C*01:02:01:19'].names(), '\n')         #4.
```

Example output:
```
Total number of alleles: 30160 

Total number of locus A alleles: 6766 

ID: A*01:01:01:01
Name: A*01:01:01:01
Description: A*01:01:01:01
Number of features: 0
/weight=1.0
Seq('ATGGCCGTCATGGCGCCCCGAACCCTCCTCCTGCTACTCTCGGGGGCCCTGGCC...TGA') 

ID: A*01:01:01:01
Name: A*01:01:01:01
Description: A*01:01:01:01
Number of features: 0
/weight=1.0
Seq('ATGGCCGTCATGGCGCCCCGAACCCTCCTCCTGCTACTCTCGGGGGCCCTGGCC...TGA') 

['A*01:01:01:01', 'C*01:02:01:14', 'C*01:02:01:15', 'C*01:02:01:16', 'C*01:02:01:17', 'C*01:02:01:18', 'C*01:02:01:19'] 
```

#### Slicing alleles of the main HLA alignment object
Using the `<LocusAlignment>` objects in `<HLA>.alignments`, you can retrieve a subsequence as a Biopython alignment from the main locus alignment by specifying codon positions., simular to how the command line program works. The first codon's position is derived from the integer stored in the `codonStarts` dictionary when creating a `<HLA>` object. Any subsequent codons are numbered starting from that position. Below are listed several types of aligned subsequences that can be retrieved:

1. `<LocusAlignment>.slice[]`: Retrieve a nucleotide subsequence including gaps. This operation is relatively fast.
1. `<LocusAlignment>.cod[]`: Retrieve a nucleotide subsequence excluding gaps. When an allele contains a lot of insertions or deletions compared to the reference allele the result of this operation may be unpredictable. This operation can be very slow for long subsequences or lots of alleles.
1. `<LocusAlignment>.amino[]`: Retrieve an amino acid subsequence. The same precautions as for `2.` apply here.

The rules for specifying codons are similar to `<HLA>[]`:
1. `<LocusAlignmen>.X[codonPos]`: By specifying a single codon position, you can retrieve the alignment of a single codon or amino acid.
1. `<LocusAlignment>.X[codonPos1, codonPos2:codonPos3, codonPos4:, :]`: You can also specify multiple codons separated by commas which are eventually all appended in that order to form a single subsequence. You can see the start and end positions of the codons in an alignent by doing `<LocusAlignment>.range()`.
 1. '`codonPos1`' returns a single codon at that position.
 1. '`codonPos2:codonPos3`' returns all codons between these two positions. The start (`codonPos2`) is included but the end (`codonPos3`) is not. 
 1. '`codonPos4:`' returns all codons starting from `codonPos4` up till the end of the sequence. Contrarily, when only the end position is given all codons from the start up till that position are returned (noninclusive).
 1. '`:`' returns the whole sequence. In this case, the whole sequence is thus appended after all the previous mentioned codons.

The Biopython alignments that are subsequently returned can be supplied to the function `hla.countAlignment(alignment, charLen)` to count the number of nucleotides (`charLen=1`), codons (`charLen=3`) or amino acids (`charLen=1`) at every position. It then returns a list with a [Python counter object](https://pymotw.com/2/collections/counter.html) for every position. In addition, these alignments can also be saved with `hla.saveAlignment(alignment, fileName, outputFolder, format)` instead of using biopython's `AlignIO` module.

Example dode:
```
locusCAlignment = HLAObject.align('C')

#Print the default codon start positions of all loci
print('Codon starts:', defaultCodonStarts, '\n')

#Print the codon range of this alignment
print('Start and end codon positions of C locus:', locusCAlignment.range(), '\n')

#The difference between slicing with slice[], cod[] and amino[]
print('Codon 3 with gaps:\n', locusCAlignment.slice[3], '\n')
print('Codon 3 without gaps:\n', locusCAlignment.cod[3], '\n')
print('Codon 3 as amino acid:\n', locusCAlignment.amino[3], '\n')

#Slicing with multiple codon positions
print('Multiple codons:\n', locusCAlignment.slice[3, 5:8, 10:, :], '\n')

#Count the kind of codons at position 1
print('Count codons at position 1:\n', countAlignment(locusCAlignment.cod[3], charLen=3))

#Save the alignment
saveAlignment(locusCAlignment.alignment, 'C', 'alignments', 'fasta')
```

Example output:
```
Start and end codon positions of C locus: (-24, 343) 

Codon 3 with gaps:
 Alignment with 6621 rows and 4 columns
C-AC C*01:02:01:01
...
C-AC C*18:14 

Codon 3 without gaps:
 Alignment with 6621 rows and 3 columns
CAC C*01:02:01:01
...
CAC C*18:14 

Codon 3 as amino acid:
 Alignment with 6621 rows and 1 columns
H C*01:02:01:01
...
H C*18:14 

Multiple codons:
 Alignment with 6621 rows and 2905 columns
C-ACATGAAGTATACATCCGTGTCCCGGCCTGGCCGCGGAGAGC...--- C*01:02:01:01
...
C-ACATGAGGTATACCGCCGTGTCCCGGCCCGGCCGCGGAGAGC...--- C*18:14 
```

#### Saving the main HLA alignment object
All of the alleles in the HLA alignment object can be saved with `<HLA>.save(outputFolder, outputFolder, format, sliceDict, sliceType`. The output can either be saved as a python pickle object (`format=pickle`) to preserve the functionality of the object or as multiple sequence alignment files per locus (again see Biopython's documentation for the supported formats). 

The pickle object can be loaded in any script with `loadAlignments(input, dbVersion, noNull)`, as long as that script contains all the classes etc. within within the module (`from hladownload.hla import *`). When a folder is given for `input` instead of a file, it is checked for pickle objects by their version number without dots ('`3.43.0`' becomes '`3430`') and whether the alignment contains null alleles or not (`all`or `nonull`) (as pickle objects are normally saved as '`<dbVersion>-all.pickle`' or '`<dbVersion>-nonull.pickle`'). In case no pickle object can be found in the folder, a new main HLA alignment object is created anew and saved into this folder.

Similarly to the command line program, `sliceDict` and `sliceType` can be used to specify which allele codons to preserve per locus and whether they should include gaps/should be translated to amino acids. However, this only works when saving the object as alignment files per locus. `sliceDict` must be contain the loci as keys and the string describing the codons to preserve as values: `{'A': '1-10,15,20-25', 'B': ...}`. `sliceType` can be `gaps`, `no_gap` or `amino` to include gaps introduced into the reference sequence, exclude gaps or translate codons, respectively.

Example code:
```
#Save the HLA alignment object as a pickle
HLAObject.save('pickles', format = 'pickle')

#Save the alignments of the HLA alignment object as fasta files but slice all locus A alleles and remove any gaps introduced in the reference sequence
HLAObject.save('alignments', format = 'fasta', sliceDict = {'A': '1-10,15,20-25'}, sliceType = 'no_gaps')
```

### Eplet submodule
The `eplet` submodule contains several functions for retrieving and parsing eplets. The input for these functions should be a `.csv` table per locus containing all eplets. This data can be retrieved from the eplet registry database for each of the following loci: [A/B/C](https://www.epregistry.com.br/index/databases/database/ABC/) (ABC), [DRB1/DRB345](https://www.epregistry.com.br/index/databases/database/DRB/) (DR), [DQA1/DQB1](https://www.epregistry.com.br/index/databases/database/DQ/) (DQ) and [DPA1/DPB1](https://www.epregistry.com.br/index/databases/database/DP/) (DP). In between parentheses is something that will be referred to as the '__locus group__' from now on, as eplets are separated based on that. The file names should reflect which locus group it belongs to, for example `ABC-eplets.csv` (for A/B/C) or `DRB-eplets.csv` (for DRB1/DRB345). The tables themselves should at least include the following columns:

eplet | polymorphic_residues | antibody_reactivity
--- | --- | ---
9Y | 8V9Y  | (Empty or 'Yes')

Alternatively, the function `retrieveEplets(output)` (explained below) can be used instead to download all of these files automatically:
1. `retrieveEplets()`: This function is used to parse the HTML pages per loci of the [eplet registry database](https://www.epregistry.com.br/index/databases/database/ABC/) into pandas dataframes per locus group which are then returned in a dictionary.
1. `parseResidues(residues)`: This function parses the string of polymorphic residues as found in the column table and returns a dictionary of their positions as keys and their possible amino acids in a list as values.
1. `readEplets(input, delimiter)`: This function either reads a table for a single locus group from a file or from a pandas dataframe in order to construct a dictionary with the eplet names as keys. As values there is another dictionary indicating whether the eplet has been verified to have antibody reactivity ('`verified`') and which polymorphic residues it encompasses ('`residues`'). These residues are contained in the dictionary that is returned by `parseResidues()`.
1. `mapEpletsToAlleles(eplets, alignmentObject, locusGroup, loci)`: This function requires the eplet dictionary for a single locus group returned by `readEplets()` and an HLA alignment object (which should preferably not contain null alleles). These are then used to determine which eplets are contained in which alleles, but only for those alleles that belong to a locus that is included in `locusGroup` and, if specified, is present in `loci`. A similar dictionary as `eplets` is returned, with the addition of the '`alleles`' key which contains allele name lists mapped by their locus in a dictionary.
4. `getMappedEpletTables(epletMaps, outputFolder, allEplets)`: This function takes a dictionary with as keys the locus groups and as values the different dictionaries generated by `mapEpletsToAlleles()`. It then converts all of the information into three pandas tables. It can be specified with 'allEplets' whether eplets without the alleles associated with them should be removed from these tables:
    1. `eplet_summary.csv`: Contains general information of the eplet on each row.
        1. `eplet`: The eplet name. This is included in all tables in order to identify an eplet.
        1. `locus_group`: The locus group the eplet belongs to. This is defined on [the eplet registry website](https://www.epregistry.com.br/index/databases/database/ABC) and is necessary  in order to distinguish eplets with the same name but from different loci (groups).

        1. `verified`: This shows whether the eplet has been verified to have antibody reactivity.
        1. `residue_number`: The number of residues that are represented by the eplet.
        1. `loci`: Comma separated loci names. These loci contain alleles with the eplet.
        
      eplet | locus_group | loci | verified | residue_number | allele_number
      --- | --- | --- | --- | --- | ---
      9Y | ABC | A,B | True | 2 | 2
    
    1. `eplet_residues.csv`: Contains information on the residue positons and amino acids per eplet. Each eplet and residue is a separate row.
        1. `position`: Amino acid position of the residue.
        1. `residue`: One letter code for the amino acid residue. Sporadically, this can feature multiple amino acids separated by a '`,`'.
    
      eplet | locus_group | position | residue
      --- | --- | --- | ---
      9Y | ABC | 8 | V
      9Y | ABC | 9 | Y
    
    1. `eplet_alleles.csv`: Contains information on which alleles are associated with which eplets. Each row is repesented by an eplet/allele association.
        1. `locus`: The locus the allele belong to.
        1. `allele`: The allele name of the allele that contains the eplet.
    
      eplet | locus_group | locus | allele
      --- | --- | --- | ---
      9Y | ABC | A | A\*01:01
      9Y | ABC | B | B\*01:01
    
1. `saveEplets(outputFolder, alignmentObject, extension, allEplets, loci)`: This function uses all of the functions above to directly save the resulting eplet tables as ',' separated `.csv` files by default. You can optionally restrict for which loci eplets should be reported with 'loci'.

Example code
```
#Retrieve and save the eplets that are associated with the alleles of only the MHCI loci
saveEplets('eplets', HLAObject, allEplets=False, loci=['A', 'B', 'C'])
```

### Frequency submodule
The `frequency` submodule contains three functions which in order retrieve, filter and save HLA allele frequencies from different global regions:
1. `retrieveFrequencies(regions, loci)`: This function requests `html` pages with the allele frequency table from the [allele frequencies website](http://www.allelefrequencies.net/) and converts them into dataframes. A new request (e.g. `http://www.allelefrequencies.net/hla6006a_scr.asp?hla_region=Europe&hla_locus=A&hla_show=>`) and thus dataframe is created per region and locus that is specified in the list parameters and will be converted to the following request parameters (see [the command line program for the available regions](#`-f`-additional-arguments):
 - '`hla_region`': Specificy from which region the allele frequencies should be displayed.
 - '`hla_locus`': Specificy from which locus the allele frequencies should be displayed.
 - '`hla_show`': Specificy whether only `negative`, `positive` or all allele frequencies (`all`) should be displayed (`all` by default, but can be changed by providing `positive` or `negative` as arguments).
 The average of the allele frequencies is then taken to obtain the region average for every allele. The resulting dataframes with region averages are then concatenated into one big dataframe with the following structure:
 
 allele | avg_frequency | total_sample_size | locus | region
 --- | --- | --- | --- | ---
 A\*01:01 | 0.07 | 503 | A | Australia

 In order to obtain global allele frequencies, the region averages are again averaged. This is to ensure that the frequencies of every region contribute evenly to the global average. If only a selection of regions is given, only frequencies of these regions will be averaged. The resulting dataframe is structured as follows:

 allele | avg_frequency | total_sample_size | locus | num_regions
 --- | --- | --- | --- | ---
 A\*01:01 | 0.06 | 3389704 | A | 11

 - `allele`: The allele name.
 - `avg_frequency`: The average frequency over either a single region or all regions. 
 - `total_sample_size`: The total number of individuals a given (average) frequency is based on. This is calculated by summing all the sample sizes of the individual studies. 
 - `locus`: The locus the allele belongs to.
 - `region`: The region over which all frequencies were averaged.
 - `num_region`: The number of average region frequencies that were used to calculate the global average.
 Both groups of dataframes are eventually returned as a tuple `(regionMean, globalMean)` with the regionMean consisting of a dictionary with per region a dataframe and globalMean being one big dataframe.
 
3. `filterFrequencies(regionMean, globalMean, specifiedAlleles, addHigher, addLower, noNull)`: This optional function takes both the region and global dataframes and removes frequencies of alleles which are not included in the list `specifiedAlleles`. In case either `addHigher` or `addLower` is set to true frequencies of either higher or lower resolution alleles of the `specifiedAlleles` will be preserved. If `noNull` is set to true, frequencies of null alleles are removed by default. The dataframes are again returned in the same type of tuple as by `retrieveFrequencies()`.
4. `saveFrequencies(regionMean, globalMean, outputFolder, outputDelimiter, outputExtension)`: This function takes the dataframes of either of the above functions and saves them to an output folder as `region_frequencies.<extension>` and `global_frequencies.<extension>`.

Example code:
```
#Retrieve the frequencies averaged per region and over all regions (globally)
regionMean, globalMean = retrieveFrequencies(regions, loci)

#Filter them and remove all null alleles, while keeping any higher resolution alleles
regionMean, globalMean = filterFrequencies(regionMean, globalMean, specifiedAlleles, addHigher=True, noNull=True)

#Save the tables
saveFrequencies(regionMean, globalMean, 'frequencies')
```
