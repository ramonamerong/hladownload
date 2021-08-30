#Import general modules
import argparse, pandas, requests, re, sys
from os import listdir
from os.path import join, isdir, isfile

#Import local modules
sys.path.append('')
try:
    from _1_alignment_utility.scripts.alignment import *
except:
    pass

#Define constants
testRun = False

#Construct dictionary with urls per locus group to the eplet registry
epletLinks = {
    'ABC': "https://www.epregistry.com.br/index/databases/database/ABC/",
    'DRB': "https://www.epregistry.com.br/index/databases/database/DRB/",
    'DQ': "https://www.epregistry.com.br/index/databases/database/DQ/",
    'DP': "https://www.epregistry.com.br/index/databases/database/DP/",
}
#Request header to prevent 403 error
requestHeader = {
  "User-Agent": "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/50.0.2661.75 Safari/537.36",
  "X-Requested-With": "XMLHttpRequest"
}
#Table header to construct the input dataframes with. 
#If an column name is not included, it is removed from the dataframe.
inputTableHeader = {
    'Eplet': 'eplet',
    'ElliProScore': 'ellipro_score',
    'PolymorphicResidues': 'polymorphic_residues',
    'AntibodyReactivity': 'antibody_reactivity',
}
#Construct a dictionary with locus group keys pointing to a list of loci names
#and loci keys pointing to the locus group
#DRB345 is not included in the former dictionary but it is in the latter
locusGroupToLocus = {
    'ABC': ['A', 'B', 'C'],
    'DRB': ['DRB1', 'DRB3', 'DRB4', 'DRB5'],
    'DQ': ['DQA1', 'DQB1'],
    'DP': ['DPA1', 'DPB1'],
}
locusToLocusGroup = {
    'A': 'ABC',
    'B': 'ABC',
    'C': 'ABC',
    'DRB1': 'DRB',
    'DRB3': 'DRB',
    'DRB4': 'DRB',
    'DRB5': 'DRB',
    'DRB345': 'DRB',
    'DQA1': 'DQ',
    'DQB1': 'DQ',
    'DPA1': 'DP',
    'DPB1': 'DP'
}
epletLoci = list(locusToLocusGroup.keys())

#Table headers for the output files
summaryTableHeader = ['eplet', 'locus_group', 'loci', 'verified', 'residue_number', 'allele_number' ]
residueTableHeader = ['eplet', 'locus_group', 'position', 'residue']
alleleTableHeader = ['eplet', 'locus_group', 'locus', 'allele']

#Function to retrieve eplets from the urls in 'epletLinks` 
#and return a dictionary with the loci groups ABC, DRB, DQ and DP as keys and dataframes containing the eplets in the values.
#You can also specify optional loci for which the eplets must only be retrieved and returned
def retrieveEplets(locusGroups = None):
    tables = {}
    print('Retrieving eplets from eplet registry database:')

    if locusGroups == None:
        locusGroups = list(epletLinks.keys)
    for locusGroup in locusGroups:
        url = epletLinks[locusGroup]
        result = requests.get(url, headers=requestHeader)
        df = pandas.read_html(result.text, attrs={'class':'table table-bordered table-hover'})[0]
        
        #Delete any columns note in tableHeader and change the column names to its values
        for columnName in list(df.columns):
            if re.sub(r'\s', '', columnName) not in inputTableHeader:
                df.drop(columns=columnName, inplace=True)
        df.columns = inputTableHeader.values()

        tables[locusGroup] = df
        print(f'{locusGroup} eplets succesfully retrieved from {url}')
    
    return tables

#Function to save the dataframes retrieved by 'retrieveEplets' into files
def saveEplets(outputFolder, extension = 'csv', delimiter = ','):
    tables = retrieveEplets()
    for locusGroup, dataframe in tables.items():
        outputPath = join(outputFolder, f'{locusGroup}-eplets.{extension}')
        dataframe.to_csv(outputPath, sep=delimiter, index=False)
        print(f'{locusGroup} eplets saved to {outputPath}.')

#Function to parse a string of eplet residues and return a dictionary of them in a list mapped by their position.
def parseEpletPositions(residues):
    epletPositions= {}
    residues = re.sub(r'\s|\(|\)|', '', residues)

    #Basics: It is a number and an amino acid. 44R is an 'R' at position 44.
    #Complications:
    #  There are () which may be linkage disequilbrium(?) amino acids
    #  Sometimes they're separated by spaces.
    #    Strategy: Assume it's 100% linkage. Ignore/strip () and spaces in all cases. Assume they mean nothing.
    #  In exactly one case it's listed as  163L-167G/S So this one needs to be parsed specifically. it can be a G or S at 167.
    #Special cases
    epletPosition = ''
    for i in range(0, len(residues)):
        residue = residues[i]
        if residue.isnumeric():
            epletPosition += residues[i]
        elif residue != '/':
            epletPosition = int(epletPosition) 
            if epletPosition not in epletPositions:
                epletPositions[epletPosition] = []
            epletPositions[epletPosition].append(residue)
            epletPosition = ''

    return epletPositions

#Function to read eplets into a dictionary with per eplet key residues mapped to their positon as values in another dictionary.
def readEplets(input, delimiter=','):
    eplets = {}

    #When the input is a dataframe use it
    if isinstance(input, pandas.DataFrame):
        df = input
    #When the input is a file open and read it
    elif isfile(input):
        df = pandas.read_csv(input, sep=delimiter)

    #Otherwise raise an error
    else:
        raise Exception('Invalid input. It should either be a path to a file or a pandas dataframe with the correct columns.')

    #Fill nan values with 'NaN'
    df.fillna('NaN', inplace = True)

    #Create eplet dictionary to return
    for index, row in df.iterrows():
        epletName = row['eplet']
        epletResidues = parseEpletPositions(row['polymorphic_residues'])
        verified = row['antibody_reactivity'].strip() == 'Yes'
        eplets[epletName] = {'residues': epletResidues, 'verified': verified}
        
        #When this is a test run, only process the first eplet
        if testRun:
            break

    return eplets

#Function to map eplets to alleles.
#A eplet mapped dictionary is returned containing other locus mapped dictionaries with lists of associated alleles.
#Optionally restrict which eplets should be reported for each locus of a locus group. If not, the 'DRB345' locus is separated by default.
def mapEpletsToAlleles(eplets, alignmentObject, locusGroup, loci = None):
    if loci is None:
        loci = locusGroupToLocus[locusGroup]
    else:
        loci = [locus for locus in loci if locusToLocusGroup[locus] == locusGroup]
    print(f'Mapping eplets to alleles for loci {", ".join(loci)}')

    #Loop over every eplet
    numEplets = len(eplets.keys())
    for epletIndex, epletName in enumerate(eplets.keys()):

        #Print a progress message for every 5 eplets
        if (epletIndex + 1) % 5 == 0 :
            print(f'Mapping eplet {epletIndex+1} of {numEplets}...')
        eplets[epletName]['alleles'] = {}

        #Consider eplets for every locus separatily to identify alleles based on their alignment
        for locus in loci:

            #Make an sliced alignment of the locus containing only the amino acids of the eplet residue positions.
            epletPositions = eplets[epletName]['residues']
            locusAlignment = alignmentObject.align(locus)
            residuesAlignment = locusAlignment.amino[list(epletPositions.keys())]

            #Consider all the alleles to see if they are valid
            for allele in residuesAlignment:
                residues = str(allele.seq)
                alleleValid = True

                #An allele contains an eplet, whenever all the allele residues are the same as the eplet residues.
                for i, (position, epletResidues) in enumerate(epletPositions.items()):
                    alleleResidue = residues[i]
                    if alleleResidue not in epletResidues:
                        alleleValid = False
                        break

                #Add a valid allele
                if alleleValid:
                    if locus not in eplets[epletName]['alleles']:
                        eplets[epletName]['alleles'][locus] = []
                    eplets[epletName]['alleles'][locus].append(allele.id)
            
        #When this is a test run, only process the first eplet
        if testRun:
            break

    return eplets

#Function to save the eplets and their mapped alleles into different tables.
#The summary table describes general information for each eplet on each row.
#The residue table has one residue per row for each eplet.
#The allele table has one row per eplet / allele association.
#It can be specified whether eplets without the alleles associated with them should be removed from the table
def saveParsedEplets(epletMaps, outputFolder, extension = 'csv', delimiter = ',', allEplets = True):
    
    #Check if the output folder is indeed a folder
    if not isdir(outputFolder) and outputFolder != '':
        raise FileNotFoundError(f'{outputFolder} is not a valid directory.')
    
    #Loop over all the eplets of every locus group to create the tables
    summaryTable = []
    residueTable = []
    alleleTable = []
    for locusGroup, epletMap in epletMaps.items():
        for epletName, eplet in epletMap.items():
            
            #Retrieve eplet information
            verified = eplet['verified']
            epletPositions = eplet['residues']
            alleles = eplet['alleles']

            #Add row to the summary table for each eplet after determining the total number of associated alleles
            alleleNum = 0
            for locusAlleles in alleles.values():
                alleleNum += len(locusAlleles)
            if alleleNum > 0 or allEplets:
                summaryTable.append([epletName, locusGroup, ','.join(list(alleles.keys())), verified, len(epletPositions), alleleNum])
            
            #Add seperate rows for each residue position contained in an eplet
            for epletPosition, epletResidue in epletPositions.items():
                residueTable.append([epletName, locusGroup, epletPosition, ','.join(epletResidue)])

            #Add seperate rows for each allele/eplet association to the allele table
            for locus, alleleList in alleles.items():
                for alleleName in alleleList:
                    alleleTable.append([epletName, locusGroup, locus, alleleName])
            
    #Save the tables to the output folder
    tables = {
        'summary_table': pandas.DataFrame(data=summaryTable, columns=summaryTableHeader),
        'residue_table': pandas.DataFrame(data=residueTable, columns=residueTableHeader),
        'allele_table':  pandas.DataFrame(data=alleleTable, columns=alleleTableHeader)
    }
    for name, table in tables.items():
        table.to_csv(join(outputFolder, f'{name}-eplets.{extension}'), index=False, sep=delimiter)

    return tables

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="This program is intended to retrieve and parse eplets from the eplet registry database (https://www.epregistry.com.br) and output them into three table files containing a summary, allele-eplet associations and residues per eplet.")
    parser.add_argument("-I", "--input_directory", help="Directory with eplet input files. They should be named like 'ABC-eplets.csv'. When it is not provided, the files will automatically be retrieved from the internet.")
    parser.add_argument("-i", "--input_delimiter", help="Delimiter that is used in the input file. Defaults to ','.", default = ',')
    parser.add_argument("-O", "--output_directory", help="Directory to save the eplet output files. Defaults to the current folder.", default = '')
    parser.add_argument("-o", "--output_delimiter", help="Delimiter that should be used for the output files. Defaults to ','.", default = ',')
    parser.add_argument("-e", "--output_extension", help="Extension that should be used for the output files. Defaults to 'csv'.", default = 'csv')
    parser.add_argument("-A", "--alignment_script_directory", help="Directory with the 'alignment.py' script. When that script is in the same folder as this one, this argument can be left empty.")
    parser.add_argument("-a", "--alignment_directory", help="Directory to search for alignment pickles generated by the '<HLA>' class from 'alignment.py'. When no alignment can be found a new one will be generated instead and also saved in the provided folder. Defaults to the current folder.", default = '')
    parser.add_argument("-D", "--database_version", help="The IMGT/HLA database version to use to determine which alleles contain which eplets.", default = 'Latest')
    parser.add_argument("-t", "--test_run", action='store_true', help="Flag to indicate that only the first eplet from each file should be processed for testing purposes.")
    args = parser.parse_args()
    
    testRun = args.test_run

    #Try to import the HLA script
    if args.alignment_script_directory != None:
        sys.path.append(args.alignment_script_directory)
        try:
            from alignment import *
        except:
            pass

    #Try to retrieve an alignment object:
    try:
        alignmentObject = loadAlignments(args.alignment_directory, noNull=True, dbVersion=args.database_version)
    except Exception as e:
        raise FileNotFoundError('The alignment.py script could not be found in the current folder or the provided directory.')

    #Retrieve the eplet files from the internet and parse them if an input folder is not given
    epletMaps = {}
    if args.input_directory == None:
    
        #Parse the eplet dataframes
        for locusGroup, df in retrieveEplets().items():

            #Load the eplets
            eplets = readEplets(df)

            #Add the epletMap as a value under the locus group
            epletMaps[locusGroup] = mapEpletsToAlleles(eplets, alignmentObject, locusGroup)

    #otherwise retrieve the file names from the folder and parse the eplets
    else:

        #Parse the eplet files
        for fileName in listdir(args.input_directory):
            try:

                #Determine the locus group from the file name
                locusGroup = re.match(r'[^-]+', fileName).group(0)
                if locusGroup is None:
                    raise Exception

                #Load the eplets
                eplets = readEplets(join(args.input_directory, fileName), delimiter = args.input_delimiter)

                #Add the epletMap as a value under the locus group
                epletMaps[locusGroup] = mapEpletsToAlleles(eplets, alignmentObject, locusGroup)

            except:
                raise FileNotFoundError('The input files don\'t have the appropiate columns or are not in the correct format: <locus group>-eplet.<extension>.')
    
    #Save the tables from the eplet maps
    saveParsedEplets(epletMaps, args.output_directory, args.output_extension, args.output_delimiter)


