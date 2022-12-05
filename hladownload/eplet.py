#Import general modules
import argparse, pandas, requests, re, sys
from os import listdir, getcwd
from os.path import join, isdir, isfile

#Import local modules
from .hla import *

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
    'Name': 'eplet',
    'Exposition': 'ellipro_score',
    'Description': 'polymorphic_residues',
    'Confirmation': 'antibody_reactivity',
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
        df.rename(columns=inputTableHeader, inplace=True)

        tables[locusGroup] = df
        print(f'{locusGroup} eplets succesfully retrieved from {url}')
    
    return tables

#Function to parse a string of eplet residues and return a dictionary of them in a list mapped by their position.
def parseEpletPositions(residues):
    epletPositions= {}
    residues = re.sub(r'\s|\(|\)|\.', '', residues)

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
            if epletPosition in epletPositions :
                epletPosition = ''
            epletPosition += residues[i]
        elif residue != '/':
            epletPosition = int(epletPosition) 
            if epletPosition not in epletPositions:
                epletPositions[epletPosition] = []
            epletPositions[epletPosition].append(residue)
            
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
#Optionally restrict for which locus of the locus groups eplets should be reported. 
#If not, the 'DRB345' locus is separated into DRB3/4/5 by default.
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

#Function to retrieve the eplets and their mapped alleles into different tables.
#The summary table describes general information for each eplet on each row.
#The residue table has one residue per row for each eplet.
#The allele table has one row per eplet / allele association.
#It can be specified with 'allEplets' whether eplets without the alleles associated with them should be removed from the table
def getMappedEpletTables(epletMaps, outputFolder, allEplets = True):
    
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

    return tables


#Overall function which uses the above ones to directly save the eplets
#When 'allEplets' is set to True, also eplets without associated alleles will be saved.
#Optionally restrict for which loci eplets should be reported with 'loci'.
def saveEplets(outputFolder, alignmentObject, extension = 'csv', delimiter = ',', allEplets = True, loci = None):
    epletMaps = {}    

    #Determine eplets of which locusgroup should be retrieved
    locusGroups = set()
    for locus in loci:
        locusGroups.add(locusToLocusGroup[locus])
    
    #Parse the eplet dataframes
    for locusGroup, df in retrieveEplets(list(locusGroups)).items():

        #Load the eplets
        eplets = readEplets(df)

        #Add the epletMap as a value under the locus group
        epletMaps[locusGroup] = mapEpletsToAlleles(eplets, alignmentObject, locusGroup, loci)

    #Save the tables from the eplet maps
    tables = getMappedEpletTables(epletMaps, outputFolder, allEplets = allEplets)
    for name, table in tables.items():
        table.to_csv(join(outputFolder, f'{name}-eplets.{extension}'), index=False, sep=delimiter)