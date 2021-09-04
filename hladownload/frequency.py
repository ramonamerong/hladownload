#Import general modules
import argparse, requests, re
import html
import pandas as pd
from os.path import join, isdir

#Define constants
testRun = False
regionOutputFile = 'region_frequencies'
globalOutputFile = 'global_frequencies'

#Setup url and parameters
frequencyURL = 'http://www.allelefrequencies.net/hla6006a_scr.asp'
frequencyLoci = ['A', 'B', 'C', 'DRB1', 'DQA1', 'DQB1', 'DPA1', 'DPB1']
regions = ['Australia', 'Europe', 'North Africa', 'North America', 'North-East Asia', 'Oceania', 'South and Central America', 'South Asia', 'South-East Asia', 'Sub-Saharan Africa', 'Western Asia']
params = {'hla_region': None, 'hla_locus': None, 'hla_show': '%3E'}

#Request header to prevent 403 error
requestHeader = {
  "User-Agent": "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/50.0.2661.75 Safari/537.36",
  "X-Requested-With": "XMLHttpRequest"
}

#Function to retrieve the frequencies from the allele frequencies website (http://www.allelefrequencies.net/).
#A tuple of dataframes is returned with all of the frequencies averaged per region and over all regions respectively.
#The samples sizes are also summed up. You can specify whether you want only 'positive', 'negative' or 'all' types of frequencies.
def retrieveFrequencies(regions, loci, hlaShow = 'all'):

    #When hlaShow is '<' or '>' assign it to the request parameters
    #when it is 'both' set it to None, otherwise throw an error
    if hlaShow == 'negative':
        params['hla_show'] = '='
    elif hlaShow == 'positive':
        params['hla_show'] = '>'
    elif hlaShow == 'all':
        params['hla_show'] = None
    else:
        raise ValueError("'hlaShow' should either be '<', '>' or 'both'.")

    #Loop over ever region and locus to retrieve all allele frequencies
    dfs = []
    for region in regions:
        for locus in loci:
            print(f'Retrieving allele frequencies of locus {locus} from region {region}...')

            #Retrieve the hmtl and parse it into a pandas dataframe
            params['hla_locus'] = locus
            params['hla_region'] = region
            result = requests.get(frequencyURL, headers=requestHeader, params=params, timeout=None)
            print(f'Frequencies succesfully retrieved from {result.url}')
            df = pd.read_html(result.text, attrs={'class':'tblNormal'})[0]

            #Filter the table to drop the 'Line' column (so it get not averaged)
            #and remove any rows with Nan values in the allele frequency column
            #df = df[df['Allele Frequency'] != 0]
            df.drop(columns='Line', inplace=True)
            df.dropna(axis = 0, how = 'any', subset = ['Allele Frequency'], inplace=True)
            
            #Average the frequencies per allele per region
            #Skip if no allele frequencies had been returned
            try:
                df = df.groupby('Allele').agg({'Allele Frequency': 'mean', 'Sample Size': 'sum'})
            except pd.core.base.DataError:
                continue

            #Expand the allele index into a column
            df.reset_index(level=['Allele'], inplace=True)

            #Rename the 'Allele', 'Allele Frequency' and 'Sample Size' columns
            df.rename(columns = {'Allele': 'allele', 'Allele Frequency': 'avg_frequency', 'Sample Size': 'total_sample_size'}, inplace=True)

            #Add locus and region columns
            df['locus'] = locus
            df['region'] = region

            #Append it to the list of dataframes
            dfs.append(df)

            #If this is a test run, break
            if testRun:
                break
        
        if testRun:
            break

    #Also average the allele frequencies over all regions, so each regions weighs the same
    regionMean = pd.concat(dfs)
    regionMean.reset_index(inplace=True, drop=True)
    try:
        globalMean = regionMean.groupby(['allele']).agg({'avg_frequency': 'mean', 'total_sample_size': 'sum', 'locus': lambda x: x.iloc[0], 'region': lambda x: len(x)})
        globalMean.rename(columns = {'region': 'num_regions'}, inplace = True)
    except pd.core.base.DataError:
        globalMean = regionMean

    #Expand the allele index into a column
    globalMean.reset_index(level=['allele'], inplace=True)

    return regionMean, globalMean


#Function to only preserve allele frequencies for alleles which can also be found in the alignment object. 
#When enabling 'addHigher' and/or 'addLower' any lower/higher resolution alleles with frequencies
#will also be reported. When noNull is set to true, any higher/lower resolution null alleles will be excluded
def filterFrequencies(regionMean, globalMean, specifiedAlleles, addHigher = False, addLower = False, noNull = False):

    #Create indexes for the tables by looping over them
    regionIndex = {}
    for i, row in regionMean.iterrows():

        #Skip any null alleles if noNull is True
        if noNull and row['allele'].endswith('N'):
            continue

        region = row['region']
        if region not in regionIndex:
            regionIndex[region] = {}

        fields = re.split(r'\*|:', row['allele'])
        currentDict = regionIndex[region]
        for j, field in enumerate(fields):

            if field not in currentDict:
                currentDict[field] = {'higher': {}, 'index': None}

            if j < len(fields) - 1:
                currentDict = currentDict[field]['higher']
            else:
                currentDict[field]['index'] = i

        currentDict = None

    # def printRes(allele):
    #     if allele['index'] is not None:
    #         print(regionMean.loc[allele['index']])
    #     for higherAllele in allele['higher'].values():
    #         printRes(higherAllele)

    # for allele in regionIndex['Australia'].values():
    #     printRes(allele)

    globalIndex = {}
    for i, row in globalMean.iterrows():

        #Skip any null alleles if noNull is True
        if noNull and row['allele'].endswith('N'):
            continue

        fields = re.split(r'\*|:', row['allele'])
        currentDict = globalIndex
        for j, field in enumerate(fields):

            if field not in currentDict:
                    currentDict[field] = {'higher': {}, 'index': None}

            if j < len(fields) - 1:
                currentDict = currentDict[field]['higher']
            else:
                currentDict[field]['index'] = i

    #Specify a function to look using the fields of a specified allele to look in the index dictionaries 
    #for all indices of higher/lower resulution alleles
    filteredRegion = set()
    filteredGlobal = set()

    def searchAlleles(fields, filteredSet, currDict):
        resolution = 0

        def searchHigherAlleles(allele):
            if allele['index'] is not None and addHigher:
                filteredSet.add(allele['index'])
            for higherAllele in allele['higher'].values():
                searchHigherAlleles(higherAllele)

        while len(currDict) > 0:

            #Try to retrieve an allele of an higher resolution
            try:
                allele = currDict[fields[resolution]]
            #If it fails, stop searching for indexes
            except KeyError:
                break
            
            #When the current resolution is that of the number of fields - 1, the allele has been found
            #add it then iterate over all it's higher resolution children to add them
            if resolution == len(fields) - 1:
                if allele['index'] is not None:
                    filteredSet.add(allele['index'])
                searchHigherAlleles(allele)
                break

            #Otherwise add the lower resolution's allele's index to the filtered set 
            #and continue seaching
            if allele['index'] is not None and addLower:
                filteredSet.add(allele['index'])
            currDict = allele['higher']
            resolution += 1
                

    #For every specified allele...
    for specifiedAllele in specifiedAlleles:

        #Split the allele into fields
        fields = re.split(r'\*|:', specifiedAllele)

        #Add indices higher/lower resolution allele frequencies to the filtered region set
        for region, currentDict in regionIndex.items():
            searchAlleles(fields, filteredRegion, currentDict)
               
        #Add indices higher/lower resolution allele frequencies to the filtered global set
        searchAlleles(fields, filteredGlobal, globalIndex)

    regionMean = regionMean.iloc[list(filteredRegion)]
    globalMean = globalMean.iloc[list(filteredGlobal)]

    return regionMean, globalMean


#Function to save the frequencies to the output folder into two files, representing the region and global average allele frequencies
def saveFrequencies(regionMean, globalMean, outputFolder, outputDelimiter, outputExtension):
    regionOutput = join(outputFolder, f'{regionOutputFile}.{outputExtension}')
    globalOutput = join(outputFolder, f'{globalOutputFile}.{outputExtension}')
    regionMean.to_csv(regionOutput, index=False, sep=outputDelimiter)
    globalMean.to_csv(globalOutput, index=False, sep=outputDelimiter)