#Import general modules
import argparse, requests
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
#The samples sizes are also summed up.
def retrieveFrequencies(regions, loci):
    
    #Loop over ever region and locus to retrieve all allele frequencies
    dfs = []
    for region in regions:
        for locus in loci:
            print(f'Retrieving allele frequencies of locus {locus} from region {region}...')

            #Retrieve the hmtl and parse it into a pandas dataframe
            params['hla_locus'] = locus
            params['hla_region'] = region
            result = requests.get(frequencyURL, headers=requestHeader, params=params, timeout=None)
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


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="This program is intended to retrieve and parse allele frequencies from the allele frequency database (http://www.allelefrequencies.net) and output them into two table files: One with the average allele frequency per region and the other globally. The total sample size it is also calculated per row.")
    parser.add_argument("-o", "--output_folder", help="File or directory to write the output to. Defaults to the current folder.", default = '')
    parser.add_argument("-e", "--output_extension", help="Extension that should be used for the output files. Defaults to 'csv'.", default = 'csv')
    parser.add_argument("-d", "--output_delimiter", help="Delimiter that should be used for the output file. Defaults to ','.", default = ',')
    parser.add_argument("-t", "--test_run", action='store_true', help=f"Flag to indicate that only allele frequencies of the first region ('{regions[0]}'') and locus ('{frequencyLoci[0]}') should be retrieved and saved for testing purposes.")
    args = parser.parse_args()
    testRun = args.test_run
    
    #Check if the output folder is indeed a folder
    if not isdir(args.output_folder) and args.output_folder != '':
        raise FileNotFoundError(f'{args.output_folder} is not a valid directory.')

    #Retrieve the frequencies averaged per region and over all regions (globally)
    regionMean, globalMean = retrieveFrequencies(regions, frequencyLoci, testRun)

    #Save the frequencies to the output folder into two files
    regionOutput = join(args.output_folder, f'{regionOutputFile}.{args.output_extension}')
    globalOutput = join(args.output_folder, f'{globalOutputFile}.{args.output_extension}')

    regionMean.to_csv(regionOutput, index=False, sep=args.output_delimiter)
    globalMean.to_csv(globalOutput, index=False, sep=args.output_delimiter)

    print(f'Allele frequencies saved to {args.output_folder}.')



