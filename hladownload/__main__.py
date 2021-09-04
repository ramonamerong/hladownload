from pickle import load
import sys, os, re
from os import mkdir
from os.path import join, isdir, isfile
from Bio import AlignIO
from .hla import *
from .eplet import *
from .frequency import *

#https://packaging.python.org/tutorials/packaging-projects/
#https://gehrcke.de/2014/02/distributing-a-python-command-line-application/
#https://guides.github.com/features/mastering-markdown/
#setup.py: https://godatadriven.com/blog/a-practical-guide-to-using-setup-py/
#To do
# - Check whether alignmentCache should be deleted when updating the package

def main():
    parser = argparse.ArgumentParser(description="This tool can be used to automatically download HLA allele sequencies (flag -H), eplets/epitopes (flag -E) and frequencies (flag -F) from various online sources.")
    alleleGroup = parser.add_argument_group(f'HLA allele alignment output arguments (-H)\nSupported loci: {", ".join(HLALoci)}')
    alleleGroup.add_argument("-H", "--HLA", action='store_true', help="Flag to indicate that HLA alleles should be saved as multiple sequence alignments in fasta format per locus. These allele frequencies are retrieved from the IMGT/HLA database github: https://github.com/ANHIG/IMGTHLA")
    alleleGroup.add_argument("-t", "--slice_type", choices=['no_gaps', 'gaps', 'amino'], help="Can only be used if '-H' is given as a flag. Indicate whether the sliced codons should include gaps, no gaps or should be translated to amino acids. These gaps are only referring to those that are introduced into the reference sequence due to the alignment. Gaps are included by default.", default = 'gaps')        
    alleleGroup.add_argument("-s", "--slice_range", metavar=('SLICE_LOCI', 'SLICE RANGES'), action='append', nargs=2, help='''Can only be used if '-H' is given as a flag. Indicate the codons of the HLA alleles per locus that should be preserved in the outputted alignment files. The arguments should consist of a comma-separated list of loci, 
    followed by a combination of ranges and comma separated numbers for the colons to preserve: e.g. '-s A,B,C 1-5,10' means codons 1 up to (but not including) 5 and 10 of loci A, B and C should be preserved. You can also specify different codons for different loci by repeating the argument, e.g. '-s A,B 1-5 -s DRB1,DRB345 5-8 etc.'.
    By default all codons of the aligned alleles of a locus will be preserved if that locus is not specified.''')    

    epletGroup = parser.add_argument_group(f'Eplet output arguments (-E)\nSupported loci: {", ".join(epletLoci)}')
    epletGroup.add_argument("-E", "--eplet", action='store_true', help="Flag to indicate that eplets should retrieved from the eplet registry database (https://www.epregistry.com.br) and outputted them into three csv files containing a summary (summary_table-eplets.csv), allele-eplet associations (allele_table-eplets.csv) and residues per eplet (residue_table-eplets.csv).")
    epletGroup.add_argument("-l", "--all_eplets", action='store_true', help="Flag to indicate that all eplets should be included, as any eplets which are not associated with any (specified) alleles are removed from the output by default.")

    frequencyGroup = parser.add_argument_group(f'HLA frequency output arguments (-F)\nSupported loci: {", ".join(frequencyLoci)}')
    frequencyGroup.add_argument("-F", "--frequency", action='store_true', help="Flag to indicate that HLA allele frequencies should be retrieved from the allele frequency database (http://www.allelefrequencies.net) and outputted into two table files: One with the average allele frequency per region and the other globally (averaged over all (specified) regions). The total sample size is also calculated per row.")
    frequencyGroup.add_argument("-r", "--region", nargs='*', metavar='REGION', choices=regions, default=regions, help=f'Indicate the regions separated by spaces for which the HLA frequencies should be retrieved. A region consisting of multiple words should be enclosed in quotes. If no regions are specified all of the regions will be reported. If it is specified the global frequencies will only be based on the average of those regions. Usable regions are: {", ".join(regions)}')
    frequencyGroup.add_argument("->", "--higher_resolution", action='store_true', default=False, help=f'Flag to indicate that frequencies of higher resolution alleles of the specified alleles (see \'-a\' or \'-A\') should also be reported.')
    frequencyGroup.add_argument("-<", "--lower_resolution", action='store_true', default=False,  help=f'Flag to indicate that frequencies of lower resolution alleles of the specified alleles (see \'-a\' or \'-A\') should also be reported.')

    otherGroup = parser.add_argument_group('Other arguments')
    otherGroup.add_argument("-o", "--output_directory", help="Directory to save the output files to. For every output flag (-H, -E or -F) that is enabled, a separate directory (/alleles, /eplets or /frequencies) will be created in the output folder with the output files. The output folder defaults to the current folder.", default = '')
    otherGroup.add_argument("-d", "--database_version", help="The IMGT/HLA database version github (https://github.com/ANHIG/IMGTHLA) to use to retrieve HLA alleles from. Defaults to the latest branch.", default = 'Latest')
    otherGroup.add_argument("-a", "--alleles", help="Indicate the HLA alleles for which an alignment, frequencies and/or eplets should be outputted. The alleles can be a combination of ranges and comma separated allele names without spaces: 'A*01:01:01:01-A*01:01:01:09,B*01:01:01:01'. The default is all alleles for all loci that are included in the given IMGT/HLA database version.")    
    otherGroup.add_argument("-A", "--alleles_file", help="Indicate in a file the HLA alleles for which an alignment, frequencies and/or eplets should be outputted. The allele names should be separated by new lines in the file. This argument overrides the alleles specified with '-a'.")        
    otherGroup.add_argument("-n", "--no_null", action='store_true', help="Flag to indicate whether HLA null alleles should be excluded.")        
    otherGroup.add_argument("-g", "--group_DRB345", action='store_true', help="Flag to indicate whether HLA alleles of the DRB3, -4 and -5 loci should be grouped under a single locus 'DRB345'.")        

    args = parser.parse_args()

    outputExtension = 'csv'
    outputDelimiter = ','


    #Only continue if at least one output type is given and the output directory exists
    if not args.HLA and not args.eplet and not args.frequency:
        exit('Specify at least one output type with the flags -H, -E or -F. Use -h for additional help.')
    if args.output_directory != '' and not isdir(args.output_directory):
        exit(f"'{args.output_directory}' is not an existing output folder.")

    #Get the list of alleles
    alleles = None
    if args.alleles_file is not None and isfile(args.alleles_file):
        with open(args.alleles_file) as file:
            alleles = file.read().splitlines()
    elif args.alleles is not None:
        alleles = args.alleles.split(',')
        for i, allele in enumerate(alleles):
            if '-' in allele:
                alleles[i] = tuple(allele.split('-'))

    #Create a folder for the alignment pickles if it does not already exist
    scriptFolder = os.path.abspath(os.path.dirname(__file__))
    alignmentCacheFolder = join(scriptFolder, 'alignmentCache')
    if not isdir(alignmentCacheFolder):
        mkdir(alignmentCacheFolder)

    #Load an alignment picle from the cache folder or create a new one
    print('\nLoading HLA locus alignments...')
    alignmentObject = loadAlignments(alignmentCacheFolder, dbVersion=args.database_version, noNull=args.no_null)
    
    #Print an error when no HLA alleles are present
    if len(alignmentObject.names()) == 0:
        exit(f'No HLA alleles could be retrieved for the the IMGT/HLA database version \'{args.database_version}\'.')

    #When alleles are specified, trim down the alignment to only those alleles
    if alleles is not None:
        print('\nFiltering specified alleles...')
        try:
            alignmentObject = alignmentObject[alleles]
        except ValueError as e:
            exit('The supplied allele names are not in the correct format, e.g: A*01:01:01:01')
        except KeyError as e:
            exit(str(e).strip('"'))

    #Make a list of all loci which contain alleles but do not consider locus DRB345 but only the DRB3, -4 and -5 loci
    loci = []
    for locus, locusAlignment in alignmentObject.alignments.items():
        if len(locusAlignment.alignment) > 0 and ( (locus != 'DRB345' and not args.group_DRB345) or (locus not in ['DRB3', 'DRB4', 'DRB5'] and args.group_DRB345) ):
            loci.append(locus)

    #Generate allele output if it should be given
    if args.HLA:
        alleleOutputFolder = join(args.output_directory, 'alleles')
        print(f'\nSaving HLA allele locus alignments to \'{alleleOutputFolder}\'...')

        #Create a directory for the allele output
        if not isdir(alleleOutputFolder):
            mkdir(alleleOutputFolder)

        #Determine which alignments need to be sliced
        sliceDict = {}
        if args.slice_range is not None:
            for sliceArgs in args.slice_range:
                try:
                    sliceLoci = sliceArgs[0].split(',')
                    sliceString = sliceArgs[1]
                    for sliceLocus in sliceLoci:
                        sliceDict[sliceLocus] = sliceString
                except:
                    exit('The arguments for the --slice_range parameter are incorrect.')

        #Save the aligments as fasta and slice it if necessary 
        alignmentObject.save(alleleOutputFolder, 'fasta', sliceDict, args.slice_type)


    #Generate eplet output if it should be given
    if args.eplet:
        epletOutputFolder = join(args.output_directory, 'eplets')
        print(f'\nSaving eplets that are present in the specified alleles to \'{epletOutputFolder}\' (this can take a while)...')

        #Determine which of the provided loci supports eplets.
        #Print a message for unsupported ones.
        epletSupportedLoci = []
        for locus in loci:
            if locus in epletLoci:
                epletSupportedLoci.append(locus)
            else:
                print(f'Locus {locus} does not support eplet output.')

        #Only continue if there are loci available
        if len(epletSupportedLoci) > 0:

            #Create a directory for the eplet output
            if not isdir(epletOutputFolder):
                mkdir(epletOutputFolder)

            #Save the eplets for the specified and valid loci
            saveEplets(epletOutputFolder, alignmentObject, outputExtension, outputDelimiter, args.all_eplets, epletSupportedLoci)

    #Generate frequency output if it should be given
    if args.frequency:
        frequencyOutputFolder = join(args.output_directory, 'frequencies')
        print(f'\nSaving global and region frequencies of the specified alleles to \'{frequencyOutputFolder}\' (this can take a while)...')

        #Determine which of the provided loci supports frequencies.
        #Print a message for unsupported ones.
        frequencySupportedLoci = []
        for locus in loci:
            if locus in frequencyLoci:
                frequencySupportedLoci.append(locus)
            else:
                print(f'Locus {locus} does not support frequency output.')

        #Only continue if there are loci available
        if len(frequencySupportedLoci) > 0:

            #Create a directory for the frequency output
            if not isdir(frequencyOutputFolder):
                mkdir(frequencyOutputFolder)

            #Retrieve the frequencies averaged per region and over all regions (globally)
            regionMean, globalMean = retrieveFrequencies(args.region, frequencySupportedLoci)

            #Only preserve frequencies of the specified alleles or any of lower/higher resolution
            if alleles is not None:
                regionMean, globalMean = filterFrequencies(regionMean, globalMean, alleles, args.higher_resolution, args.lower_resolution, args.no_null)
                
            #Save the frequencies to the output folder into two files
            saveFrequencies(regionMean, globalMean, frequencyOutputFolder, outputDelimiter, outputExtension)

if __name__ == '__main__':
    main()