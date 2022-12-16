#Import modules
import json, Bio, re, os, requests, pickle
from os.path import join, isdir
from collections import Counter
from copy import copy
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq, MutableSeq, UnknownSeq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Data.CodonTable import TranslationError
from io import StringIO

#Define constants (DPA2/B2 is left out because they do not produce a well defined amino acid sequence, as they contain stop codons).
defaultLoci = ['A', 'B', 'C', 'DRA', 'DRB1', 'DRB345', 'DQA1', 'DQA2', 'DQB1', 'DPA1', 'DPB1', 'MICA', 'MICB']
defaultUrl = 'https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/msf/'
defaultCodonStarts = {
    "A": -24,
    "B": -24,
    "C": -24,
    "DRA": -25,
    "DRB1": -29,
    "DRB3": -29,
    "DRB4": -29,
    "DRB5": -29,
    "DRB345": -29,
    "DQA1": -23,
    "DQA2": -23,
    "DQB1": -32,
    "DPA1": -31,
    "DPB1": -29,
    "MICA": -23,
    "MICB
}
HLALoci = list(defaultCodonStarts.keys())

#Function to retrieve the IMGT/HLA latest database version from their github
latestDbVersion = None
def getLatestDbVersion():
    
    #Retrieve, the latest db version if it has not already been determined
    global latestDbVersion
    if latestDbVersion is None:
        completeUrl = defaultUrl.replace('msf/', 'Allele_status.txt')
        text = requests.get(completeUrl).text
        latestDbVersion = re.search(r'IPD-IMGT/HLA ([\d\.]+)', text[0:300]).group(1).replace('.', '')

    #Print and return the latest dbVersion
    print(f'Looking for latest IMGT/HLA database version number: {latestDbVersion}')
    return latestDbVersion

#Function to retrieve the alignment files from the github of IMGTHLA and save it locally when a outFolder is given.
#A list of loci can either be provided or the default is used and a database version can also be given (default is 'Latest').
#In any case a dictionary of file strings is returned
def downloadAlignments(outFolder = None, molecule = 'nt', loci = defaultLoci, databaseVersion = 'Latest'):
    alignments = {}
    databaseVersion = str(databaseVersion).replace('.', '')
    if databaseVersion != 'Latest' and getLatestDbVersion() == databaseVersion:
        databaseVersion = 'Latest'

    for locus in loci:

        #Recontruct the file name
        if molecule in ['nt', 'nuc', 'nucleotide']:
            inbetween = 'nuc'
        elif molecule in ['aa', 'prot', 'amino acid']:
            inbetween = 'prot'
        else:
            raise ValueError("The molecule argument must be 'nt', 'nuc' or 'nucleotide' for nucleotide alignments or 'aa', 'prot' or 'amino acid' for protein alignments.")
        fileName = f'{locus}_{inbetween}.msf'
        completeUrl = defaultUrl.replace('Latest', databaseVersion.replace('.','')) + fileName
        print(f'Retrieve alignment for locus {locus} from:\n{completeUrl}')

        #Save the file per locus
        text = requests.get(completeUrl).text
        if '404: Not Found' in text:
            print(f'Locus {locus} alignment not found.')
        else:
            if outFolder is not None:
                with open(os.path.join(outFolder, fileName), 'w') as f:
                    f.write(text)
            alignments[locus] = text

    return alignments

#Function to load an object of the class below from a pickle file.
#When a file location is given, the pickle is loaded directly
#If a folder name is given, a dbVersion must also be given in order to download the pickle of the same name
#If it is not available it is downloaded from the internet and saved to the input location
#When noNull is enabled a pickle object with no null alleles is searched for instead.
#The files should be/are named the following way:
# - Alignment with null alleles: <dbVersion>-all.pickle
# - Alignment without null alleles: <dbVersion>-nonull.pickle
def loadAlignments(input, dbVersion = 'Latest', noNull = False):
    if os.path.isfile(input):
        with open(input, 'rb') as f:
            return pickle.load(f)
    elif os.path.isdir(input) or input == '':
        
        #If the latest one should be returned, retrieve the latest db version number
        if dbVersion == 'Latest':
            dbVersion = getLatestDbVersion()

        #Check if an HLA object already exists as a .pickle.
        #If so, load and return it
        whichAlleles = 'nonull' if noNull else 'all'
        picklePath = os.path.join(input, f'{dbVersion}-{whichAlleles}.pickle')
        if os.path.isfile(picklePath):
            with open(picklePath, 'rb') as f:
                return pickle.load(f)
        #Otherwise create a new one to be saved as a file and into the HLAObjects dictionary
        else:
            HLAObject = HLA(databaseVersion=dbVersion, noNull = noNull)
            HLAObject.save(input)
            return HLAObject
    else:
        raise FileNotFoundError(f"'{input}' is neither a file or folder.")

#Function to calculate frequences of characters at every position in a BioPython multiple sequence alignment (not locusAlignment!).
#It returns an array of dictionaries for every position with:
# -'count': a counter object: https://pymotw.com/2/collections/counter.html 
# -'alleles': a dictionary with codons mapped to allele names
def countAlignment(alignment, charLen = 1):
    
    #charLen must have a minimum of 1 and the total alignment length must be dividable by it
    alignmentLen = alignment.get_alignment_length()
    if charLen < 1:
        raise ValueError('charLen must have a value of at least 1.')
    if alignmentLen % charLen != 0:
        raise ValueError('The length of the alignment must be a multiple of charLen.')

    counts = []
    for charIndex in range(0, alignmentLen, charLen):
        count = Counter()
        alleles = {}
        for seqRec in alignment:
            chars = seqRec.seq[charIndex : charIndex + charLen]
            count.update({chars: 1})
            
            if chars not in alleles:
                alleles[chars] = []
            alleles[chars].append(seqRec.id)

        counts.append({'count': count, 'alleles': alleles})
    return counts

#Function to save Biopython's alignment objects
def saveAlignment(alignment, fileName, outputFolder, format):
    
    #Check if the output folder is indeed a folder
    if not isdir(outputFolder) and outputFolder != '':
        raise FileNotFoundError(f'{outputFolder} is not a valid directory.')
    
    #Save the alignment
    AlignIO.write(alignment, join(outputFolder, f'{fileName}.{format}'), format)

#Define class to represent all aligned HLA alleles, the used references and useful methods/operators
#When no alignment input is given, the alignment files per locus are retrieved from the internet, a databaseVersion can then be given to be used (the default is 'Latest').
#Otherwise it should be either a list of file names or a folder of files whose name starts with the locus name 
#(and either followed by an underscore '_' or dot '.')
#Whenever noNull is True, any null alleles are left out of the alignment.
class HLA:
    def __init__(self, alignmentInput = None, format = 'msf', codonStarts = defaultCodonStarts, databaseVersion = 'Latest', noNull = False):
        self.codonStart = defaultCodonStarts
        if databaseVersion == 'Latest':
            self.databaseVersion = getLatestDbVersion()
        else:
            self.databaseVersion = str(databaseVersion).replace('.', '')
        self.noNull = noNull

        #Setup two datastructures: One for the alignments per locus
        #and one with the name pointing to the SeqRec object within the above structure
        self.alignments = {}
        self.IDIndex = {}

        #If no input is given, retrieve the alignments from the internet
        if alignmentInput == None:
            for locus, text in downloadAlignments(databaseVersion=databaseVersion).items():
                handle = StringIO(text)
                self.alignments[locus] = AlignIO.read(handle, format = 'msf')

        #If it is a string to a folder, load a list of all files in the folder first
        #and assign it to 'AlignmentInput'
        elif isinstance(alignmentInput, str):
            flist = []
            for f in os.listdir(alignmentInput):
                path = os.path.join(alignmentInput, f)
                flist.append(path)
                if os.path.isfile(path):
                    locus = re.search(r'^(\w{1,6}?)(_|\.)', f).group(1)
                    alignment = AlignIO.read(path, format = format)
                    self.alignments[locus] = alignment

        #If it is a list of paths, loop over all of them and load the alignments
        elif isinstance(alignmentInput, list):
            for f in alignmentInput:
                locus = re.search(r'^(\w{1,6}?)(_|\.)', os.path.basename(f)).group(1)
                alignment = AlignIO.read(f, format = format)
                self.alignments[locus] = alignment

        #Loop over every alignment while replacing it with an locusAlignment class which support codon indexing
        #while saving the locus and ids of every sequence to the IDIndex (except for null alleles if noNull is True)
        for locus, alignment in list(self.alignments.items()):

            #Do save null alleles
            if not self.noNull:

                #Loop over all sequences in the alignment to assign and index to their IDs.
                for i, rec in enumerate(alignment):
                    self.IDIndex[rec.id] = (locus, i)
                self.alignments[locus] = LocusAlignment(alignment, locus, codonStarts[locus])

            #Do not save null alleles
            else:

                #Loop over all sequences in the alignment to assign index their IDs and filter out any null alleles (if noNull is True)
                offset = 0
                filteredAlignment = []
                for i, rec in enumerate(alignment):
                    if not rec.id.endswith('N'):
                        self.IDIndex[rec.id] = (locus, i+offset)
                        filteredAlignment.append(rec)
                    else:
                        offset -= 1
                self.alignments[locus] = LocusAlignment(MultipleSeqAlignment(filteredAlignment), locus, codonStarts[locus])

            #When the locus is 'DRB345' additionaly create new locusAlignments for every locus and add the appropiate alleles
            if locus == 'DRB345':
                for i in range(3, 6):
                    self.alignments[f'DRB{i}'] = self.alignments[locus].copy(MultipleSeqAlignment([]))

                for rec in alignment:
                    subLocus = re.match(r'\w{4}', rec.id).group(0)
                    self.alignments[subLocus].alignment.append(rec)

    #Private method to add an allele (not supported)
    def __setitem__(self, key, rec):
        raise TypeError('You cannot assign new alleles to this object.')

    #Private method to check the type of the key
    def __check_key__(self, key):
        if isinstance(key, int):
            return 'HLAIndex'
        elif isinstance(key, str):
            if re.match(r'HLA\d{5}', key):
                return 'HLAID'
            elif re.match(r'[\w]{1,6}\*\d{2,4}:\d{2,4}:\d{2,4}:\d{2,4}[^\W\d_]?|[\w]{1,6}\*\d{2,4}:\d{2,4}:\d{2,4}[^\W\d_]?|[\w]{1,6}\*\d{2,4}:\d{2,4}[^\W\d_]?|[\w]{1,6}\*\d{2,4}[^\W\d_]?', key):
                return 'HLAName'
            elif re.search(r'^[\w]{1,6}$', key):
                return 'HLALocus'
            #Raise a ValueError when the string is not a valid key
            else:
                raise ValueError('The string is not in a supported format.')
        elif isinstance(key, slice) or isinstance(key, tuple) or isinstance(key, list):
            return 'HLASlice'
        elif key is None:
            return None
        #For any other key type, raise a TypeError
        else:
            raise TypeError('HLA alleles can only be accessed by strings in the appropiate format.')
    
    #Private slice method ([]) to return the allele based on either 
    #the locus name (reference allele) or allele name (get first match when it does not fully match)
    #or a range of alleles either by comma separated names ([name1, name2]) or in between two of them ([name1:name2])
    #When index is enabled it will return the index instead for HLAName, but not the HLA locus.
    def __getitem__(self, key, returnIndex = False):
        #Return...
        keyType = self.__check_key__(key)
        if keyType == 'HLAIndex':                   #an allele with matching index (not supported)
            raise TypeError('HLA alleles can not be accessed by indices, only by key strings.')
        if keyType == 'HLAID':                      #an allele with matching ID (not supported)
            raise ValueError('HLAIDs cannot be used. Use the allele name/locus instead.')        
        elif keyType == 'HLAName':                #an allele with matching name
            try:
                locus, index = self.IDIndex[key]
                return self.alignments[locus][index] if not returnIndex else (locus, index)
            #If no exact name match is present, search for other HLA alleles instead which are either of higher or at least lower resolution 
            #and return those in a list or return just the indices
            except KeyError:
                locus = re.match('\w{1,6}', key).group(0)
                higher = []
                lower = []
                for i, rec in enumerate(self.align(locus).alignment):
                    #Search for higher res
                    if rec.id.startswith(key) and rec.id[len(key)] == ':':
                        if not returnIndex:
                            higher.append(rec)
                        else:
                            higher.append((locus, i))
                    #Search for lower res
                    elif key.startswith(rec.id) and key[len(rec.id)] == ':':
                        if not returnIndex:
                            lower.append(rec)
                        else:
                            lower.append((locus, i))                        

                #Preferably return the higher res, otherwise the lower res alleles.
                #If nothing was found, raise an error after all.
                if len(higher) > 0:
                    print(f'The allele name {key} did not return anything, so a list of higher resolution alleles is returned instead.')
                    return higher
                elif len(lower) > 0:
                    print(f'The allele name {key} did not return anything, so a list of lower resolution alleles is returned instead.')
                    return lower
                else:
                    raise KeyError(f'There is no HLA allele with the name \'{key}\' present, either of the same, higher or lower resolution.')
        elif keyType == 'HLALocus':                 #an reference allele with matching locus name
            try:
                if not returnIndex:
                    return self.alignments[key].ref
                else:
                    raise KeyError(f'You cannot retrieve the index for a locus. This is contained in <LocusAlignment>.ref.')
            except:
                raise KeyError(f'There is no {key} locus in the current allele reference list.')
        elif keyType == 'HLASlice':                 #multiple alleles in a HLASlice object when a slice is given
            return HLASlice(self, key)
        elif keyType is None:
            raise KeyError(f'You have to supply a allele name for a single key.')

    #Private method to return the number of HLA alleles present
    def __len__(self):
        return sum([len(al.alignment) for al in self.alignments.values()])

    #Shorthand method to return the aligment of a locus
    def align(self, locus):
        return self.alignments[locus]

    #Method to only return the allele names. When a locus is given, return the names of only that locus.
    #When returnDict is True, return a dictionary with loci as keys and lists of the names as values (giving a locus name will not do anything then).
    def names(self, optionalLocus = None, returnDict = False):
        if not returnDict:
            return [rec.id for locus, al in self.alignments.items() if locus not in ['DRB3', 'DRB4', 'DRB5'] for rec in al.alignment] if optionalLocus == None else [rec.id for rec in self.alignments[optionalLocus].alignment]
        else:
            alleleDict = {}
            for locus, al in self.alignments.items():
                for rec in al.alignment:
                    if locus in alleleDict:
                        alleleDict[locus].append(rec.id)
                    else:
                        alleleDict[locus] = [rec.id]
            return alleleDict

    #Method to save the current instance of the class as a pickle or alignments to the output folder
    #When specifying an alignment format, the alignments can be sliced by giving slice string as values 
    #for loci in sliceDict, such as {'A': '1-10,15,20-25', 'B': ...}.
    #The type of slicing can also be specified by sliceType, which must either be 'gaps', 'no_gaps' or 'amino',
    #to preserve gaps introduced in the reference sequence, remove them or translate everything to amino acids.
    def save(self, outputFolder, format = 'pickle', sliceDict = {}, sliceType = 'gaps'):

        #Check if the output folder is indeed a folder
        if not isdir(outputFolder) and outputFolder != '':
            raise FileNotFoundError(f'{outputFolder} is not a valid directory.')
        
        #Save the current object as a pickle
        if format == 'pickle':
            whichAlleles = 'nonull' if self.noNull else 'all'
            outputFile = os.path.join(outputFolder, f'{self.databaseVersion}-{whichAlleles}.pickle')
            with open(outputFile, 'wb') as f:
                pickle.dump(self, f)
                print(f'Pickle of HLA object succesfully saved to \'{outputFile}\'.')

        #Save the current object as separate alignment files
        else:

            #Convert the slice string to slice objects
            for locus, sliceString in sliceDict.items():
                try:
                    sliceObject = [int(rng) if not '-' in rng else slice(*[int(i) for i in rng.split('-')]) for rng in sliceString.split(',')]
                    sliceDict[locus] = sliceObject
                except:
                    raise ValueError("The slice strings in sliceDict are in an incorrect format, it should be: {'A': '1-10,15,20-25', 'B': ...}")

            #When at least a single allele is present in a locus alignment, slice the alignment if necessary and save it
            for locus, locusAlignment in self.alignments.items():
                if locus in sliceDict:
                    if sliceType == 'gaps':
                        alignment = locusAlignment.slice[sliceDict[locus]]
                    elif sliceType == 'no_gaps':
                        alignment = locusAlignment.cod[sliceDict[locus]]
                    elif sliceType == 'amino':
                        alignment = locusAlignment.amino[sliceDict[locus]]
                    else:
                        raise ValueError("The sliceType is not 'gaps', 'no_gaps' or 'amino'.")
                else:
                    alignment = locusAlignment.alignment

                saveAlignment(alignment, locus, outputFolder , format)


#Class for wrapping around a locus alignment.
#Slice operators can be used as a shortcut of <locusAlignment>.alignment[].
#When doing <locusAlignment>.cod[] specific codons in an aligment will be returned.
class LocusAlignment:

    def __init__(self, alignment, locus, codonStart):
        self.alignment = alignment
        self.ref = alignment[0]
        self.locus = locus
        self.codonStart = codonStart
        self.codonIndex = {}

        #Save the positions of the nucleotides of every codon into the codonIndex
        pos = 0
        codonPos = codonStart
        while pos < len(self.ref.seq):
            codon = []
            while len(codon) < 3 and pos < len(self.ref.seq):
                if self.ref.seq[pos] != '-':
                    codon.append(pos)
                pos += 1
            if len(codon) == 3:
                self.codonIndex[codonPos] = codon
            elif len(codon) in [1,2]:
                raise Exception(f'The sequence for the reference allele for locus {locus} is not a multiple of tree.')
            codonPos += 1 if codonPos != -1 else 2

        #Create a codon (fast) and aminoacid (slow) slicer as a property
        self.slice = self.__AlignmentSlicer__(alignment, self.codonIndex)
        self.cod = self.__CodonSlicer__(alignment, self.codonIndex, False)
        self.amino = self.__CodonSlicer__(alignment, self.codonIndex, True)

    #Pass down the use of slice operators to the alignment
    def __getitem__(self, key):
        return self.alignment[key]

    #Return the number of alleles in the alignment
    def __len__(self):
        return len(self.alignment)

    #Method to return the range of the codons
    def range(self):
        if len(self.codonIndex) > 0:
            return list(self.codonIndex.keys())[0], list(self.codonIndex.keys())[-1]
        else:
            return None

    #Return a copy of the locus alignment with an updated alignment
    def copy(self, newAlignment):
        newLocusAlignment = copy(self)
        newLocusAlignment.alignment = newAlignment
        newLocusAlignment.slice = copy(self.slice)
        newLocusAlignment.slice.alignment = newAlignment
        newLocusAlignment.cod = copy(self.cod)
        newLocusAlignment.cod.alignment = newAlignment
        newLocusAlignment.amino = copy(self.amino)
        newLocusAlignment.amino.alignment = newAlignment
        return newLocusAlignment

    #Create a general parent class for the slicers
    class __Slicer__:
        def __init__(self, alignment, codonIndex):
            self.alignment = alignment
            self.ref = alignment[0]
            self.codonIndex = codonIndex

    #Create an internal class which slices based on codon positions but preserves gaps (fast)
    class __AlignmentSlicer__(__Slicer__):

        #Return the sliced alignment based on the codon positions
        def __getitem__(self, sliceObject):

            #When the sliceObject is an actual slice, return an alignment with only the codons (with gaps) in between these two values
            if isinstance(sliceObject, slice):
                start = sliceObject.start
                stop = sliceObject.stop
                if 0 in [start, stop]:
                    raise ValueError('No codon 0 exists: Codons are numbered starting from 1 or -1.')

                if len(self.codonIndex) > 0:
                    
                    try:
                        startIndex = self.codonIndex[start][0] if start is not None else None
                    except:
                        startCodonPos = list(self.codonIndex.keys())[0]
                        stopCodonPos = list(self.codonIndex.keys())[-1]
                        raise KeyError(f'Codon {start} does not exist for this alignment. The codons range from {startCodonPos} to {stopCodonPos}')
                    
                    try:
                        stopIndex = self.codonIndex[stop-1][2] + 1 if stop is not None else None
                    except:
                        startCodonPos = list(self.codonIndex.keys())[0]
                        stopCodonPos = list(self.codonIndex.keys())[-1]
                        raise KeyError(f'Codon {stop} does not exist for this alignment. The codons range from {startCodonPos} to {stopCodonPos}')
               
                else:
                    raise KeyError(f'There are no valid codons in this alignment.')

                return self.alignment[:, startIndex:stopIndex]

            #When it is a tuple or list, call the method recursively on each subSliceObject
            elif isinstance(sliceObject, tuple) or isinstance(sliceObject, list):
                for i, subSliceObject in enumerate(sliceObject):
                    
                    subAlignment = self.__getitem__(subSliceObject)
                    if i == 0:
                        alignmentSum = subAlignment
                    else:
                        alignmentSum += subAlignment

                return alignmentSum

            #When it is an int it refers to a single codon, so return only this codon
            elif isinstance(sliceObject, int):
                codon = self.codonIndex[sliceObject]
                start = codon[0]
                stop = codon[2] + 1
                return self.alignment[:, start:stop]

    #Create an internal class for a codon slicer which can be translated or not (takes very long for long slices or lots of sequences)
    class __CodonSlicer__(__Slicer__):
        def __init__(self, alignment, codonIndex, translate = False):
            super().__init__(alignment, codonIndex)
            self.translate = translate

        #Return the the sliced codons (translated to aminoacids if self.translate is True) based on a slice object
        def __getitem__(self, sliceObject):

            #When the sliceObject is an actual slice, return the (translated) codons in between these two values
            if isinstance(sliceObject, slice):
                start = sliceObject.start
                stop = sliceObject.stop
                if 0 in [start, stop]:
                    raise ValueError('No codon 0 exists: Codons are numbered starting from 1 or -1.')

                if len(self.codonIndex) > 0:
                    try:

                        #Determine first and last codon positions
                        startCodonPos = list(self.codonIndex.keys())[0]
                        stopCodonPos = list(self.codonIndex.keys())[-1] + 1

                        #Assign these to start/stop when they are None
                        if start == None:
                            start = startCodonPos
                        if stop == None:
                            stop = stopCodonPos

                        #Retrieve (and translate) all the indices of all codons between start and stop
                        newAlignment = []
                        for allele in self.alignment:
                            sequence = []
                            for i in range(start, stop):
                                #Skip position 0
                                if i == 0:
                                    continue

                                #Try to translate and append an aminoacid/codon.
                                try:
                                    sequence.append(self.__get__codon__(i, allele, self.translate))
                                #When an aminoacid could not be translated, the sequence from that point onwards is just '-'
                                except TranslationError:
                                    sequence.extend(list((stop - i) * '-'))
                                    break

                            #Add the (translated) allele to the new alignment
                            newAlignment.append(SeqRecord(  
                                id = allele.id,
                                seq = Seq(''.join(sequence)),
                                description=''
                                ))       

                        return MultipleSeqAlignment(newAlignment)

                    except KeyError:
                        raise KeyError(f'Codon {i} does not exist for this alignment. The codons range from {startCodonPos} to {stopCodonPos}.')
                else:
                    raise KeyError(f'There are no valid codons in this alignment.')

            #When it is a tuple, call the method recursively on each subSliceObject
            elif isinstance(sliceObject, tuple) or isinstance(sliceObject, list):
                for i, subSliceObject in enumerate(sliceObject):
                    
                    subAlignment = self.__getitem__(subSliceObject)
                    if i == 0:
                        alignmentSum = subAlignment
                    else:
                        alignmentSum += subAlignment

                return alignmentSum

            #When it is an int it refers to a single codon, so return only this translated codon
            elif isinstance(sliceObject, int):
                return self.__getitem__(slice(sliceObject, sliceObject+1))
                        
        #Method to get a single codon (which can be translated)
        def __get__codon__(self, codonPos, allele, translate):
            #When a codon is not present raise a KeyError
            try:
                codon = self.codonIndex[codonPos]
            except:
                raise KeyError(f'Codon {codonPos} does not exist for this alignment.')

            #Try to construct the codon and translate it if enabled, if it fails raise an error
            codonStr = allele[codon[0]]  +  allele[codon[1]]  +  allele[codon[2]]
            if translate:
                try:
                    codonSeq = Seq(codonStr)
                    aminoacid = codonSeq.translate()
                except:
                    raise TranslationError(f'The codon {codonPos} could not be translated.')
                return str(aminoacid)
            else:
                return codonStr

#Class for slicing the main HLA class object into a smaller object containing a subset of HLA alleles
class HLASlice(HLA):
    def __init__(self, HLAObject, sliceObject):

        #Setup empty variables
        self.alignments = {}
        self.IDIndex = {}

        #Run a recursive method to add all of the HLA alleles based on the slice object
        self.__parse_slice__(HLAObject, sliceObject)

    #Private recursive method which splits tuples of slice objects and keys and adds any alleles based on the keys it comes across
    def __parse_slice__(self, HLAObject, sliceObject):

        #When the sliceObject is an actual slice, retrieve all the alleles in between these two values 
        if isinstance(sliceObject, slice):

            #Check the type of the slice, it may only be names and no locus. 
            #Optionally, only one of them can also be None to mean the first/last allele of the locus. 
            start = sliceObject.start
            startType = self.__check_key__(start)
            stop = sliceObject.stop
            stopType = self.__check_key__(stop)
            if startType in ['HLAID', 'HLAIndex', 'HLALocus'] or stopType in ['HLAID', 'HLAIndex', 'HLALocus']:
                raise ValueError(f'HLA IDs, indices or locus names cannot be used for slicing. Use the allele name instead.')

            #Retrieve the indices. Choose the first/last allele of a locus if one of them is None
            #Additionaly when the index returned is a list, so choose the first / last item as index for the start / stop, respectively.
            if startType is not None and stopType is not None:
                startReturn = HLAObject.__getitem__(start, True)
                stopReturn = HLAObject.__getitem__(stop, True)
                startLocus, startIndex = startReturn[0] if isinstance(startReturn, list) else startReturn
                stopLocus, stopIndex = stopReturn[-1] if isinstance(stopReturn, list) else stopReturn
            elif startType is None and stopType is not None:
                stopReturn = HLAObject.__getitem__(stop, True)
                startReturn = None
                stopLocus, stopIndex = stopReturn[-1] if isinstance(stopReturn, list) else stopReturn
                startLocus, startIndex = (stopLocus, 0)
            elif startType is not None and stopType is None:
                startReturn = HLAObject.__getitem__(start, True)
                stopReturn = None
                startLocus, startIndex = startReturn[0] if isinstance(startReturn, list) else startReturn
                stopLocus, stopIndex = (startLocus, len(HLAObject.align(startLocus)) - 1)
            else:
                raise ValueError('When slicing, at least either a start or end allele name must be provided.')

            #when the start and stop index are not in the same locus, raise an error
            if startLocus != stopLocus:
                raise ValueError('When slicing alleles they must be present in the same locus.')

            #When the index of the start index is bigger than the stop index, switch them around
            if startIndex > stopIndex:
                temp = startIndex
                startIndex = stopIndex if not isinstance(stopReturn, list) else stopReturn[0]
                stopIndex = temp if not isinstance(startReturn, list) else startReturn[-1]

            #Retrieve and save all alleles between these two indices, but extend the current alignment if it is already saved
            for i in range(startIndex, stopIndex + 1):
                self.__add_allele__(startLocus, i, HLAObject)

        #When it is a tuple or list, call the method recursively on each subSliceObject
        elif isinstance(sliceObject, tuple) or isinstance(sliceObject, list):
            for subSliceObject in sliceObject:
                self.__parse_slice__(HLAObject, subSliceObject)
        
        #Otherwise it is a single key, so save it
        else:
            #Retrieve the locus and allele differently depending on whether it is a name or a locus (retrieve reference).
            keyType = self.__check_key__(sliceObject)

            if keyType == 'HLAName':
                #When only one allele is returned, add it to the alignment
                indices = HLAObject.__getitem__(sliceObject, True)
                if not isinstance(indices, list):
                    locus, index = indices
                    self.__add_allele__(locus, index, HLAObject)
                #When multiple alleles are returned, also add them all to the alignment
                else:
                    locus = indices[0][0]
                    for index in indices:
                        self.__add_allele__(locus, index[1], HLAObject)
            elif keyType == 'HLALocus':
                locus = sliceObject
                self.__add_alignment__(locus, HLAObject.align(locus).ref, HLAObject)
            else:
                raise ValueError('HLA IDs or indices cannot be used for access. Use the allele name instead.')

    #Method to add an allele from the old HLAObject
    def __add_allele__(self, locus, index, HLAObject):
        oldLocusAlignment = HLAObject.alignments[locus]
        allele = oldLocusAlignment.alignment[index]

        #Add the allele and create a new alignment from the old one if not already present
        self.__add_alignment__(locus, allele, HLAObject)
        
        #Update the IDIndex
        self.IDIndex[allele.id] = (locus, len(self.alignments[locus]) - 1)

        #When the locus is 'DRB345' additionaly add them to either the 'DRB3, -4 or -5' alignment
        if locus == 'DRB345':
            subLocus = allele.id[:4]
            self.__add_alignment__(subLocus, allele, HLAObject)

    #Method to only create a new alignment and add an allele to it without modifying the self.IDIndex from the old HLAObject
    def __add_alignment__(self, locus, allele, HLAObject):
        oldLocusAlignment = HLAObject.alignments[locus]
        if locus not in self.alignments:
            newLocusAlignment = oldLocusAlignment.copy(MultipleSeqAlignment([allele]))
            self.alignments[locus] = newLocusAlignment
        else:
            self.alignments[locus].alignment.append(allele)
