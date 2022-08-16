#!/usr/bin/env python3
from collections import defaultdict
from . import common, data_structures
import re
import os
import numpy as np
import sys
import statistics
import time
from Bio.Blast.Applications import NcbiblastnCommandline, NcbimakeblastdbCommandline
import textwrap
import argparse
from collections import defaultdict



def startScanPre(seq):
    CODON = "n"
    for i in range(len(seq), -3, -3):
        if i != len(seq):
            codon = seq[i:i+3]
            if codon in ["ATG", "GTG", "TTG"]:
                CODON = codon
    return CODON


def startScanPost(seq):
    CODON = "n"
    for i in range(0, len(seq), 3):
        if i != len(seq):
            codon = seq[i:i + 3]
            if codon in ["ATG", "GTG", "TTG"]:
                CODON = codon
                break
    return CODON


def stopScanPre(seq):
    CODON = "n"
    for i in range(len(seq), -3, -3):
        if i != len(seq):
            codon = seq[i:i+3]
            if codon in ["TAG", "TAA", "TGA"]:
                CODON = codon
    return CODON


def stopScanPost(seq):
    CODON = "n"
    for i in range(0, len(seq), 3):
        if i != len(seq):
            codon = seq[i:i + 3]
            if codon in ["TAG", "TAA", "TGA"]:
                CODON = codon
                break
    return CODON


def alnCheckN(seq1, seq2, slack):
    count = 0
    countGaps = 0
    for i in range(0, slack):
        try:
            if "-" not in [seq1[i], seq2[i]]:
                if seq1[i] == seq2[i]:
                    count += 1
            else:
                countGaps += 1
        except IndexError:
            break
    return count/slack, countGaps/slack


def alnCheckC(seq1, seq2, slack):
    count = 0
    countGaps = 0
    for i in range(len(seq1)-1, len(seq1)-slack, -1):
        if "-" not in [seq1[i], seq2[i]]:
            if seq1[i] == seq2[i]:
                count += 1
        else:
            countGaps += 1
    return count/slack, countGaps/slack


def AAI(seq1, seq2):
    counter = 0
    for i in range(len(seq1)):
        if "-" not in [seq1[i], seq2[i]]:
            if seq1[i] == seq2[i]:
                counter += 1
    return (counter/len(seq1))*100


def numerateList(ls):
    total = 0
    for i in ls:
        total += len(i)
    return total


def trimmer(seq1, seq2):
    ''' 
    seq1 needs to be the reference
    '''
    counter = 0
    newSeq1 = ""
    newSeq2 = ""
    for i in range(0, len(seq1), 3):
        if re.findall(r'-', seq1[i:i+3]) or re.findall(r'-', seq2[i:i+3]):
            if counter == 0:
                pass
            else:
                newSeq1 += seq1[i:i + 3]
                newSeq2 += seq2[i:i + 3]
        else:
            counter += 1
            newSeq1 += seq1[i:i+3]
            newSeq2 += seq2[i:i + 3]
    return newSeq1, newSeq2


def backTrim(seq1, seq2):
    ''' 
        seq1 needs to be the reference
        '''
    newSeq1 = seq1
    newSeq2 = seq2
    for i in range(len(seq1), 0, -3):
        if lastItem(newSeq2) == "-":
            newSeq1 = seq1[0:i]
            newSeq2 = seq2[0:i]
        else:
            break

    return newSeq1, newSeq2


def gapper(seq):
    gapList = []
    gapListFrame = []
    gaps = ''
    for i in seq:
        if i != '-':
            if len(gaps) > 0:
                if len(gaps) % 3:
                    gapList.append(gaps)
                    gaps = ''
                else:
                    gapListFrame.append(gaps)
                    gaps = ''
        else:
            gaps += i
    return gapList, gapListFrame


def degapper(seq1, seq2):

    """
    Seqs need to be the same length, 
    from an alignment
    """

    gaps1 = ''
    gaps2 = ''
    newSeq1 = ''
    newSeq2 = ''
    for i in range(len(seq1)):
        if "-" not in [seq1[i], seq2[i]]:
            if len(gaps1) > 0:
                if len(gaps1) % 3 == 0:
                    gaps1 = ''
                    gaps2 = ''
                    newSeq1 += seq1[i]
                    newSeq2 += seq2[i]
                else:
                    newSeq1 += gaps1
                    newSeq2 += gaps2
                    gaps1 = ''
                    gaps2 = ''
                    newSeq1 += seq1[i]
                    newSeq2 += seq2[i]
            else:
                newSeq1 += seq1[i]
                newSeq2 += seq2[i]
        else:
            gaps1 += seq1[i]
            gaps2 += seq2[i]
    return newSeq1, newSeq2


def lastItem(ls):
    x = ''
    for i in ls:
        if i != "":
            x = i
    return x


def startFinder(seq):
    start = ''
    for i in range(len(seq)):
        codon = (seq[i:i + 3])
        if not re.findall(r'-', codon):
            start = i
            break
    return start


def stopfinder(seq):
    newseq = ''
    gaps = ''
    for i in seq:
        if i != '-':
            newseq += gaps
            newseq += i
            gaps = ''
        else:
            gaps += i
    return newseq, len(gaps)


def diff(num1, num2):
    lower = sorted([num1, num2])[0]
    higher = sorted([num1, num2])[1]
    return higher - lower


def reverseComplement(seq):
    out = []
    for i in range(len(seq) - 1, -1, -1):
        nucleotide = seq[i]
        if nucleotide == "C":
            nucleotide = "G"
        elif nucleotide == "G":
            nucleotide = "C"
        elif nucleotide == "T":
            nucleotide = "A"
        elif nucleotide == "A":
            nucleotide = "T"
        out.append(nucleotide)
    outString = "".join(out)
    return outString


def Complement(seq):
    out = []
    for i in range(0, len(seq)):
        nucleotide = seq[i]
        if nucleotide == "C":
            nucleotide = "G"
        elif nucleotide == "G":
            nucleotide = "C"
        elif nucleotide == "T":
            nucleotide = "A"
        elif nucleotide == "A":
            nucleotide = "T"
        out.append(nucleotide)
    outString = "".join(out)
    return outString


def ribosome(seq):
    NTs = ['T', 'C', 'A', 'G']
    stopCodons = ['TAA', 'TAG', 'TGA']
    Codons = []
    for i in range(4):
        for j in range(4):
            for k in range(4):
                codon = NTs[i] + NTs[j] + NTs[k]
                # if not codon in stopCodons:
                Codons.append(codon)

    CodonTable = {}
    AAz = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    AAs = list(AAz)
    k = 0
    for base1 in NTs:
        for base2 in NTs:
            for base3 in NTs:
                codon = base1 + base2 + base3
                CodonTable[codon] = AAs[k]
                k += 1

    prot = []
    for j in range(0, len(seq), 3):
        codon = seq[j:j + 3]
        try:
            prot.append(CodonTable[codon])
        except KeyError:
            prot.append("")
    protein = ("".join(prot))
    return protein


def SeqCoord(seq, start, end):
    return seq[start:end]


def howMany(ls, exclude):
    counter = 0
    for i in ls:
        if i != exclude:
            counter += 1
    return counter


def stabilityCounter(int):
    if len(str(int)) == 1:
        string = (str(0) + str(0) + str(0) + str(0) + str(int))
        return (string)
    if len(str(int)) == 2:
        string = (str(0) + str(0) + str(0) + str(int))
        return (string)
    if len(str(int)) == 3:
        string = (str(0) + str(0) + str(int))
        return (string)
    if len(str(int)) == 4:
        string = (str(0) + str(int))
        return (string)
    if len(str(int)) > 4:
        string = str(int)
        return (string)


def sum(ls):
    count = 0
    for i in ls:
        count += float(i)
    return count


def ave(ls):
    count = 0
    for i in ls:
        try:
            count += float(i)
        except ValueError:
            pass
    return count / len(ls)


def derep(ls):
    outLS = []
    for i in ls:
        if i not in outLS:
            outLS.append(i)
    return outLS


def cluster(data, maxgap):
    '''Arrange data into groups where successive elements
       differ by no more than *maxgap*

        #->>> cluster([1, 6, 9, 100, 102, 105, 109, 134, 139], maxgap=10)
        [[1, 6, 9], [100, 102, 105, 109], [134, 139]]

        #->>> cluster([1, 6, 9, 99, 100, 102, 105, 134, 139, 141], maxgap=10)
        [[1, 6, 9], [99, 100, 102, 105], [134, 139, 141]]

    '''
    # data = sorted(data)
    data.sort(key=int)
    groups = [[data[0]]]
    for x in data[1:]:
        if abs(x - groups[-1][-1]) <= maxgap:
            groups[-1].append(x)
        else:
            groups.append([x])
    return groups


def GCcalc(seq):
    count = 0
    for i in seq:
        if i == "G" or i == "C":
            count += 1
    return count / len(seq)


def reject_outliers(data):
    m = 2
    u = np.mean(data)
    s = np.std(data)
    filtered = [e for e in data if (u - 2 * s < e < u + 2 * s)]
    return filtered


def RemoveDuplicates(ls):
    empLS = []
    counter = 0
    for i in ls:
        if i not in empLS:
            empLS.append(i)
        else:
            pass
    return empLS


def allButTheLast(iterable, delim):
    x = ''
    length = len(iterable.split(delim))
    for i in range(0, length - 1):
        x += iterable.split(delim)[i]
        x += delim
    return x[0:len(x) - 1]


def secondToLastItem(ls):
    x = ''
    for i in ls[0:len(ls) - 1]:
        x = i
    return x


def pull(item, one, two):
    ls = []
    counter = 0
    for i in item:
        if counter == 0:
            if i != one:
                pass
            else:
                counter += 1
                ls.append(i)
        else:
            if i != two:
                ls.append(i)
            else:
                ls.append(i)
                counter = 0
    outstr = "".join(ls)
    return outstr


def replace(stringOrlist, list, item):
    emptyList = []
    for i in stringOrlist:
        if i not in list:
            emptyList.append(i)
        else:
            emptyList.append(item)
    outString = "".join(emptyList)
    return outString


def remove(stringOrlist, list):
    emptyList = []
    for i in stringOrlist:
        if i not in list:
            emptyList.append(i)
        else:
            pass
    outString = "".join(emptyList)
    return outString


def removeLS(stringOrlist, list):
    emptyList = []
    for i in stringOrlist:
        if i not in list:
            emptyList.append(i)
        else:
            pass
    return emptyList


def fasta(fasta_file):
    count = 0
    seq = ''
    header = ''
    Dict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    for i in fasta_file:
        i = i.rstrip()
        if re.match(r'^>', i):
            count += 1
            if count % 1000000 == 0:
                print(count)

            if len(seq) > 0:
                Dict[header] = seq
                header = i[1:]
                # header = header.split(" ")[0]
                seq = ''
            else:
                header = i[1:]
                # header = header.split(" ")[0]
                seq = ''
        else:
            seq += i
    Dict[header] = seq
    # print(count)
    return Dict


def fasta2(fasta_file):
    count = 0
    seq = ''
    header = ''
    Dict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    for i in fasta_file:
        i = i.rstrip()
        if re.match(r'^>', i):
            count += 1
            if count % 1000000 == 0:
                print(count)

            if len(seq) > 0:
                Dict[header] = seq
                header = i[1:]
                header = header.split(" ")[0]
                seq = ''
            else:
                header = i[1:]
                header = header.split(" ")[0]
                seq = ''
        else:
            seq += i
    Dict[header] = seq
    # print(count)
    return Dict


def allButTheFirst(iterable, delim):
    x = ''
    length = len(iterable.split(delim))
    for i in range(1, length):
        x += iterable.split(delim)[i]
        x += delim
    return x[0:len(x) - 1]


def filter(list, items):
    outLS = []
    for i in list:
        if i not in items:
            outLS.append(i)
    return outLS


def filterRe(list, regex):
    ls1 = []
    ls2 = []
    for i in list:
        if re.findall(regex, i):
            ls1.append(i)
        else:
            ls2.append(i)
    return ls1, ls2


def delim(line):
    ls = []
    string = ''
    for i in line:
        if i != " ":
            string += i
        else:
            ls.append(string)
            string = ''
    ls = filter(ls, [""])
    return ls


def firstNonspace(ls):
    for i in ls:
        if i != "":
            break
    return i


def gc(seq):
    gc = 0
    for bp in seq:
        if bp == "C" or bp == "G":
            gc += 1
    return gc / len(seq)


def Dictparser(Dictionary):
    lowest = float(1000)
    for i in Dictionary:
        if float(Dictionary[i]) < float(lowest):
            lowest = Dictionary[i]
            key = i
    return [i, lowest]


def reverseComplement(seq):
    out = []
    for i in range(len(seq) - 1, -1, -1):
        nucleotide = seq[i]
        if nucleotide == "C":
            nucleotide = "G"
        elif nucleotide == "G":
            nucleotide = "C"
        elif nucleotide == "T":
            nucleotide = "A"
        elif nucleotide == "A":
            nucleotide = "T"
        out.append(nucleotide)
    outString = "".join(out)
    return outString


def Complement(seq):
    out = []
    for i in range(0, len(seq)):
        nucleotide = seq[i]
        if nucleotide == "C":
            nucleotide = "G"
        elif nucleotide == "G":
            nucleotide = "C"
        elif nucleotide == "T":
            nucleotide = "A"
        elif nucleotide == "A":
            nucleotide = "T"
        out.append(nucleotide)
    outString = "".join(out)
    return outString


def ribosome(seq):
    NTs = ['T', 'C', 'A', 'G']
    stopCodons = ['TAA', 'TAG', 'TGA']
    Codons = []
    for i in range(4):
        for j in range(4):
            for k in range(4):
                codon = NTs[i] + NTs[j] + NTs[k]
                # if not codon in stopCodons:
                Codons.append(codon)

    CodonTable = {}
    AAz = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    AAs = list(AAz)
    k = 0
    for base1 in NTs:
        for base2 in NTs:
            for base3 in NTs:
                codon = base1 + base2 + base3
                CodonTable[codon] = AAs[k]
                k += 1

    prot = []
    for j in range(0, len(seq), 3):
        codon = seq[j:j + 3]
        try:
            prot.append(CodonTable[codon])
        except KeyError:
            prot.append("X")
    protein = ("".join(prot))
    return protein


def SeqCoord(seq, start, end):
    return seq[start:end]


def howMany(ls, exclude):
    counter = 0
    for i in ls:
        if i != exclude:
            counter += 1
    return counter


def stabilityCounter(int):
    if len(str(int)) == 1:
        string = (str(0) + str(0) + str(0) + str(0) + str(int))
        return (string)
    if len(str(int)) == 2:
        string = (str(0) + str(0) + str(0) + str(int))
        return (string)
    if len(str(int)) == 3:
        string = (str(0) + str(0) + str(int))
        return (string)
    if len(str(int)) == 4:
        string = (str(0) + str(int))
        return (string)
    if len(str(int)) > 4:
        string = str(int)
        return (string)


def sum(ls):
    count = 0
    for i in ls:
        count += float(i)
    return count


def ave(ls):
    count = 0
    for i in ls:
        count += float(i)
    return count / len(ls)


def derep(ls):
    outLS = []
    for i in ls:
        if i not in outLS:
            outLS.append(i)
    return outLS


def cluster(data, maxgap):
    '''Arrange data into groups where successive elements
       differ by no more than *maxgap*

        #->>> cluster([1, 6, 9, 100, 102, 105, 109, 134, 139], maxgap=10)
        [[1, 6, 9], [100, 102, 105, 109], [134, 139]]

        #->>> cluster([1, 6, 9, 99, 100, 102, 105, 134, 139, 141], maxgap=10)
        [[1, 6, 9], [99, 100, 102, 105], [134, 139, 141]]

    '''
    # data = sorted(data)
    data.sort(key=int)
    groups = [[data[0]]]
    for x in data[1:]:
        if abs(x - groups[-1][-1]) <= maxgap:
            groups[-1].append(x)
        else:
            groups.append([x])
    return groups


def GCcalc(seq):
    count = 0
    for i in seq:
        if i == "G" or i == "C":
            count += 1
    return count / len(seq)


# def reject_outliers(data):
#     m = 2
#     u = np.mean(data)
#     s = np.std(data)
#     filtered = [e for e in data if (u - 2 * s < e < u + 2 * s)]
#     return filtered


def lastItem(ls):
    x = ''
    for i in ls:
        if i != "":
            x = i
    return x


def RemoveDuplicates(ls):
    empLS = []
    counter = 0
    for i in ls:
        if i not in empLS:
            empLS.append(i)
        else:
            pass
    return empLS


def allButTheLast(iterable, delim):
    x = ''
    length = len(iterable.split(delim))
    for i in range(0, length - 1):
        x += iterable.split(delim)[i]
        x += delim
    return x[0:len(x) - 1]


def secondToLastItem(ls):
    x = ''
    for i in ls[0:len(ls) - 1]:
        x = i
    return x


def pull(item, one, two):
    ls = []
    counter = 0
    for i in item:
        if counter == 0:
            if i != one:
                pass
            else:
                counter += 1
                ls.append(i)
        else:
            if i != two:
                ls.append(i)
            else:
                ls.append(i)
                counter = 0
    outstr = "".join(ls)
    return outstr


def replace(stringOrlist, list, item):
    emptyList = []
    for i in stringOrlist:
        if i not in list:
            emptyList.append(i)
        else:
            emptyList.append(item)
    outString = "".join(emptyList)
    return outString


def remove(stringOrlist, list):
    emptyList = []
    for i in stringOrlist:
        if i not in list:
            emptyList.append(i)
        else:
            pass
    outString = "".join(emptyList)
    return outString


def removeLS(stringOrlist, list):
    emptyList = []
    for i in stringOrlist:
        if i not in list:
            emptyList.append(i)
        else:
            pass
    return emptyList


def fasta(fasta_file):
    count = 0
    seq = ''
    header = ''
    Dict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    for i in fasta_file:
        i = i.rstrip()
        if re.match(r'^>', i):
            count += 1
            if count % 1000000 == 0:
                print(count)

            if len(seq) > 0:
                Dict[header] = seq
                header = i[1:]
                # header = header.split(" ")[0]
                seq = ''
            else:
                header = i[1:]
                # header = header.split(" ")[0]
                seq = ''
        else:
            seq += i
    Dict[header] = seq
    # print(count)
    return Dict


def fasta2(fasta_file):
    count = 0
    seq = ''
    header = ''
    Dict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    for i in fasta_file:
        i = i.rstrip()
        if re.match(r'^>', i):
            count += 1
            if count % 1000000 == 0:
                print(count)

            if len(seq) > 0:
                Dict[header] = seq
                header = i[1:]
                header = header.split(" ")[0]
                seq = ''
            else:
                header = i[1:]
                header = header.split(" ")[0]
                seq = ''
        else:
            seq += i
    Dict[header] = seq
    # print(count)
    return Dict


def allButTheFirst(iterable, delim):
    x = ''
    length = len(iterable.split(delim))
    for i in range(1, length):
        x += iterable.split(delim)[i]
        x += delim
    return x[0:len(x)]


def filter(list, items):
    outLS = []
    for i in list:
        if i not in items:
            outLS.append(i)
    return outLS


def filterRe(list, regex):
    ls1 = []
    ls2 = []
    for i in list:
        if re.findall(regex, i):
            ls1.append(i)
        else:
            ls2.append(i)
    return ls1, ls2


def delim(line):
    ls = []
    string = ''
    for i in line:
        if i != " ":
            string += i
        else:
            ls.append(string)
            string = ''
    ls = filter(ls, [""])
    return ls


def main():
    args = common.get_args('sleuth')

    ref_genes = args.ref_genes
    ref_gff = args.ref_gff
    target_genome = args.contigs
    outdir = args.outdir
    e = args.evalue
    threads = args.threads

    ffn = open(ref_genes)
    ffn = fasta2(ffn)

    genome = open(target_genome)
    genome = fasta2(genome)

    # SETTING LOCATION OF CODEML PARAMTER FILE
    ctl = os.path.dirname(os.path.dirname(__file__)) + "/modules/codeml-2.ctl"

    # MAKING DIRECTORIES
    os.system("mkdir -p %s" % outdir)
    os.system("mkdir -p %s/nuc_aln" % outdir)

    # READING IN RRNA AND TRNA FROM GFF FILE (FOR EXCLUSION LATER ON)
    rRNAdict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))

    if ref_gff:
        gff = open(args.ref_gff)
        for i in gff:
            if re.findall(r'FASTA', i):
                break
            else:
                if not re.match(r'#', i):
                    ls = i.rstrip().split("\t")
                    if ls[2] in ["rRNA", "tRNA"]:
                        try:
                            rRNAdict[ls[8].split("ID=")[1].split(";")[0]] = ls
                        except IndexError:
                            pass

                        try:
                            rRNAdict[ls[8].split("locus=")[1].split(";")[0]] = ls
                        except IndexError:
                            pass

                        try:
                            rRNAdict[ls[8].split("locus_tag=")[1].split(";")[0]] = ls
                        except IndexError:
                            pass

    else:
        out = open("%s-fixed.ffn" % (allButTheLast(args.ref_genes, ".")), "w")
        for i in ffn.keys():
            if not re.findall(r'tRNA-\.\.\.\(', i) and not re.findall(r'S ribosomal RNA', i):
                out.write(">" + i + "\n")
                out.write(ffn[i] + "\n")
        out.close()

        ffn = open("%s-fixed.ffn" % (allButTheLast(args.ref_genes, ".")))
        ffn = fasta2(ffn)

        ref_genes = "%s-fixed.ffn" % (allButTheLast(args.ref_genes, "."))


    # RUNNING BLAST
    print("Running BLAST")
    count = 0
    for i in genome.keys():
        if re.findall(r'\|', i):
            count += 1
    if count > 0:

        print("detected \'|\' characters in input fasta file. Creating a new version of the file %s-fixed.fna" % allButTheLast(target_genome, "."))
        out = open("%s-fixed.fna" % (allButTheLast(target_genome, ".")), "w")
        for i in genome.keys():
            out.write(">" + remove(i, ["|"]) + "\n")
            out.write(genome[i] + "\n")
        out.close()
        genome = open("%s-fixed.fna" % (allButTheLast(target_genome, ".")))
        genome = fasta2(genome)

        os.system("makeblastdb -dbtype nucl -in %s-fixed.fna -out %s-fixed.fna > /dev/null 2>&1" % (allButTheLast(target_genome, "."), allButTheLast(target_genome, ".")))
        os.system("blastn -query %s -db %s-fixed.fna -outfmt 6 -out %s/cds.genome.blast -evalue %s -perc_identity %s -qcov_hsp_perc %s -num_threads %s > /dev/null 2>&1" % (
            ref_genes, allButTheLast(target_genome, "."), outdir, str(e), str(float(args.perc_id)*100), str(float(args.perc_cov)*100), args.threads))

    else:
        os.system("makeblastdb -dbtype nucl -in %s -out %s > /dev/null 2>&1" % (target_genome, target_genome))
        os.system("blastn -query %s -db %s -outfmt 6 -out %s/cds.genome.blast -evalue %s -perc_identity %s -qcov_hsp_perc %s -num_threads %s > /dev/null 2>&1" % (
        ref_genes, target_genome, outdir, e, str(float(args.perc_id)*100), str(float(args.perc_cov)*100), args.threads))

    print("Done with BLAST\n.")
    # PARSING THE BLAST OUTPUT AND WRITING SEQUENCE FILES FOR ALIGNMENT
    diffList = []
    blastDict = defaultdict(list)
    blastDict2 = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    aniDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    counter = 0
    blast = open("%s/cds.genome.blast" % outdir)
    for i in blast:
        ls = i.rstrip().split("\t")
        if ls[0] not in rRNAdict.keys():
            alnregionlength = (int(ls[7]) - int(ls[6]) + 1)
            targetlength = diff(int(ls[8]), int(ls[9]))
            querylength = len(ffn[ls[0]])
            diffList.append(querylength - targetlength)

            slack = querylength * 0.1
            slack = round(slack)

            start = int(ls[8]) - 1
            end = int(ls[9])
            if start > end:
                seq = reverseComplement(genome[ls[1]][end - slack:start + slack])
            else:
                seq = genome[ls[1]][start - slack:end + slack]

            START = sorted([start, end])[0]
            END = sorted([start, end])[1]

            if ls[0] not in blastDict.keys():
                blastDict[ls[0]].append(ls)

                if counter > 0:
                    counter = 0
                    blastDict2.pop(key, None)

                    lengthQuery = len(ffn[key])
                    blastDict3 = defaultdict(list)
                    for k in blastDict[key]:
                        blastDict3[k[1]].append(k)

                    for l in blastDict3.keys():
                        if len(blastDict3[l]) > 1:
                            coordList = []
                            coordDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
                            for m in blastDict3[l]:
                                coordList.append(int(m[8]))
                                coordList.append(int(m[9]))
                                coordDict[int(m[8])] = m
                                coordDict[int(m[9])] = m

                            clus = (cluster(sorted(coordList), lengthQuery))
                            for j in clus:
                                target = l + "-" + str(j[0]) + "-" + str(lastItem(j))
                                start = int(coordDict[j[0]][8])
                                end = int(coordDict[j[0]][9])

                                if start > end:
                                    seq = reverseComplement(
                                        genome[l][j[0] - 1 - slack:lastItem(j) - 1 + slack])
                                else:
                                    seq = genome[l][j[0] - 1 - slack:lastItem(j) + slack]

                                aniList = []
                                for n in j:
                                    aniList.append(float(coordDict[n][2]))
                                ANI = statistics.mean(aniList)
                                blastDict2[key][target] = seq
                                aniDict[key][target] = ANI
                        else:
                            target = l + "-" + str(blastDict3[l][0][8]) + "-" + str(blastDict3[l][0][9])
                            start = int(blastDict3[l][0][8])
                            end = int(blastDict3[l][0][9])

                            if start > end:
                                seq = reverseComplement(
                                    genome[l][end - 1 - slack:start - 1 + slack])
                            else:
                                seq = genome[l][start - 1 - slack:end + slack]

                            blastDict2[key][target] = seq
                            aniDict[key][target] = str(blastDict3[l][0][2])

                else:
                    target = ls[1] + "-" + ls[8] + "-" + ls[9]
                    blastDict2[ls[0]][target] = seq
                    aniDict[ls[0]][target] = ls[2]

            else:
                counter += 1
                blastDict[ls[0]].append(ls)
                key = ls[0]

    for i in blastDict2.keys():
        for j in blastDict2[i]:
            out = open("%s/nuc_aln/%s.ffn" % (outdir, i + "__" + j), "w")
            out.write(">" + i + "\n")
            out.write(ffn[i] + "\n")
            out.write(">" + j + "\n")
            out.write(str(blastDict2[i][j]) + "\n")
            out.close()

    # RUNNING MUSCLE ON NUCLEOTIDE PAIRS
    print("starting Muscle")
    count = 0
    for file in (os.listdir(outdir + "/nuc_aln")):
        if lastItem(file.split(".")) == "ffn":
            count += 1
    total = count
    count = 0
    for file in os.listdir(outdir + "/nuc_aln"):
        if lastItem(file.split(".")) == "ffn":
            os.system("muscle -in %s/nuc_aln/%s -out %s/nuc_aln/%s.fa > /dev/null 2>&1" % (outdir, file, outdir, allButTheLast(file, ".")))
            count += 1
            perc = (count / total) * 100
            sys.stdout.write("running Muscle: %d%%   \r" % (perc))
            sys.stdout.flush()

    # PARSING THE ALIGNMENT FILES FOR DEEP ANALYSIS OF PSEUDOGENIZATION
    count = 0
    summaryDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    alnDir = "%s/nuc_aln" % outdir
    for i in os.listdir(alnDir):
        if re.findall(r'.fa', i):
            file = open("%s/%s" % (alnDir, i))
            file = fasta2(file)
            if len(file.keys()) < 2:
                pass  # wtf?
            else:

                # raw nucleotide alignment
                ref = (list(file.keys())[0])
                refSeq = file[ref]
                target = (list(file.keys())[1])
                targetSeq = file[target]

                querylength = len(remove(refSeq, ["-"]))
                slack = querylength * 0.1
                slack = round(slack)

                # cutting off the leading gaps in front of start codon
                start = (startFinder(refSeq))
                refAbridged = (refSeq[start:])
                targetAbridged = (targetSeq[start:])
                targetLeading = [targetSeq[0:start+ round(querylength*float(args.perc_cov))], targetSeq[start:start+round(querylength*float(args.perc_cov))]]

                # cutting off the trailing gaps after stop codon
                refAbridged2 = (stopfinder(refAbridged)[0])
                end = (stopfinder(refAbridged)[1])
                targetAbridged2 = targetAbridged[0:len(targetAbridged) - end]
                targetTrailing = [targetAbridged[len(targetAbridged)-end:len(targetAbridged)],
                                  targetAbridged[len(targetAbridged) - end - round(querylength*float(args.perc_cov)):len(targetAbridged)-end]]

                NterminalStats = (alnCheckN(refAbridged2, targetAbridged2, round(slack)))
                percIdentN = NterminalStats[0]
                gapsN = NterminalStats[1]

                lengthDiff = len(remove(targetAbridged2, ["-"])) / len(remove(refAbridged2, ["-"]))
                if lengthDiff > float(args.perc_cov):

                    if gapsN < 0.25 and percIdentN > 0.5 and AAI(refAbridged2, targetAbridged2) > float(args.perc_id):

                        # looking for start and stop codons in the *reference* sequence
                        START = refAbridged[0:3]
                        END = refAbridged2[len(refAbridged2) - 3:len(refAbridged2)]

                        if START in ["ATG", "GTG", "TTG"] and END in ["TAG", 'TAA', "TGA"]:
                            count += 1

                            # Nterminal =

                            # writing first setup file for apparent dN/dS using correct alignment
                            setup = open(ctl)
                            out = open("%s/%s.ctl" % (alnDir, i), "w")

                            for line in setup:
                                if re.findall('seqfile', line):
                                    out.write(
                                        '' + ' ' + '' + ' ' + '' + ' ' + '' + ' ' + '' + ' ' + '' + ' ' + 'seqfile' + ' ' + '= ' +
                                        alnDir + '/aln' + str(
                                            count) + '.ffn ' + '*' + ' ' + 'sequence' + ' ' + 'data' + ' ' + 'filename\n')

                                elif re.findall(r'outfile', line):
                                    out.write(
                                        '' + ' ' + '' + ' ' + '' + ' ' + '' + ' ' + '' + ' ' + '' + ' ' + 'outfile' + ' ' + '=' + ' '
                                        + alnDir + '/mlcTree_' + str(
                                            i) + ' ' + '' + ' ' + '' + ' ' + '' + ' ' + '' + ' ' + '' + ' ' + '' + ' ' + '' + ' '
                                        + '' + ' ' + '' + ' ' + '' + ' ' + '*' + ' ' + 'main' + ' ' + 'result' + ' ' + 'file' + ' ' + 'name\n')

                                else:
                                    out.write(line)
                            out.close()

                            # writing second setup file for no-mercy dN/dS (incorrect value, but represents reality)
                            setup = open(ctl)
                            out = open("%s/%s2.ctl" % (alnDir, i), "w")

                            for line in setup:
                                if re.findall('seqfile', line):
                                    out.write(
                                        '' + ' ' + '' + ' ' + '' + ' ' + '' + ' ' + '' + ' ' + '' + ' ' + 'seqfile' + ' ' + '= ' +
                                        alnDir + '/aln2-' + str(
                                            count) + '.ca.ffn ' + '*' + ' ' + 'sequence' + ' ' + 'data' + ' ' + 'filename\n')

                                elif re.findall(r'outfile', line):
                                    out.write(
                                        '' + ' ' + '' + ' ' + '' + ' ' + '' + ' ' + '' + ' ' + '' + ' ' + 'outfile' + ' ' + '=' + ' '
                                        + alnDir + '/mlcTree2_' + str(
                                            i) + ' ' + '' + ' ' + '' + ' ' + '' + ' ' + '' + ' ' + '' + ' ' + '' + ' ' + '' + ' '
                                        + '' + ' ' + '' + ' ' + '' + ' ' + '*' + ' ' + 'main' + ' ' + 'result' + ' ' + 'file' + ' ' + 'name\n')

                                else:
                                    out.write(line)
                            out.close()

                            Cterm = "y"
                            Nterm = "y"
                            if targetAbridged2[0] == "-":
                                # reference blast hit did not encompass entire gene...stupid BLAST...
                                # trimming alignment to only keep the inner aligned region, in frame
                                newSeqs = trimmer(refAbridged2, targetAbridged2)
                                refAbridged2 = newSeqs[0]
                                targetAbridged2 = newSeqs[1]
                                Nterm = "n"

                            if targetAbridged2[len(targetAbridged2) - 1] == "-":
                                newSeqs = backTrim(refAbridged2, targetAbridged2)
                                refAbridged2 = newSeqs[0]
                                targetAbridged2 = newSeqs[1]
                                Cterm = "n"

                            # detecting in-frame and frameshifting insertsions and deletions
                            outframeInsert = (gapper(refAbridged2)[0])
                            inframeInsert = (gapper(refAbridged2)[1])
                            outframeDel = (gapper(targetAbridged2)[0])
                            inframeDel = (gapper(targetAbridged2)[1])

                            # counting deletions and insertions
                            numInframeInserts = len(inframeInsert)
                            numOutframeInserts = len(outframeInsert)
                            totalInserts = numerateList(outframeInsert) + numerateList(inframeInsert)
                            numInframeDels = len(inframeDel)
                            numOutframeDels = len(outframeDel)
                            totalDels = numerateList(inframeDel) + numerateList(outframeDel)

                            # checking for start/stop codons
                            # start
                            if Nterm == "y":
                                if targetAbridged2[0:3] in ["ATG", "GTG", "TTG"]:
                                    startcodon = "y"
                                    realStartCodon = targetAbridged2[0:3]
                                else:
                                    anotherStartCodon = startScanPost(remove(targetLeading[1], ["-"]))
                                    if anotherStartCodon not in ["ATG", "TTG", "GTG"]:
                                        anotherStartCodon = startScanPre(remove(targetLeading[0], ["-"]))
                                        if anotherStartCodon not in ["ATG", "TTG", "GTG"]:
                                            startcodon = "n"
                                            realStartCodon = targetAbridged2[0:3]
                                        else:
                                            startcodon = "y"
                                            realStartCodon = anotherStartCodon
                                    else:
                                        startcodon = "y"
                                        realStartCodon = anotherStartCodon

                                preferredStartLoss = ''
                                preferredStartGain = ''
                                if START == "ATG":
                                    preferredStartGain = "n"
                                    if realStartCodon == "ATG":
                                        preferredStartLoss = "n"
                                    else:
                                        preferredStartLoss = "y"
                                else:
                                    preferredStartLoss = "n"
                                    if realStartCodon == "ATG":
                                        preferredStartGain = "y"
                                    else:
                                        preferredStartGain = "n"

                            else:
                                anotherStartCodon = startScanPost(remove(targetLeading[1], ["-"]))
                                if anotherStartCodon not in ["ATG", "TTG", "GTG"]:
                                    anotherStartCodon = startScanPre(remove(targetLeading[0], ["-"]))
                                    if anotherStartCodon not in ["ATG", "TTG", "GTG"]:
                                        startcodon = "n"
                                        preferredStartGain = "NA"
                                        preferredStartLoss = "NA"
                                    else:
                                        startcodon = "y"
                                        preferredStartGain = "NA"
                                        preferredStartLoss = "NA"

                                else:
                                    startcodon = "y"
                                    preferredStartGain = "NA"
                                    preferredStartLoss = "NA"

                            # stop
                            targetStraightDegapRaw = remove(targetAbridged2, ["-"])
                            if Cterm == "y":
                                putativeStopCodon = ""
                                for codon in range(0, len(targetStraightDegapRaw), 3):
                                    putativeStopCodon = targetStraightDegapRaw[codon:codon + 3]

                                if putativeStopCodon in ["TAG", 'TAA', "TGA"]:
                                    stopcodon = "y"
                                else:
                                    stopcodon = stopScanPost(remove(targetTrailing, ["-"]))
                                    if stopcodon != "n":

                                        stopcodon = "y"
                                    else:
                                        stopcodon = "n"
                            else:
                                stopcodon = stopScanPost(remove(targetTrailing, ["-"]))
                                if stopcodon != "n":
                                    stopcodon = "y"
                                else:
                                    stopcodon = "n"

                            # removing gaps from both sequences based on presence of a gap in the reference sequence.
                            # in other words, positions with insertions in the query sequence will be removed,
                            # presumably returning the alignment length to a multiple of 3!
                            refDegap1 = ''
                            targetDegap1 = ''
                            for j in range(len(refAbridged2)):
                                if refAbridged2[j] != "-":
                                    refDegap1 += refAbridged2[j]
                                    targetDegap1 += targetAbridged2[j]

                            # removing gap-containing codons from alignment - these are all codons with a gap, or deletion,
                            # in the query sequence. Removal of entire codons should keep the alignment a multiple of 3
                            refDegap2 = ''
                            targetDegap2 = ''
                            originalTargetSeq = ''
                            for j in range(0, len(targetDegap1), 3):
                                targetCodon = (targetDegap1[j:j + 3])
                                refCodon = (refDegap1[j:j + 3])
                                if not re.findall(r'-', targetCodon):  # this line is redundant and should do nothing
                                    if targetCodon not in ["TAG", 'TAA', "TGA"] and refCodon not in ["TAG", 'TAA',
                                                                                                      "TGA"]:
                                        targetDegap2 += targetCodon
                                        refDegap2 += refCodon
                                        originalTargetSeq += targetCodon

                            outFfn = open("%s/aln%s.ffn" % (alnDir, count), "w")
                            outFfn.write(">" + ref + "\n")
                            outFfn.write(refDegap2 + "\n")
                            outFfn.write(">" + target + "\n")
                            outFfn.write(targetDegap2 + "\n")
                            outFfn.close()

                            refAbridged3 = degapper(refAbridged2, targetAbridged2)[0]
                            targetAbridged3 = degapper(refAbridged2, targetAbridged2)[1]
                            refStraightDegap = remove(refAbridged3, ["-"])
                            targetStraightDegap = remove(targetAbridged3, ["-"])

                            refNoMercy = ""
                            targetNoMercy = ""
                            for j in range(0, sorted([len(refStraightDegap), len(targetStraightDegap)])[0], 3):
                                targetCodon = (targetStraightDegap[j:j + 3])
                                refCodon = (refStraightDegap[j:j + 3])
                                if targetCodon not in ["TAG", 'TAA', "TGA"] and refCodon not in ["TAG", 'TAA', "TGA"]:
                                    if len(targetCodon) == 3 and len(refCodon) == 3:
                                        targetNoMercy += targetCodon
                                        refNoMercy += refCodon

                            refNoMercyTrans = ribosome(refNoMercy)
                            targetNoMercyTrans = ribosome(targetNoMercy)

                            outFfnNoMercy = open("%s/aln2-%s.ffn" % (alnDir, count), "w")
                            outFfnNoMercy.write(">" + ref + "\n")
                            outFfnNoMercy.write(refNoMercy + "\n")
                            outFfnNoMercy.write(">" + target + "\n")
                            outFfnNoMercy.write(targetNoMercy + "\n")
                            outFfnNoMercy.close()

                            outFaaNoMercy = open("%s/aln2-%s.faa" % (alnDir, count), "w")
                            outFaaNoMercy.write(">" + ref + "\n")
                            outFaaNoMercy.write(refNoMercyTrans + "\n")
                            outFaaNoMercy.write(">" + target + "\n")
                            outFaaNoMercy.write(targetNoMercyTrans + "\n")
                            outFaaNoMercy.close()

                            stops = ribosome(targetStraightDegapRaw).count("*")
                            fragments = ribosome(targetStraightDegap).split("*")
                            firstFragment = fragments[0]
                            severity = len(fragments[0]) / len(ribosome(targetStraightDegap))

                            if stopcodon == "y":
                                stops = stops - 1

                            summaryDict[i]["numInframeInserts"] = numInframeInserts
                            summaryDict[i]["numOutframeInserts"] = numOutframeInserts
                            # summaryDict[i]["totalInserts"] = totalInserts / len(refDegap1)
                            summaryDict[i]["totalInserts"] = totalInserts
                            summaryDict[i]["numInframeDels"] = numInframeDels
                            summaryDict[i]["numOutframeDels"] = numOutframeDels
                            # summaryDict[i]["totalDels"] = totalDels / len(refDegap1)
                            summaryDict[i]["totalDels"] = totalDels
                            summaryDict[i]["startcodon"] = startcodon
                            summaryDict[i]["stopcodon"] = stopcodon
                            summaryDict[i]["stops"] = stops
                            summaryDict[i]["refAbridged2"] = refAbridged2
                            summaryDict[i]["refDegap2"] = refDegap2
                            summaryDict[i]["refNoMercy"] = refNoMercy
                            summaryDict[i]["targetAbridged2"] = targetAbridged2
                            summaryDict[i]["targetDegap2"] = targetDegap2
                            summaryDict[i]["targetNoMercy"] = targetNoMercy
                            summaryDict[i]["lengthDiff"] = lengthDiff
                            summaryDict[i]["stopSeverity"] = severity
                            summaryDict[i]["preferredStartLoss"] = preferredStartLoss
                            summaryDict[i]["preferredStartGain"] = preferredStartGain

    # # RUNNING CODEML
    count = 0
    for file in (os.listdir(outdir + "/nuc_aln")):
        if re.findall(r'aln2-', file):
            if lastItem(file.split(".")) == "ffn":
                count += 1
    total = count
    count = 0
    aaiDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    for file in os.listdir(alnDir):
        if re.findall(r'aln2-', file):
            if lastItem(file.split(".")) == "ffn":
                os.system("muscle -in %s/%s.faa -out %s/%s.faa.fa > /dev/null 2>&1" %
                          (alnDir, file.split(".")[0], alnDir, file.split(".")[0]))
                protAln = open("%s/%s.faa.fa" % (alnDir, file.split(".")[0]))
                protAln = fasta(protAln)
                ref = (list(protAln.keys())[0])
                refProt = protAln[ref]
                target = (list(protAln.keys())[1])
                targetProt = protAln[target]
                aai = AAI(refProt, targetProt)
                aaiDict[ref + "__" + target + ".fa"] = str(aai)
                os.system("pal2nal.pl %s/%s.faa.fa %s/%s -output fasta > %s/%s.ca.ffn" %
                          (alnDir, file.split(".")[0], alnDir, file, alnDir, file.split(".")[0]))
                count += 1
                perc = (count / total) * 100
                sys.stdout.write("running pal2nal: %d%%   \r" % (perc))
                sys.stdout.flush()
    print("")
    count = 0
    for file in (os.listdir(outdir + "/nuc_aln")):
        if lastItem(file.split(".")) == "ctl":
            count += 1
    total = count
    count = 0
    for file in os.listdir(alnDir):
        if lastItem(file.split(".")) == "ctl":
            os.system("codeml %s/%s > /dev/null 2>&1" % (alnDir, file))
            count += 1
            perc = (count / total) * 100
            sys.stdout.write("running codeml: %d%%   \r" % (perc))
            sys.stdout.flush()
    print("")
    print("done with codeml\n.")
    time.sleep(2)
    count = 0
    alnDir = "%s/nuc_aln" % outdir
    for i in os.listdir(alnDir):
        if re.findall(r'mlcTree_', i):
            count += 1
    print("%s out of %s input reference CDS had detectable homology (below evalue of %s)" %
          (str(count), str(len(ffn.keys())), str(e)))
    print("writing output file")


    # PARSING CODEML OUTPUT AND COMBINING WITH OTHER RESULTS
    fastaOut = open(outdir + "/ref_based_cds_predictions.ffn", "w")
    mainOut = open(outdir + "/sleuth_report.csv", "w")
    mainOut.write(",".join(data_structures.SleuthData._fields) + "\n")
    # mainOut.write("reference_locus,target_locus,ANI,AAI,ref_length,aln_query_cov,start,"
    #
    #               "loss_of_preferred_start,gain_of_preferred_start_codon,stop_codon,internal_stops,first_stop_codon,"
    #
    #               "out_of_frame_inserts,out_of_frame_dels,inframe_inserts,inframe_dels,total_inserts,total_deleted_bases,"
    #
    #               "ds,ds_no_mercy,dsds,delta_ds,dnds,dnds_no_mercy,"
    #
    #               "full_seq,mercy_aln,no_mercy_aln,full_ref_seq,mercy_aln_ref,no_mercy_aln_ref\n")

    for i in os.listdir(alnDir):
        if re.findall(r'mlcTree_', i):
            originalname = i.split("mlcTree_")[1]
            file = open("%s/%s" % (alnDir, originalname))
            file = fasta2(file)

            ref = (list(file.keys())[0])
            refSeq = file[ref]
            target = (list(file.keys())[1])
            targetSeq = file[target]

            ANI = aniDict[ref][target]

            ###############################################
            mlc = open("%s/%s" % (alnDir, i))
            for line in mlc:
                if not re.match(r'^ ', line):
                    try:
                        dn = (line.rstrip().split("dN = ")[1])
                        ds = (line.rstrip().split("dS = ")[1])
                    except IndexError:
                        pass

            ###############################################
            try:
                dn = float(dn.split(" ")[0])
                ds = float(ds.split(" ")[0])
            except AttributeError:
                dnds = "NA"

            if ds >= 0.001:
                dnds = dn / ds
            else:
                dnds = "NA"

            ###############################################
            print("%s/mlcTree2_%s" % (alnDir, i.split("Tree_")[1]))
            mlcNoMercy = open("%s/mlcTree2_%s" % (alnDir, i.split("Tree_")[1]))
            for line in mlcNoMercy:
                if not re.match(r'^ ', line):
                    try:
                        dnNoMercy = (line.rstrip().split("dN = ")[1])
                        dsNoMercy = (line.rstrip().split("dS = ")[1])
                    except IndexError:
                        pass

            ###############################################
            try:
                dnNoMercy = float(dnNoMercy.split(" ")[0])
                dsNoMercy = float(dsNoMercy.split(" ")[0])
            except AttributeError:
                dndsNoMercy = "NA"

            #########################################
            if dsNoMercy >= 0.001:
                try:
                    dndsNoMercy = dnNoMercy / dsNoMercy
                except ZeroDivisionError:
                    dndsNoMercy = "NA"
            else:
                dndsNoMercy = "NA"

            #########################################
            if ds > 0.001:
                dsds = dsNoMercy / ds
                deltads = dsNoMercy - ds
            else:
                dsds = "NA"
                deltads = dsNoMercy - ds

            numInframeInserts = summaryDict[originalname]["numInframeInserts"]
            numOutframeInserts = summaryDict[originalname]["numOutframeInserts"]
            totalInserts = summaryDict[originalname]["totalInserts"]
            numInframeDels = summaryDict[originalname]["numInframeDels"]
            numOutframeDels = summaryDict[originalname]["numOutframeDels"]
            totalDels = summaryDict[originalname]["totalDels"]
            startcodon = summaryDict[originalname]["startcodon"]
            stopcodon = summaryDict[originalname]["stopcodon"]
            stops = summaryDict[originalname]["stops"]
            refAbridged2 = summaryDict[originalname]["refAbridged2"]
            refDegap2 = summaryDict[originalname]["refDegap2"]
            refNoMercy = summaryDict[originalname]["refNoMercy"]
            targetAbridged2 = summaryDict[originalname]["targetAbridged2"]
            targetDegap2 = summaryDict[originalname]["targetDegap2"]
            targetNoMercy = summaryDict[originalname]["targetNoMercy"]
            lengthDiff = summaryDict[originalname]["lengthDiff"]
            severity = summaryDict[originalname]["stopSeverity"]
            preferredStartLoss = summaryDict[originalname]["preferredStartLoss"]
            preferredStartGain = summaryDict[originalname]["preferredStartGain"]

            aai = aaiDict[originalname]
            targetStraightDegapRaw = remove(targetAbridged2, ["-"])

            fastaOut.write(">" + target + "\n")
            fastaOut.write(targetStraightDegapRaw + "\n")
            mainOut.write(
                ref + "," + target + "," + str(ANI) + "," + str(aai) + "," + str(len(refAbridged2)) + "," + str(
                    lengthDiff) + "," +
                startcodon + "," + preferredStartLoss + "," + preferredStartGain + "," + stopcodon + "," + str(
                    stops) + "," + str(severity) + "," +
                str(numOutframeInserts) + "," + str(numOutframeDels) + "," + str(
                    numInframeInserts) + "," +
                str(numInframeDels) + "," + str(totalInserts) + "," + str(totalDels) + "," +
                str(ds) + "," + str(dsNoMercy) + "," + str(dsds) + "," + str(deltads) + "," + str(dnds) + "," + str(
                    dndsNoMercy) + "," + targetAbridged2 +
                "," + targetDegap2 + "," + targetNoMercy + "," + refAbridged2 + "," + refDegap2 + "," + refNoMercy + "\n")

    mainOut.close()
    fastaOut.close()

    os.system("rm -f 2NG.t 2NG.dN 2NG.dS rst1 rst 2ML.t 2ML.dN 2ML.dS 4fold.nuc rub")
    os.system("tar -cf %s/nuc_aln.tar %s/nuc_aln" % (outdir, outdir))
    os.system("gzip %s/nuc_aln.tar" % outdir)
    os.system("rm -r %s/nuc_aln" % outdir)
    os.system("rm %s.n*" % target_genome)


def full(args, file_dict, log_file_dict=None):
    if log_file_dict is None:
        target_genome = file_dict['contigs_filename']
    else:   # This would happen if reannotate is called and we need to find the genome from previously
        target_genome = log_file_dict['contigs_filename']

    ref_ffn = file_dict['ref_cds_filename']
    target_genome_blastdb = file_dict['temp_dir'] + "/target_contigs.fasta"
    ctl = file_dict['ctl']
    out = file_dict['sleuthDir']
    e = args.evalue
    threads = args.threads

    # MAKING DIRECTORIES
    os.system("mkdir -p %s" % out)
    os.system("mkdir -p %s/nuc_aln" % out)

    # RUNNING BLAST
    common.print_with_time("Running BLAST...")
    ffn = open(ref_ffn)
    ffn = fasta2(ffn)

    genome = open(target_genome)
    genome = fasta2(genome)
    count = 0
    for i in genome.keys():
        if re.findall(r'\|', i):
            count += 1
    if count > 0:

        common.print_with_time("detected \'|\' characters in input fasta file. Creating a new version of the file %s-fixed.fna" % allButTheLast(target_genome, "."))
        outFixed = open("%s/%s-fixed.fna" % (out, allButTheLast(target_genome, ".")), "w")
        for i in genome.keys():
            outFixed.write(">" + remove(i, ["|"]) + "\n")
            outFixed.write(genome[i] + "\n")
            outFixed.close()
        genome = open("%s/%s-fixed.fna" % (out, allButTheLast(target_genome, ".")))
        genome = fasta2(genome)

        os.system("makeblastdb -dbtype nucl -in %s/%s-fixed.fna -out %s/%s-fixed.fna > /dev/null 2>&1" % (
            out, allButTheLast(target_genome, "."), out, allButTheLast(target_genome, ".")))
        os.system("blastn -query %s -db %s/%s-fixed.fna -outfmt 6 -out %s/cds.genome.blast -evalue %s -num_threads %s -perc_identity %s -qcov_hsp_perc %s > /dev/null 2>&1" % (
            ref_ffn, out, allButTheLast(target_genome, "."), out, str(e), str(threads), str(float(args.perc_id)*100), str(float(args.perc_cov)*100)))

    else:
        os.system("makeblastdb -dbtype nucl -in %s -out %s > /dev/null 2>&1" %
                  (target_genome, target_genome_blastdb))
        os.system("blastn -query %s -db %s -outfmt 6 -out %s/cds.genome.blast -evalue %s -num_threads %s -perc_identity %s -qcov_hsp_perc %s > /dev/null 2>&1" %
                  (ref_ffn, target_genome_blastdb, out, e, str(threads), str(float(args.perc_id)*100), str(float(args.perc_cov)*100)))

    common.print_with_time("Done with BLAST.")
    # print("Processing...")
    # time.sleep(4)
    # PARSING THE BLAST OUTPUT AND WRITING SEQUENCE FILES FOR ALIGNMENT
    diffList = []
    blastDict = defaultdict(list)
    blastDict2 = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    aniDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    counter = 0
    blast = open("%s/cds.genome.blast" % out)
    redunDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    for i in blast:

        ls = i.rstrip().split("\t")
        identity = float(ls[2])
        alnregionlength = (int(ls[7]) - int(ls[6]) + 1)
        targetlength = diff(int(ls[8]), int(ls[9]))
        querylength = len(ffn[ls[0]])
        diffList.append(querylength - targetlength)

        slack = querylength * 0.1
        slack = round(slack)

        start = int(ls[8]) - 1
        end = int(ls[9])
        if start > end:
            seq = reverseComplement(genome[ls[1]][end - slack:start + slack])
        else:
            seq = genome[ls[1]][start - slack:end + slack]

        START = sorted([start, end])[0]
        END = sorted([start, end])[1]

        if ls[0] not in blastDict.keys():
            blastDict[ls[0]].append(ls)

            if counter > 0:
                counter = 0
                blastDict2.pop(key, None)

                lengthQuery = len(ffn[key])
                blastDict3 = defaultdict(list)
                for k in blastDict[key]:
                    blastDict3[k[1]].append(k)

                for l in blastDict3.keys():
                    if len(blastDict3[l]) > 1:
                        coordList = []
                        coordDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
                        for m in blastDict3[l]:
                            coordList.append(int(m[8]))
                            coordList.append(int(m[9]))
                            coordDict[int(m[8])] = m
                            coordDict[int(m[9])] = m

                        clus = (cluster(sorted(coordList), lengthQuery))
                        # frag = 0
                        for j in clus:
                            # if len(j) > 2:
                            #     frag = 1
                            targetStart = l + "-" + str(j[0])
                            if targetStart not in redunDict.keys():

                                target = l + "-" + str(j[0]) + "-" + str(lastItem(j))
                                start = int(coordDict[j[0]][8])
                                end = int(coordDict[j[0]][9])

                                if start > end:
                                    seq = reverseComplement(
                                        genome[l][j[0] - 1 - slack:lastItem(j) - 1 + slack])
                                else:
                                    seq = genome[l][j[0] - 1 - slack:lastItem(j) + slack]

                                aniList = []
                                for n in j:
                                    aniList.append(float(coordDict[n][2]))
                                ANI = statistics.mean(aniList)
                                blastDict2[key][target] = seq
                                aniDict[key][target] = ANI

                                redunDict[targetStart]["ID"] = identity
                                redunDict[targetStart]["target"] = target
                            else:
                                if identity > redunDict[targetStart]["ID"]:
                                    blastDict2[key].pop(redunDict[targetStart]["target"], None)

                                    target = l + "-" + str(j[0]) + "-" + str(lastItem(j))
                                    start = int(coordDict[j[0]][8])
                                    end = int(coordDict[j[0]][9])

                                    if start > end:
                                        seq = reverseComplement(
                                            genome[l][j[0] - 1 - slack:lastItem(j) - 1 + slack])
                                    else:
                                        seq = genome[l][j[0] - 1 - slack:lastItem(j) + slack]

                                    aniList = []
                                    for n in j:
                                        aniList.append(float(coordDict[n][2]))
                                    ANI = statistics.mean(aniList)
                                    blastDict2[key][target] = seq
                                    aniDict[key][target] = ANI

                                    redunDict[targetStart]["ID"] = identity
                                    redunDict[targetStart]["target"] = target
                                else:
                                    pass

                    else:
                        targetStart = l + "-" + str(blastDict3[l][0][8])
                        if targetStart not in redunDict.keys():

                            target = l + "-" + str(blastDict3[l][0][8]) + "-" + str(blastDict3[l][0][9])
                            start = int(blastDict3[l][0][8])
                            end = int(blastDict3[l][0][9])

                            if start > end:
                                seq = reverseComplement(
                                    genome[l][end - 1 - slack:start - 1 + slack])
                            else:
                                seq = genome[l][start - 1 - slack:end + slack]

                            blastDict2[key][target] = seq
                            aniDict[key][target] = str(blastDict3[l][0][2])

                            redunDict[targetStart]["ID"] = identity
                            redunDict[targetStart]["target"] = target
                        else:
                            if identity > redunDict[targetStart]["ID"]:
                                blastDict2[key].pop(redunDict[targetStart]["target"], None)

                                target = l + "-" + str(blastDict3[l][0][8]) + "-" + str(blastDict3[l][0][9])
                                start = int(blastDict3[l][0][8])
                                end = int(blastDict3[l][0][9])

                                if start > end:
                                    seq = reverseComplement(
                                        genome[l][end - 1 - slack:start - 1 + slack])
                                else:
                                    seq = genome[l][start - 1 - slack:end + slack]

                                blastDict2[key][target] = seq
                                aniDict[key][target] = str(blastDict3[l][0][2])

                                redunDict[targetStart]["ID"] = identity
                                redunDict[targetStart]["target"] = target
                            else:
                                pass

            else:
                targetStart = ls[1] + "-" + ls[8]
                if targetStart not in redunDict.keys():

                    target = ls[1] + "-" + ls[8] + "-" + ls[9]
                    blastDict2[ls[0]][target] = seq
                    aniDict[ls[0]][target] = ls[2]

                    redunDict[targetStart]["ID"] = identity
                    redunDict[targetStart]["target"] = target
                else:
                    if identity > redunDict[targetStart]["ID"]:
                        blastDict2[ls[0]].pop(redunDict[targetStart]["target"], None)

                        target = ls[1] + "-" + ls[8] + "-" + ls[9]
                        blastDict2[ls[0]][target] = seq
                        aniDict[ls[0]][target] = ls[2]

                        redunDict[targetStart]["ID"] = identity
                        redunDict[targetStart]["target"] = target
                    else:
                        pass

        else:
            counter += 1
            blastDict[ls[0]].append(ls)
            key = ls[0]

    for i in blastDict2.keys():
        for j in blastDict2[i]:
            outffn = open("%s/nuc_aln/%s.ffn" % (out, i + "__" + j), "w")
            outffn.write(">" + i + "\n")
            outffn.write(ffn[i] + "\n")
            outffn.write(">" + j + "\n")
            outffn.write(str(blastDict2[i][j]) + "\n")
            outffn.close()

    # RUNNING MUSCLE ON NUCLEOTIDE PAIRS
    common.print_with_time("Starting Muscle.")
    count = 0
    for file in (os.listdir(out + "/nuc_aln")):
        if lastItem(file.split(".")) == "ffn":
            count += 1
    total = count
    count = 0
    # os.system("for i in %s/nuc_aln/*ffn; do"
    #           " muscle -in $i -out $i.fa > /dev/null 2>&1;"
    #           " done" % out)
    for file in os.listdir(out + "/nuc_aln"):
        if lastItem(file.split(".")) == "ffn":
            os.system("muscle -in %s/nuc_aln/%s -out %s/nuc_aln/%s.fa > /dev/null 2>&1" % (out, file, out, allButTheLast(file, ".")))
            count += 1
            perc = (count / total) * 100
            sys.stdout.write(f"{common.current_time()}\tRunning Muscle: {round(perc)}%\r")
            sys.stdout.flush()
    print("")
    # for i in range(10):
    #     print(".")
    #     time.sleep(i / 25)
    common.print_with_time("Done with Muscle.")
    time.sleep(1)
    common.print_with_time("Preparing for CodeML...")

    # PARSING THE ALIGNMENT FILES FOR DEEP ANALYSIS OF PSEUDOGENIZATION
    count = 0
    summaryDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    alnDir = "%s/nuc_aln" % out
    for i in os.listdir(alnDir):
        if re.findall(r'.fa', i):
            file = open("%s/%s" % (alnDir, i))
            file = fasta2(file)
            if len(file.keys()) < 2:
                pass  # wtf?
            else:

                # raw nucleotide alignment
                ref = (list(file.keys())[0])
                refSeq = file[ref]
                target = (list(file.keys())[1])
                targetSeq = file[target]

                querylength = len(remove(refSeq, ["-"]))
                slack = querylength * 0.1
                slack = round(slack)

                # cutting off the leading gaps in front of start codon
                start = (startFinder(refSeq))
                refAbridged = (refSeq[start:])
                targetAbridged = (targetSeq[start:])
                targetLeading = [targetSeq[0:start + round(querylength * float(args.perc_cov))],
                                 targetSeq[start:start + round(querylength * float(args.perc_cov))]]

                # cutting off the trailing gaps after stop codon
                refAbridged2 = (stopfinder(refAbridged)[0])
                end = (stopfinder(refAbridged)[1])
                targetAbridged2 = targetAbridged[0:len(targetAbridged) - end]
                targetTrailing = [targetAbridged[len(targetAbridged) - end:len(targetAbridged)],
                                  targetAbridged[
                                  len(targetAbridged) - end - round(querylength * float(args.perc_cov)):len(targetAbridged) - end]]

                lengthDiff = len(remove(targetAbridged2, ["-"])) / len(remove(refAbridged2, ["-"]))

                if lengthDiff > float(args.perc_cov):
                    NterminalStats = (alnCheckN(refAbridged2, targetAbridged2, round(slack)))
                    percIdentN = NterminalStats[0]
                    gapsN = NterminalStats[1]

                    if gapsN < 0.25 and percIdentN > 0.5 and AAI(refAbridged2, targetAbridged2) > float(args.perc_id):

                        # looking for start and stop codons in the *reference* sequence
                        START = refAbridged[0:3]
                        END = refAbridged2[len(refAbridged2) - 3:len(refAbridged2)]

                        if START in ["ATG", "GTG", "TTG"] and END in ["TAG", 'TAA', "TGA"]:
                            count += 1

                            # Nterminal =

                            # writing first setup file for apparent dN/dS using correct alignment
                            setup = open(ctl)
                            outctl = open("%s/%s.ctl" % (alnDir, i), "w")

                            for line in setup:
                                if re.findall('seqfile', line):
                                    outctl.write(
                                        '' + ' ' + '' + ' ' + '' + ' ' + '' + ' ' + '' + ' ' + '' + ' ' + 'seqfile' + ' ' + '= ' +
                                        alnDir + '/aln' + str(
                                            count) + '.ffn ' + '*' + ' ' + 'sequence' + ' ' + 'data' + ' ' + 'filename\n')

                                elif re.findall(r'outfile', line):
                                    outctl.write(
                                        '' + ' ' + '' + ' ' + '' + ' ' + '' + ' ' + '' + ' ' + '' + ' ' + 'outfile' + ' ' + '=' + ' '
                                        + alnDir + '/mlcTree_' + str(
                                            i) + ' ' + '' + ' ' + '' + ' ' + '' + ' ' + '' + ' ' + '' + ' ' + '' + ' ' + '' + ' '
                                        + '' + ' ' + '' + ' ' + '' + ' ' + '*' + ' ' + 'main' + ' ' + 'result' + ' ' + 'file' + ' ' + 'name\n')

                                else:
                                    outctl.write(line)
                            outctl.close()

                            # writing second setup file for no-mercy dN/dS (incorrect value, but represents reality)
                            setup = open(ctl)
                            outctl = open("%s/%s2.ctl" % (alnDir, i), "w")

                            for line in setup:
                                if re.findall('seqfile', line):
                                    outctl.write(
                                        '' + ' ' + '' + ' ' + '' + ' ' + '' + ' ' + '' + ' ' + '' + ' ' + 'seqfile' + ' ' + '= ' +
                                        alnDir + '/aln2-' + str(
                                            count) + '.ca.ffn ' + '*' + ' ' + 'sequence' + ' ' + 'data' + ' ' + 'filename\n')

                                elif re.findall(r'outfile', line):
                                    outctl.write(
                                        '' + ' ' + '' + ' ' + '' + ' ' + '' + ' ' + '' + ' ' + '' + ' ' + 'outfile' + ' ' + '=' + ' '
                                        + alnDir + '/mlcTree2_' + str(
                                            i) + ' ' + '' + ' ' + '' + ' ' + '' + ' ' + '' + ' ' + '' + ' ' + '' + ' ' + '' + ' '
                                        + '' + ' ' + '' + ' ' + '' + ' ' + '*' + ' ' + 'main' + ' ' + 'result' + ' ' + 'file' + ' ' + 'name\n')

                                else:
                                    outctl.write(line)
                            outctl.close()

                            Cterm = "y"
                            Nterm = "y"
                            if targetAbridged2[0] == "-":
                                # reference blast hit did not encompass entire gene...stupid BLAST...
                                # trimming alignment to only keep the inner aligned region, in frame
                                newSeqs = trimmer(refAbridged2, targetAbridged2)
                                refAbridged2 = newSeqs[0]
                                targetAbridged2 = newSeqs[1]
                                Nterm = "n"

                            if targetAbridged2[len(targetAbridged2) - 1] == "-":
                                newSeqs = backTrim(refAbridged2, targetAbridged2)
                                refAbridged2 = newSeqs[0]
                                targetAbridged2 = newSeqs[1]
                                Cterm = "n"

                            # detecting in-frame and frameshifting insertsions and deletions
                            outframeInsert = (gapper(refAbridged2)[0])
                            inframeInsert = (gapper(refAbridged2)[1])
                            outframeDel = (gapper(targetAbridged2)[0])
                            inframeDel = (gapper(targetAbridged2)[1])

                            # counting deletions and insertions
                            numInframeInserts = len(inframeInsert)
                            numOutframeInserts = len(outframeInsert)
                            totalInserts = numerateList(outframeInsert) + numerateList(inframeInsert)
                            numInframeDels = len(inframeDel)
                            numOutframeDels = len(outframeDel)
                            totalDels = numerateList(inframeDel) + numerateList(outframeDel)

                            # checking for start/stop codons
                            # start
                            if Nterm == "y":
                                if targetAbridged2[0:3] in ["ATG", "GTG", "TTG"]:
                                    startcodon = "y"
                                    realStartCodon = targetAbridged2[0:3]
                                else:
                                    anotherStartCodon = startScanPost(remove(targetLeading[1], ["-"]))
                                    if anotherStartCodon not in ["ATG", "TTG", "GTG"]:
                                        anotherStartCodon = startScanPre(remove(targetLeading[0], ["-"]))
                                        if anotherStartCodon not in ["ATG", "TTG", "GTG"]:
                                            startcodon = "n"
                                            realStartCodon = targetAbridged2[0:3]
                                        else:
                                            startcodon = "y"
                                            realStartCodon = anotherStartCodon
                                    else:
                                        startcodon = "y"
                                        realStartCodon = anotherStartCodon

                                preferredStartLoss = ''
                                preferredStartGain = ''
                                if START == "ATG":
                                    preferredStartGain = "n"
                                    if realStartCodon == "ATG":
                                        preferredStartLoss = "n"
                                    else:
                                        preferredStartLoss = "y"
                                else:
                                    preferredStartLoss = "n"
                                    if realStartCodon == "ATG":
                                        preferredStartGain = "y"
                                    else:
                                        preferredStartGain = "n"

                            else:
                                anotherStartCodon = startScanPost(remove(targetLeading[1], ["-"]))
                                if anotherStartCodon not in ["ATG", "TTG", "GTG"]:
                                    anotherStartCodon = startScanPre(remove(targetLeading[0], ["-"]))
                                    if anotherStartCodon not in ["ATG", "TTG", "GTG"]:
                                        startcodon = "n"
                                        preferredStartGain = "NA"
                                        preferredStartLoss = "NA"
                                    else:
                                        startcodon = "y"
                                        preferredStartGain = "NA"
                                        preferredStartLoss = "NA"

                                else:
                                    startcodon = "y"
                                    preferredStartGain = "NA"
                                    preferredStartLoss = "NA"

                            # stop
                            targetStraightDegapRaw = remove(targetAbridged2, ["-"])
                            if Cterm == "y":
                                putativeStopCodon = ""
                                for codon in range(0, len(targetStraightDegapRaw), 3):
                                    putativeStopCodon = targetStraightDegapRaw[codon:codon + 3]

                                if putativeStopCodon in ["TAG", 'TAA', "TGA"]:
                                    stopcodon = "y"
                                else:
                                    stopcodon = stopScanPost(remove(targetTrailing, ["-"]))
                                    if stopcodon != "n":

                                        stopcodon = "y"
                                    else:
                                        stopcodon = "n"
                            else:
                                stopcodon = stopScanPost(remove(targetTrailing, ["-"]))
                                if stopcodon != "n":
                                    stopcodon = "y"
                                else:
                                    stopcodon = "n"

                            # removing gaps from both sequences based on presence of a gap in the reference sequence.
                            # in other words, positions with insertions in the query sequence will be removed,
                            # presumably returning the alignment length to a multiple of 3!
                            refDegap1 = ''
                            targetDegap1 = ''
                            for j in range(len(refAbridged2)):
                                if refAbridged2[j] != "-":
                                    refDegap1 += refAbridged2[j]
                                    targetDegap1 += targetAbridged2[j]

                            # removing gap-containing codons from alignment - these are all codons with a gap, or deletion,
                            # in the query sequence. Removal of entire codons should keep the alignment a multiple of 3
                            refDegap2 = ''
                            targetDegap2 = ''
                            originalTargetSeq = ''
                            for j in range(0, len(targetDegap1), 3):
                                targetCodon = (targetDegap1[j:j + 3])
                                refCodon = (refDegap1[j:j + 3])
                                if not re.findall(r'-',
                                                  targetCodon):  # this line is redundant and should do nothing
                                    if targetCodon not in ["TAG", 'TAA', "TGA"] and refCodon not in ["TAG", 'TAA',
                                                                                                     "TGA"]:
                                        targetDegap2 += targetCodon
                                        refDegap2 += refCodon
                                        originalTargetSeq += targetCodon

                            outffn = open("%s/aln%s.ffn" % (alnDir, count), "w")
                            outffn.write(">" + ref + "\n")
                            outffn.write(refDegap2 + "\n")
                            outffn.write(">" + target + "\n")
                            outffn.write(targetDegap2 + "\n")
                            outffn.close()

                            refAbridged3 = degapper(refAbridged2, targetAbridged2)[0]
                            targetAbridged3 = degapper(refAbridged2, targetAbridged2)[1]
                            refStraightDegap = remove(refAbridged3, ["-"])
                            targetStraightDegap = remove(targetAbridged3, ["-"])

                            refNoMercy = ""
                            targetNoMercy = ""
                            for j in range(0, sorted([len(refStraightDegap), len(targetStraightDegap)])[0], 3):
                                targetCodon = (targetStraightDegap[j:j + 3])
                                refCodon = (refStraightDegap[j:j + 3])
                                if targetCodon not in ["TAG", 'TAA', "TGA"] and refCodon not in ["TAG", 'TAA',
                                                                                                 "TGA"]:
                                    if len(targetCodon) == 3 and len(refCodon) == 3:
                                        targetNoMercy += targetCodon
                                        refNoMercy += refCodon

                            refNoMercyTrans = ribosome(refNoMercy)
                            targetNoMercyTrans = ribosome(targetNoMercy)

                            outFfnNoMercy = open("%s/aln2-%s.ffn" % (alnDir, count), "w")
                            outFfnNoMercy.write(">" + ref + "\n")
                            outFfnNoMercy.write(refNoMercy + "\n")
                            outFfnNoMercy.write(">" + target + "\n")
                            outFfnNoMercy.write(targetNoMercy + "\n")
                            outFfnNoMercy.close()

                            outFaaNoMercy = open("%s/aln2-%s.faa" % (alnDir, count), "w")
                            outFaaNoMercy.write(">" + ref + "\n")
                            outFaaNoMercy.write(refNoMercyTrans + "\n")
                            outFaaNoMercy.write(">" + target + "\n")
                            outFaaNoMercy.write(targetNoMercyTrans + "\n")
                            outFaaNoMercy.close()

                            stops = ribosome(targetStraightDegapRaw).count("*")
                            fragments = ribosome(targetStraightDegap).split("*")
                            firstFragment = fragments[0]
                            severity = len(fragments[0]) / len(ribosome(targetStraightDegap))

                            if stopcodon == "y":
                                stops = stops - 1

                            summaryDict[i]["numInframeInserts"] = numInframeInserts
                            summaryDict[i]["numOutframeInserts"] = numOutframeInserts
                            # summaryDict[i]["totalInserts"] = totalInserts / len(refDegap1)
                            summaryDict[i]["totalInserts"] = totalInserts
                            summaryDict[i]["numInframeDels"] = numInframeDels
                            summaryDict[i]["numOutframeDels"] = numOutframeDels
                            # summaryDict[i]["totalDels"] = totalDels / len(refDegap1)
                            summaryDict[i]["totalDels"] = totalDels
                            summaryDict[i]["startcodon"] = startcodon
                            summaryDict[i]["stopcodon"] = stopcodon
                            summaryDict[i]["stops"] = stops
                            summaryDict[i]["refAbridged2"] = refAbridged2
                            summaryDict[i]["refDegap2"] = refDegap2
                            summaryDict[i]["refNoMercy"] = refNoMercy
                            summaryDict[i]["targetAbridged2"] = targetAbridged2
                            summaryDict[i]["targetDegap2"] = targetDegap2
                            summaryDict[i]["targetNoMercy"] = targetNoMercy
                            summaryDict[i]["lengthDiff"] = lengthDiff
                            summaryDict[i]["stopSeverity"] = severity
                            summaryDict[i]["preferredStartLoss"] = preferredStartLoss
                            summaryDict[i]["preferredStartGain"] = preferredStartGain

    # # RUNNING CODEML
    count = 0
    for file in (os.listdir(out + "/nuc_aln")):
        if re.findall(r'aln2-', file):
            if lastItem(file.split(".")) == "ffn":
                count += 1
    total = count
    count = 0
    aaiDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    for file in os.listdir(alnDir):
        if re.findall(r'aln2-', file):
            if lastItem(file.split(".")) == "ffn":
                os.system("muscle -in %s/%s.faa -out %s/%s.faa.fa > /dev/null 2>&1" %
                          (alnDir, file.split(".")[0], alnDir, file.split(".")[0]))
                protAln = open("%s/%s.faa.fa" % (alnDir, file.split(".")[0]))
                protAln = fasta(protAln)
                ref = (list(protAln.keys())[0])
                refProt = protAln[ref]
                target = (list(protAln.keys())[1])
                targetProt = protAln[target]
                aai = AAI(refProt, targetProt)
                aaiDict[ref + "__" + target + ".fa"] = str(aai)
                os.system("pal2nal.pl %s/%s.faa.fa %s/%s -output fasta > %s/%s.ca.ffn" %
                          (alnDir, file.split(".")[0], alnDir, file, alnDir, file.split(".")[0]))
                count += 1
                perc = (count / total) * 100
                #sys.stdout.write("running pal2nal: %d%%   \r" % (perc))
                sys.stdout.write(f"{common.current_time()}\tRunning PAL2NAL: {round(perc)}%\r")
                sys.stdout.flush()
    print("")
    count = 0
    for file in (os.listdir(out + "/nuc_aln")):
        if lastItem(file.split(".")) == "ctl":
            count += 1
    total = count
    count = 0
    for file in os.listdir(alnDir):
        if lastItem(file.split(".")) == "ctl":
            os.system("codeml %s/%s > /dev/null 2>&1" % (alnDir, file))
            count += 1
            perc = (count / total) * 100
            #sys.stdout.write("running codeml: %d%%   \r" % (perc))
            sys.stdout.write(f"{common.current_time()}\tRunning CodeML: {round(perc)}%\r")
            sys.stdout.flush()
    print("")
    common.print_with_time("Done with CodeML.")
    time.sleep(2)
    count = 0
    alnDir = "%s/nuc_aln" % out
    for i in os.listdir(alnDir):
        if re.findall(r'mlcTree_', i):
            count += 1
    common.print_with_time("%s out of %s input reference CDS had detectable homology (below evalue of %s)" %
          (str(count), str(len(ffn.keys())), str(e)))
    common.print_with_time("Writing output file...")

    # PARSING CODEML OUTPUT AND COMBINING WITH OTHER RESULTS
    fastaOut = open(out + "/ref_based_cds_predictions.ffn", "w")
    mainOut = open(out + "/sleuth_report.csv", "w")
    mainOut.write(",".join(data_structures.SleuthData._fields) + "\n")
    # mainOut.write("reference_locus,target_locus,ANI,AAI,ref_length,aln_query_cov,start,"
    #
    #               "loss_of_preferred_start,gain_of_preferred_start_codon,stop_codon,internal_stops,first_stop_codon,"
    #
    #               "out_of_frame_inserts,out_of_frame_dels,inframe_inserts,inframe_dels,total_inserts,total_deleted_bases,"
    #
    #               "ds,ds_no_mercy,dsds,delta_ds,dnds,dnds_no_mercy,"
    #
    #               "full_seq,mercy_aln,no_mercy_aln,full_ref_seq,mercy_aln_ref,no_mercy_aln_ref\n")

    for i in os.listdir(alnDir):
        if re.findall(r'mlcTree_', i):
            originalname = i.split("mlcTree_")[1]
            file = open("%s/%s" % (alnDir, originalname))
            file = fasta2(file)

            ref = (list(file.keys())[0])
            refSeq = file[ref]
            target = (list(file.keys())[1])
            targetSeq = file[target]

            ANI = aniDict[ref][target]

            #####################################################
            #####################################################
            mlc = open("%s/%s" % (alnDir, i))
            for line in mlc:
                if not re.match(r'^ ', line):
                    try:
                        dn = (line.rstrip().split("dN = ")[1])
                        ds = (line.rstrip().split("dS = ")[1])
                    except IndexError:
                        pass

            ######################################
            try:
                dn = float(dn.split(" ")[0])
                ds = float(ds.split(" ")[0])
            except AttributeError:
                dnds = "NA"

            ############################################
            if ds >= 0.001:
                dnds = dn / ds
            else:
                dnds = "NA"

            ####################################################
            ####################################################
            mlcNoMercy = open("%s/mlcTree2_%s" % (alnDir, i.split("Tree_")[1]))
            for line in mlcNoMercy:
                if not re.match(r'^ ', line):
                    try:
                        # print(line.rstrip())
                        dnNoMercy = (line.rstrip().split("dN = ")[1])
                        dsNoMercy = (line.rstrip().split("dS = ")[1])
                    except IndexError:
                        pass

            #####################################
            try:
                dnNoMercy = float(dnNoMercy.split(" ")[0])
                dsNoMercy = float(dsNoMercy.split(" ")[0])
            except AttributeError:
                dndsNoMercy = "NA"

            #########################################
            if dsNoMercy >= 0.001:
                try:
                    dndsNoMercy = dnNoMercy / dsNoMercy
                except ZeroDivisionError:
                    dndsNoMercy = "NA"
            else:
                dndsNoMercy = "NA"

            #########################################
            if ds > 0.001:
                dsds = dsNoMercy / ds
                deltads = dsNoMercy - ds
            else:
                dsds = "NA"
                deltads = dsNoMercy - ds

            numInframeInserts = summaryDict[originalname]["numInframeInserts"]
            numOutframeInserts = summaryDict[originalname]["numOutframeInserts"]
            totalInserts = summaryDict[originalname]["totalInserts"]
            numInframeDels = summaryDict[originalname]["numInframeDels"]
            numOutframeDels = summaryDict[originalname]["numOutframeDels"]
            totalDels = summaryDict[originalname]["totalDels"]
            startcodon = summaryDict[originalname]["startcodon"]
            stopcodon = summaryDict[originalname]["stopcodon"]
            stops = summaryDict[originalname]["stops"]
            refAbridged2 = summaryDict[originalname]["refAbridged2"]
            refDegap2 = summaryDict[originalname]["refDegap2"]
            refNoMercy = summaryDict[originalname]["refNoMercy"]
            targetAbridged2 = summaryDict[originalname]["targetAbridged2"]
            targetDegap2 = summaryDict[originalname]["targetDegap2"]
            targetNoMercy = summaryDict[originalname]["targetNoMercy"]
            lengthDiff = summaryDict[originalname]["lengthDiff"]
            severity = summaryDict[originalname]["stopSeverity"]
            preferredStartLoss = summaryDict[originalname]["preferredStartLoss"]
            preferredStartGain = summaryDict[originalname]["preferredStartGain"]

            aai = aaiDict[originalname]
            targetStraightDegapRaw = remove(targetAbridged2, ["-"])

            fastaOut.write(">" + target + "\n")
            fastaOut.write(targetStraightDegapRaw + "\n")
            mainOut.write(ref + "," + target + "," + str(ANI) + "," + str(aai) + "," + str(len(refAbridged2)) + "," + str(lengthDiff) + "," +
                          startcodon + "," + preferredStartLoss + "," + preferredStartGain + "," + stopcodon + "," + str(
                stops) + "," + str(severity) + "," +
                          str(numOutframeInserts) + "," + str(numOutframeDels) + "," + str(
                numInframeInserts) + "," +
                          str(numInframeDels) + "," + str(totalInserts) + "," + str(totalDels) + "," +
                          str(ds) + "," + str(dsNoMercy) + "," + str(dsds) + "," + str(deltads) + "," + str(dnds) + "," + str(dndsNoMercy) + "," + targetAbridged2 +
                          "," + targetDegap2 + "," + targetNoMercy + "," + refAbridged2 + "," + refDegap2 + "," + refNoMercy + "\n")

    mainOut.close()
    fastaOut.close()
    common.print_with_time("Sleuth finished.")

    os.system("rm -f 2NG.t 2NG.dN 2NG.dS rst1 rst 2ML.t 2ML.dN 2ML.dS 4fold.nuc rub")
    os.system("tar -cf %s/nuc_aln.tar %s/nuc_aln" % (out, out))
    os.system("gzip %s/nuc_aln.tar" % out)


def relate_sleuth_data_to_locus_tags(args, file_dict: dict, log_file_dict=None) -> dict:
    """
    Compares sleuth_data.full_seq to file_dict['cds_filename'] via blast, creating a link between the two data types.
    Returns a dictionary, where the keys are locus tags and the values are annotate.sleuth_data NamedTuples.
    """
    if log_file_dict is None:
        cds_annotations = file_dict['cds_filename']
    else:
        cds_annotations = log_file_dict['cds_filename']

    sleuth_data = sleuth_report_reader(file_dict['sleuth_report'])
    cds_predictions = file_dict['cds_predictions']
    db_name = file_dict['temp_dir'] + '/cds_db'
    blast_out = file_dict['sleuthDir'] + '/sleuth.finder.blast'

    makeblastdb_dict = {'dbtype': 'nucl',
                        'input_file': cds_annotations,
                        'out': db_name}

    blastn_dict = {'query': cds_predictions,
                   'db': db_name,
                   'outfmt': 6,
                   'out': blast_out,
                   'num_threads': args.threads,
                   'perc_identity': 99.9,
                   'evalue': args.evalue}

    NcbimakeblastdbCommandline(**makeblastdb_dict)()
    NcbiblastnCommandline(**blastn_dict)()

    blastdict = dict()
    with open(blast_out, 'r') as blast_out:
        for line in blast_out.readlines():
            fields = line.split("\t")
            query, ref = fields[0], fields[1]
            blastdict[ref] = sleuth_data[query]

    return blastdict


def sleuth_report_reader(sleuth_report: str) -> dict:
    """
    Reads in the sleuth report csv file and returns a dictionary, where the keys are sleuth_data.target_locus and the values are
    annotate.sleuth_data NamedTuples.
    """
    sleuthDict = dict()
    with open(sleuth_report, 'r') as sleuth_report:
        for line in sleuth_report.readlines()[1:]:  # Read each line, skipping the first (header line)
            fields = [common.literal_eval(x.rstrip()) for x in line.split(sep=',')]  # Convert each field from string to appropriate atomic data type
            data = data_structures.SleuthData(*fields)    # Unpack the list of fields into the sleuth_data NamedTuple
            sleuthDict[data.target_locus] = data    # Add to dictionary

    return sleuthDict


def merge(args, file_dict, log_file_dict=None):
    # CONSOLIDATING PSEUDOFINDER AND SLEUTH OUTPUT
    out = file_dict["sleuthDir"]
    cds = file_dict["cds_filename"]
    e = args.evalue
    base = file_dict["base_filename"]

    os.system("makeblastdb -dbtype nucl -in %s -out %s" % (cds, cds))
    os.system("blastn -db %s -query %s/ref_based_cds_predictions.ffn -outfmt 6 -out %s/sleuth.finder.blast -num_threads 4 -perc_identity 99.9 -evalue %s" % (cds, out, out, e))

    blastDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    blast = open("%s/sleuth.finder.blast" % out)
    for i in blast:
        ls = i.rstrip().split("\t")
        blastDict[ls[1]] = ls[0]

    pseudoDict = defaultdict(list)
    intactDict = defaultdict(list)
    pseudoSeqDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    summary = open("%s/sleuth_report.csv" % out)
    for i in summary:
        ls = i.rstrip().split(",")
        if ls[2] != "ANI":
            aai = float(ls[3])
            cov = float(ls[5])
            start = (ls[6])
            stop = (ls[9])
            instop = int(ls[10])
            stopSeverity = float(ls[11])
            fsIn = int(ls[12])
            fsDel = int(ls[13])

            #dN/dS
            if ls[22] != "NA":
                dnds = float(ls[22])
                if dnds > 0.3:
                    if dnds < 3:
                        string = "Elevated dN/dS: " + str(round(dnds, 4)) + "."
                        pseudoDict[ls[1]].append(string)
                    else:
                        string = "Elevated dN/dS: " + str(round(dnds, 4)) + " (this exceptionally high dN/dS is likely " \
                                                                            "caused by a poor alignment, which may be " \
                                                                            "indiciative of a true pseudogene, or a false " \
                                                                            "positive BLAST hit)."
                        pseudoDict[ls[1]].append(string)

            # dS/dS
            if ls[20] != "NA":
                dsds = float(ls[20])
                if dsds > 15:
                    string = str(round(fsIn + fsDel)) + " significant frameshift-inducing indel(s)."
                    pseudoDict[ls[1]].append(string)
            else:
                deltaDS = float(ls[21])
                if deltaDS > 0.1:
                    string = str(round(fsIn + fsDel)) + " significant frameshift-inducing indel(s)."
                    pseudoDict[ls[1]].append(string)

            # START
            start = ls[6]
            if start == "n":
                string = "Missing start codon."
                pseudoDict[ls[1]].append(string)

            # STOP
            stop = ls[9]
            internal = int(ls[10])
            position = float(ls[11])
            if stop == "n" and internal == 0:
                string = "Missing stop codon."
                pseudoDict[ls[1]].append(string)

            if internal != 0 and position < 0.75:
                string = "Internal stop codon at " + str(round(position * 100)) + "% expected length."
                pseudoDict[ls[1]].append(string)

            pseudoSeqDict[ls[1]] = ls[24]


            # ds = float(ls[16])
            # dsNoMercy = float(ls[18])
            # dnds = (ls[17])
            # try:
            #     bias = sorted([ds, dsNoMercy])[1] / sorted([ds, dsNoMercy])[0]
            # except ZeroDivisionError:
            #     bias = sorted([ds, dsNoMercy])[1]
            #
            # if fsIn + fsDel > 0:
            #     if bias > 1.25:
            #         string = str(round(fsIn + fsDel)) + " significant frameshift-inducing indel(s)."
            #         pseudoDict[ls[1]].append(string)
            #
            # if instop > 0:
            #     if stopSeverity < args.length_pseudo:
            #         string = "Internal stop codon at " + str(round(stopSeverity * 100)) + "% expected length."
            #         pseudoDict[ls[1]].append(string)
            #
            # if dnds != "NA":
            #     dnds = float(dnds)
            #     if dnds > args.max_dnds:
            #         if dnds < 3:
            #             string = "Elevated dN/dS: " + str(round(dnds, 4)) + "."
            #             pseudoDict[ls[1]].append(string)
            #         else:
            #             string = "Elevated dN/dS: " + str(round(dnds, 4)) + " (this exceptionally high dN/dS is likely caused by a poor alignment, which may be indiciative of a true pseudogene, or a false positive BLAST hit)."
            #             pseudoDict[ls[1]].append(string)
            #
            # if bias < 1.10:
            #     if stopSeverity > 0.95:
            #         if cov > args.length_pseudo:
            #             if dnds != "NA":
            #                 dnds = float(dnds)
            #                 if dnds < 0.1:
            #                     intactDict[ls[1]] = ls
            #             else:
            #                 intactDict[ls[1]] = ls
            # pseudoSeqDict[ls[1]] = ls[20]

    count = 0
    cds = open(file_dict["cds_filename"])
    cds = fasta2(cds)
    for i in cds:
        if i in blastDict.keys():
            count += 1
    # print(count)
    # print(len(cds.keys()))

    count = 0
    intactToPseudo = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    intact = open(file_dict["intact_gff"])
    outNewIntact = open(file_dict["new_intact_gff"], 'w')
    for i in intact:
        if not re.match(r'#', i):
            ls = i.rstrip().split("\t")
            locus = ls[8].split("=")[1]
            if locus in blastDict.keys():
                pseudoLocus = blastDict[locus]
                if pseudoLocus in pseudoDict.keys():
                    count += 1
                    intactToPseudo[locus]["reason"] = pseudoDict[pseudoLocus]
                    intactToPseudo[locus]["stats"] = ls
                else:
                    outNewIntact.write(i.rstrip() + "\n")
        else:
            outNewIntact.write(i.rstrip() + "\n")

    # print(count)

    count = 0
    counter = 0
    # REMOVING GENES THAT LOOK REALLY INTACT FROM PSEUDOS.GFF
    # ADDING MORE REASONS FOUND FOR CONFIRMED PSEUDOGENES TO PSEUDOS.GFF
    PseudoToPseudo = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    PseudoToIntact = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    locusDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    dndsDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    redundList = []
    pseudos = open(file_dict["pseudos_gff"])
    outNewPseudo = open(file_dict["new_pseudos_gff"], "w")
    for i in pseudos:
        if not re.match(r'#', i):
            ls = i.rstrip().split("\t")
            oldLocus = (ls[8].split("old_locus_tag=")[1])
            newLocus = (ls[8].split("locus_tag=")[1].split(";")[0])
            counter += 1
            locusDict[newLocus] = oldLocus
            locusDict[oldLocus] = newLocus
            if not re.findall(r'dN/dS', i):
                if not re.findall(r'ign', oldLocus):
                    if len(oldLocus.split(",")) > 1:
                        for j in oldLocus.split(","):
                            locus = j
                            if locus in blastDict.keys():
                                pseudoLocus = blastDict[locus]
                                if pseudoLocus in pseudoDict.keys():
                                    # PseudoToPseudo[locus] = pseudoDict[pseudoLocus]
                                    reason = (ls[8].split(";")[0].split("Reason: ")[1]) + " "
                                    for k in pseudoDict[pseudoLocus]:
                                        reason += k + " "
                                    newReason = ("note=Pseudogene candidate. Reason: " + reason + "".join(
                                        allButTheFirst(ls[8], ";")))
                                    count += 1
                                    outNewPseudo.write(
                                        ls[0] + "\t" + ls[1] + "\t" + ls[2] + "\t" + ls[3] + "\t" + ls[4] + "\t" + ls[
                                            5] + "\t" + ls[6] + "\t" + ls[7] + "\t" + str(newReason) + "\n")
                                    break

                                elif pseudoLocus in intactDict.keys():
                                    PseudoToIntact[locus] = ls

                                else:
                                    if oldLocus not in redundList:
                                        redundList.append(oldLocus)
                                        outNewPseudo.write(i.rstrip() + "\n")
                                        count += 1

                            else:
                                if oldLocus not in redundList:
                                    redundList.append(oldLocus)
                                    outNewPseudo.write(i.rstrip() + "\n")
                                    count += 1

                    else:
                        locus = oldLocus
                        if locus in blastDict.keys():
                            pseudoLocus = blastDict[locus]
                            if pseudoLocus in pseudoDict.keys():
                                # PseudoToPseudo[locus] = pseudoDict[pseudoLocus]
                                reason = (ls[8].split(";")[0].split("Reason: ")[1]) + " "
                                for k in pseudoDict[pseudoLocus]:
                                    reason += k + " "
                                newReason = (
                                "note=Pseudogene candidate. Reason: " + reason + "".join(allButTheFirst(ls[8], ";")))
                                count += 1
                                outNewPseudo.write(
                                    ls[0] + "\t" + ls[1] + "\t" + ls[2] + "\t" + ls[3] + "\t" + ls[4] + "\t" + ls[
                                        5] + "\t" + ls[6] + "\t" + ls[7] + "\t" + str(newReason) + "\n")

                            elif pseudoLocus in intactDict.keys():
                                PseudoToIntact[locus] = ls

                            else:
                                count += 1
                                outNewPseudo.write(i.rstrip() + "\n")

                        else:
                            count += 1
                            outNewPseudo.write(i.rstrip() + "\n")

                else:
                    outNewPseudo.write(i.rstrip() + "\n")
                    count += 1
            else:
                if oldLocus in blastDict.keys():
                    pseudoLocus = blastDict[oldLocus]
                    if pseudoLocus in pseudoDict.keys():
                        reason = ""
                        count += 1
                        for k in pseudoDict[pseudoLocus]:
                            reason += k + " "
                        newReason = (
                        "note=Pseudogene candidate. Reason: " + reason + "".join(allButTheFirst(ls[8], ";")))
                        outNewPseudo.write(
                            ls[0] + "\t" + ls[1] + "\t" + ls[2] + "\t" + ls[3] + "\t" + ls[4] + "\t" + ls[
                                5] + "\t" + ls[6] + "\t" + ls[7] + "\t" + str(newReason) + "\n")

                    else:
                        dndsDict[newLocus] = ls
                else:
                    dndsDict[newLocus] = ls

        else:
            outNewPseudo.write(i.rstrip() + "\n")

    # REMOVING GENES THAT LOOK REALLY INTACT FROM PSEUDOS.FASTA
    # REMEMBERING THE LAST INDEX NUMBER FOR NEW PSEUDOGENES (USED IN NEXT CODE BLOCK TO GENERATE NEW PSEUDOGENE NAMES)
    pseudoSeqs = open(file_dict["pseudos_fasta"])
    pseudoSeqs = fasta(pseudoSeqs)
    pseudoSeqsOut = open(file_dict["new_pseudos_fasta"], "w")
    num = 0
    for i in pseudoSeqs.keys():
        newLocus = i.split(" ")[0]
        oldLocus = locusDict[i.split(" ")[0]]
        if oldLocus in PseudoToIntact.keys():
            pass
        elif newLocus in dndsDict.keys():
            pass
        else:
            num = (i.split(" ")[0].split("_")[3])
            pseudoSeqsOut.write(">" + i + "\n")
            pseudoSeqsOut.write(pseudoSeqs[i] + "\n")

    # ADDING NEWLY DISCOVERED PSEUDOGENES TO PSEUDOS.GFF
    # ADDING NEWLY DISCOVERED PSEUDOGENES TO PSEUDOS.FASTA
    num = int(num)
    for i in intactToPseudo.keys():
        stats = intactToPseudo[i]["stats"]
        reason = intactToPseudo[i]["reason"]
        newReason = ""
        for j in reason:
            newReason += j + " "
        num += 1
        newLocus = stats[0] + "_pseudo_" + str(stabilityCounter(int(num)))
        newLine = stats[0] + "\t" + stats[1] + "\t" + stats[2] + "\t" + stats[3] + "\t" + stats[4] + "\t" + stats[5] + \
                  "\t" + stats[6] + "\t" + stats[7] + "\t" + "note=Pseudogene candidate. Reason: " + newReason[0:len(
            newReason) - 1] + "; locus_tag=" + newLocus + "; old_locus_tag=" + i
        outNewPseudo.write(newLine + "\n")
        count += 1
        pseudoSeqsOut.write(
            ">" + newLocus + " " + stats[0] + "[" + stats[3] + ":" + stats[4] + "](" + stats[6] + ")" + "\n")
        pseudoSeqsOut.write(pseudoSeqDict[blastDict[i]] + "\n")

    pseudoSeqsOut.close()
    outNewPseudo.close()
    print("\n\n\n\n")
    for i in PseudoToIntact.keys():
        ls = (PseudoToIntact[i])
        oldLocus = (ls[8].split("old_locus_tag=")[1])
        outNewIntact.write(
            ls[0] + "\t" + ls[1] + "\t" + ls[2] + "\t" + ls[3] + "\t" + ls[4] + "\t" + ls[5] + "\t" + ls[6] + "\t" + ls[
                7] + "\tlocus_tag=" + oldLocus + '\n')
    outNewIntact.close()

    CDS = open(file_dict["cds_filename"])
    CDS = fasta(CDS)
    prots = open(file_dict["proteome_filename"])
    prots = fasta(prots)
    outNewIntactProts = open(file_dict["new_intact_faa"], 'w')
    outNewIntactSeqs = open(file_dict["new_intact_ffn"], 'w')

    intactGffDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    intactGff = open(file_dict["new_intact_gff"])
    for i in intactGff:
        if not re.match(r'#', i):
            ls = i.rstrip().split("\t")
            locus = ls[8].split("=")[1]
            intactGffDict[locus] = ls

    for i in CDS.keys():
        if (i.split(" ")[0]) in intactGffDict.keys():
            outNewIntactSeqs.write(">" + i + "\n")
            outNewIntactSeqs.write(CDS[i] + "\n")
    outNewIntactSeqs.close()

    for i in prots.keys():
        if (i.split(" ")[0]) in intactGffDict.keys():
            outNewIntactProts.write(">" + i + "\n")
            outNewIntactProts.write(prots[i] + "\n")

    outNewIntactProts.close()

    blastDict = defaultdict(list)
    blast = open("%s/sleuth.finder.blast" % out)
    for i in blast:
        ls = i.rstrip().split("\t")
        blastDict[ls[0]].append(ls[1])

    blastDict2 = defaultdict(lambda: defaultdict(lambda: 'NA'))
    for i in blastDict.keys():
        LS = []
        for j in blastDict[i]:
            LS.append(j)
        loci = ";".join(LS)
        blastDict2[i] = loci

    summary = open("%s/sleuth_report.csv" % out)
    outfile = open("%s/sleuth_report-2.csv" % out, "w")
    for i in summary:
        ls = i.rstrip().split(",")
        if ls[0] != "reference_locus":
            try:
                bias = float(sorted([ls[16], ls[18]])[1]) / float(sorted([ls[16], ls[18]])[0])
            except ZeroDivisionError:
                bias = "NA"

            outfile.write(
                ls[0] + "," + str(len(ls[23])) + "," + ls[1] + "," + str(blastDict2[ls[1]]) + "," + ls[2] + "," + ls[3] + "," + ls[4] + "," +
                ls[5] + "," + ls[6] + "," + ls[7] + "," + ls[8] + "," + ls[9] + "," + ls[10] + "," +
                ls[11] + "," + ls[12] + "," + ls[13] + "," + ls[14] + "," + ls[15] + "," + ls[16] + "," +
                ls[17] + "," + ls[18] + "," + ls[19] + "," + str(bias) + "," + ls[20] + "," + ls[21] + "," + ls[22] + "," +
                ls[23] + "," + ls[24] + "," + ls[25] + "\n")
        else:
            outfile.write(ls[0] + "," + "ref_length" + "," + ls[1] + ",target_gbk_locus_tags," + ls[2] + "," + ls[3] + "," + ls[4] + "," +
                          ls[5] + "," + ls[6] + "," + ls[7] + "," + ls[8] + "," + ls[9] + "," + ls[10] + "," +
                          ls[11] + "," + ls[12] + "," + ls[13] + "," + ls[14] + "," + ls[15] + "," + ls[16] + "," +
                          ls[17] + "," + ls[18] + "," + ls[19] + ",ds_ratio," + "," + ls[20] + "," + ls[21] + "," + ls[22] + "," +
                          ls[23] + "," + ls[24] + "," + ls[25] + "\n")

    os.system("mv %s/sleuth_report.csv ./" % out)
    os.system("mv %s/sleuth_report-2.csv ./sleuth_report.csv" % out)
    os.system("tar -cf %s.tar %s" % (out, out))
    os.system("gzip %s.tar" % out)
    os.system("mv %s %s" % (file_dict["new_intact_faa"], file_dict["intact_faa"]))
    os.system("mv %s %s" % (file_dict["new_intact_ffn"], file_dict["intact_ffn"]))
    os.system("mv %s %s" % (file_dict["new_intact_gff"], file_dict["intact_gff"]))
    os.system("mv %s %s" % (file_dict["new_pseudos_gff"], file_dict["pseudos_gff"]))
    os.system("mv %s %s" % (file_dict["new_pseudos_fasta"], file_dict["pseudos_fasta"]))
    os.system("rm *.fasta.n*")
    os.system("rm -r %s" % out)
    os.system("mv sleuth_report.csv %ssleuth_report.csv" % base)









