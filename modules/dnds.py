#!/usr/bin/env python3
from collections import defaultdict
import re
import os
import textwrap
import argparse
import numpy as np
import sys


def localize(item, ls):
    count = 0
    for i in ls:
        if i == item:
            break
        else:
            count += 1
    return count


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


def reject_outliers(data):
    m = 2
    u = np.mean(data)
    s = np.std(data)
    filtered = [e for e in data if (u - 2 * s < e < u + 2 * s)]
    return filtered


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
    parser = argparse.ArgumentParser(
        prog="dnds.py",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent('''
        ************************************************************************
    
        Developed by Arkadiy Garber and John McCutcheon
        University of Montana, Biological Sciences
        Please send comments and inquiries to arkadiy.garber@mso.umt.edu
    
        ************************************************************************
        '''))


    parser.add_argument('-n', type=str, help="ORFs in nucleic acid format", default="NA")
    parser.add_argument('-a', type=str, help="ORFs in amino acid format", default="NA")
    parser.add_argument('-rn', type=str, help="Reference ORFs in nucleic acid format", default="NA")
    parser.add_argument('-ra', type=str, help="Reference ORFs in amino acid format", default="NA")
    parser.add_argument('-r', type=str, help="Reference contigs (provide this if you do not have ORFs for the reference genome)",
                        default="NA")

    parser.add_argument('-ctl', type=str, help="template control file for codeml", default="NA")

    parser.add_argument('-out', type=str, help="name output directory", default="PseudoHunter_output")

    parser.add_argument('-d', type=float, help="maximum dN/dS value for gene too be considered \'intact\' (default = 0.3)",
                        default=0.3)
    parser.add_argument('-M', type=float, help="maximum dS value for dN/dS calculation (default = 3)", default=3)
    parser.add_argument('-m', type=float, help="minimum dS value for dN/dS calculation (default = 0.001)", default=0.001)
    parser.add_argument('-t', type=int, help="number of threads to use for BLAST or DIAMOND", default=1)
    parser.add_argument('-s', type=str, help="search engine to use (blast/diamond). Default = blast", default="blast")


    args = parser.parse_args()

    cwd = os.getcwd()

    os.system("echo ${ctl} > ctl.txt")
    file = open("ctl.txt")
    for i in file:
        ctl = (i.rstrip())
    os.system("rm ctl.txt")

    if not re.findall(r'codeml', ctl):
        if args.ctl == "NA":
            print("\nCannot locate the necessary \ncodeml control file. Please double check that "
                  "installation\n with the \'setup.sh\' script worked without any errors. If issues \npersist, please "
                  "report an issue on GitHub\n")
            raise SystemExit
        else:
            ctl = args.ctl

    faa = args.a
    fna = args.n

    if args.ra != "NA" and args.rn != "NA":
        refFaa = args.ra
        refFna = args.rn
    else:
        if args.r != "NA":

            os.system("prodigal -i %s -a %s-proteins.faa -d %s-proteins.fna" % (args.r, args.r, args.r))

            refFna = args.r + "-proteins.fna"
            refFaa = args.r + "-proteins.faa"
        else:
            print("Did not find reference datasets. Please provided these using the \'-ra\' and \'rn\', or \'r\' flag")


    print("Starting pipeline...")

    os.system("mkdir " + args.out)
    os.system('mkdir ' + args.out + "/dnds-analysis")

    if args.s == "blast":
        print("Running BLAST")
        os.system("makeblastdb -dbtype prot -in %s -out %s" % (refFaa, refFaa))
        # os.system("rm makeblastdb.out")
        os.system("blastp -query %s -db %s "
                  "-outfmt 6 -out %s/pseudogene.blast -evalue 1E-6 -num_threads %s -max_target_seqs 1" % (
                      faa, refFaa, args.out, args.t))

        os.system("rm %s.psq" % refFaa)
        os.system("rm %s.phr" % refFaa)
        os.system("rm %s.pin" % refFaa)

    elif args.s == "diamond":
        print("Running DIAMOND")
        os.system(
            "diamond makedb --in %s -d %s" % (refFaa, refFaa))

        os.system("diamond blastp --db %s.dmnd --query %s --outfmt 6 --out %s/pseudogene.blast "
                  "--max-target-seqs 1 --evalue 1E-6 --threads %d" % (refFaa, faa, args.out, args.t))


    ####################################################################################################################
    faaDict = open(faa)
    faaDict = fasta2(faaDict)

    fnaDict = open(fna)
    fnaDict = fasta2(fnaDict)

    refFaaDict = open(refFaa)
    refFaaDict = fasta2(refFaaDict)

    refFnaDict = open(refFna)
    refFnaDict = fasta2(refFnaDict)

    prescreened = []
    alnLengthDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    blast = open("%s/pseudogene.blast" % args.out)
    for i in blast:
        ls = i.rstrip().split("\t")
        if ls[0] not in prescreened:
            alnLengthDict[ls[0]] = ls[3]

            outNUC = open(args.out + "/dnds-analysis/%s.faa.fna" % ls[0], "w")
            outNUC.write(">" + ls[1] + "\n")
            outNUC.write(refFnaDict[ls[1]] + "\n")
            outNUC.write(">" + ls[0] + "\n")
            outNUC.write(fnaDict[ls[0]] + "\n")
            outNUC.close()

            outAA = open(args.out + "/dnds-analysis/%s.faa" % ls[0], "w")
            outAA.write(">" + ls[1] + "\n")
            outAA.write(refFaaDict[ls[1]] + "\n")
            outAA.write(">" + ls[0] + "\n")
            outAA.write(faaDict[ls[0]] + "\n")
            outAA.close()


    # ALIGNING PROTEIN SEQUENCES AND CREATING A CODON ALIGNMENT
    print("aligning files...")
    DIR = args.out + "/dnds-analysis"
    os.system("for i in %s/*faa; do"
              " muscle -in $i -out $i.aligned.fa;"
              # " rm muscle.out;"
              " pal2nal.pl $i.aligned.fa $i.fna -output fasta > $i.codonalign.fa;"
              " done" % DIR)

    # BUILDING CONTROL FILES
    print("preparing for codeml analysis")
    DIR = args.out + "/dnds-analysis"
    codealign = os.listdir(DIR)
    count = 0
    for file in codealign:
        if re.findall(r'codonalign', file):
            count += 1
            clu = file.split(".faa")[0]
            setup = open(ctl)
            out = open("%s/%s.ctl" % (DIR, str(clu)), "w")

            for i in setup:
                if re.findall('seqfile', i):
                    out.write(
                        '' + ' ' + '' + ' ' + '' + ' ' + '' + ' ' + '' + ' ' + '' + ' ' + 'seqfile' + ' ' + '= ' + args.out + '/dnds-analysis/' + file + ' ' + '*' + ' ' + 'sequence' + ' ' + 'data' + ' ' + 'filename\n')

                elif re.findall(r'outfile', i):
                    out.write(
                        '' + ' ' + '' + ' ' + '' + ' ' + '' + ' ' + '' + ' ' + '' + ' ' + 'outfile' + ' ' + '=' + ' '
                        + args.out + '/dnds-analysis/mlcTree_' + str(
                            clu) + ' ' + '' + ' ' + '' + ' ' + '' + ' ' + '' + ' ' + '' + ' ' + '' + ' ' + '' + ' '
                        + '' + ' ' + '' + ' ' + '' + ' ' + '*' + ' ' + 'main' + ' ' + 'result' + ' ' + 'file' + ' ' + 'name\n')

                else:
                    out.write(i)
            out.close()
    print("")

    # RUNNING CODEML FOR DN/DS CALCULATION
    total = 0
    for i in codealign:
        if re.findall(r'codonalign', i):
            total += 1

    count = 0
    codealign = os.listdir(DIR)
    for file in codealign:
        if lastItem(file.split(".")) == "ctl":
            count += 1
            perc = (count / total) * 100
            sys.stdout.write("running codeml: %d%%   \r" % (perc))
            sys.stdout.flush()
            os.system("codeml %s/dnds-analysis/%s" % (args.out, file))
            # os.system("rm codeml.out")

    # PARSING CODEML OUTPUT

    cwd = os.getcwd()
    DIR = args.out + "/dnds-analysis"


    alnLengthDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    alnIdDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    blast = open(args.out + "/pseudogene.blast")
    for i in blast:
        ls = i.rstrip().split("\t")
        alnLengthDict[ls[0]][ls[1]] = ls[3]
        alnIdDict[ls[0]][ls[1]] = ls[2]

    print("summarizing codeml output")
    codealign = os.listdir(DIR)
    dndsDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    for i in codealign:
        if re.findall(r'mlc', i):
            file = open(DIR + "/%s" % i, "r")
            for j in file:
                if re.search(r'#1', j):
                    ls = (j.rstrip().split(" "))
                    orf = ls[1]
                if re.search(r'#2', j):
                    ls = (j.rstrip().split(" "))
                    NODE = ls[1]
                line = j.rstrip()
            ls = line.split("  ")
            dS = remove(lastItem(ls), [" ", "=", "d", "S"])
            dN = remove(lastItem(ls[0:len(ls) - 1]), [" ", "=", "d", "N"])
            dndsDict[NODE]["orf"] = orf
            dndsDict[NODE]["dn"] = dN
            dndsDict[NODE]["ds"] = dS

    count = 0
    dndsList = []
    dndsDict2 = defaultdict(list)
    for i in sorted(dndsDict.keys()):
        count += 1
        if float(dndsDict[i]["dn"]) <= args.M and float(dndsDict[i]["ds"]) <= args.M and float(
                dndsDict[i]["ds"]) >= args.m and float(dndsDict[i]["dn"]) >= args.m:
            orf = dndsDict[i]["orf"]
            dn = dndsDict[i]["dn"]
            ds = dndsDict[i]["ds"]
            dnds = float(dndsDict[i]["dn"]) / float(dndsDict[i]["ds"])
            dndsDict2[orf].append(i)
            dndsList.append(dnds)

    dsList = []
    dndsList = []
    total = 0
    out = open(args.out + "/dnds-summary.csv", "w")
    out.write("ORF" + "," + "referenceMatch" + "," + "dN" + "," + "dS" + "," + "dN/dS" + "," + "PG" + "\n")
    for i in sorted(dndsDict.keys()):
        total += 1

        if float(dndsDict[i]["ds"]) < float(args.M) and float(dndsDict[i]["ds"]) > float(args.m):

            dnds = float(dndsDict[i]["dn"]) / float(dndsDict[i]["ds"])

            dsList.append(float(dndsDict[i]["ds"]))
            dndsList.append(dnds)

            if dnds < args.d:
                pg = "N"
            else:
                pg = "Y"

            out.write(i + "," + str(dndsDict[i]["orf"]) + "," + str(dndsDict[i]["dn"]) + "," + str(dndsDict[i]["ds"]) + "," +

                      str(dnds) + "," + str(pg) + "\n")

    out.close()

    print("")
    print("Done!")
    print("********************************************************************************")
    print("Identified a total of %d orthologs between the query and reference datasets," % total)
    print("with %d orthologs between the dS range of %d-%d" % (len(dsList), args.m, args.M))
    print("")
    print("Average dS among orthologs within the specified dS range: " + str(ave(dsList)))
    print("Average dN/dS among orthologs within the specified dS range: " + str(ave(dndsList)))
    print("")
    print("See %s/dnds-summay.csv for detailed results." % args.out)
    print("********************************************************************************")
    print("Thanks for using!")

    os.system("rm -f %s/pseudogene.blast" % args.out)
    os.system("rm -f 2NG.t 2NG.dN 2NG.dS rst1 rst 2ML.t 2ML.dN 2ML.dS 4fold.nuc rub")







