#!/usr/bin/env python3
from collections import defaultdict
from . import common
import re
import os
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
    args = common.get_args('selection')


    if not args.skip:
        ctl = os.path.dirname(os.path.dirname(__file__)) + "/codeml-2.ctl"
        faa = args.prots
        fna = args.genes

        if args.ref_prots != "NA" and args.ref_genes != "NA":
            refFaa = args.ref_prots
            refFna = args.ref_genes
        else:
            if args.r != "NA":

                os.system("prodigal -i %s -a %s-proteins.faa -d %s-proteins.fna > /dev/null 2>&1" % (args.r, args.r, args.r))

                refFna = args.reference + "-proteins.fna"
                refFaa = args.reference + "-proteins.faa"
            else:
                print("Did not find reference datasets. Please provided these using the \'-ra\' and \'rn\', or \'r\' flag")


        print("Starting pipeline...")
        print(args.outdir)
        print(args.outdir + "/dnds-analysis")
        os.system("mkdir -p " + args.outdir)
        os.system("mkdir -p " + args.outdir + "/dnds-analysis")

        if args.search_engine == "blast":
            print("Running BLAST")
            os.system("makeblastdb -dbtype prot -in %s -out %s > /dev/null 2>&1" % (refFaa, refFaa))
            # os.system("rm makeblastdb.out")
            os.system("blastp -query %s -db %s "
                      "-outfmt 6 -out %s/pseudogene.blast -evalue 1E-6 -num_threads %s -max_target_seqs 1 > /dev/null 2>&1" % (
                          faa, refFaa, args.outdir, args.threads))

            os.system("rm %s.psq" % refFaa)
            os.system("rm %s.phr" % refFaa)
            os.system("rm %s.pin" % refFaa)

        elif args.search_engine == "diamond":
            print("Running DIAMOND")
            os.system(
                "diamond makedb --in %s -d %s > /dev/null 2>&1" % (refFaa, refFaa))

            os.system("diamond blastp --db %s.dmnd --query %s --outfmt 6 --out %s/pseudogene.blast "
                      "--max-target-seqs 1 --evalue 1E-6 --threads %s > /dev/null 2>&1" % (refFaa, faa, args.outdir, args.threads))


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
        aaiDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
        blast = open("%s/pseudogene.blast" % args.outdir)
        for i in blast:
            ls = i.rstrip().split("\t")
            if ls[0] not in prescreened:
                alnLengthDict[ls[0]] = ls[3]
                aaiDict[ls[0]] = ls[2]
                outNUC = open(args.outdir + "/dnds-analysis/%s.faa.fna" % ls[0], "w")
                outNUC.write(">" + ls[1] + "\n")
                outNUC.write(refFnaDict[ls[1]] + "\n")
                outNUC.write(">" + ls[0] + "\n")
                outNUC.write(fnaDict[ls[0]] + "\n")
                outNUC.close()

                outAA = open(args.outdir + "/dnds-analysis/%s.faa" % ls[0], "w")
                outAA.write(">" + ls[1] + "\n")
                outAA.write(refFaaDict[ls[1]] + "\n")
                outAA.write(">" + ls[0] + "\n")
                outAA.write(faaDict[ls[0]] + "\n")
                outAA.close()


        # ALIGNING PROTEIN SEQUENCES AND CREATING A CODON ALIGNMENT
        print("aligning files...")
        DIR = args.outdir + "/dnds-analysis"
        os.system("for i in %s/*faa; do"
                  " muscle -in $i -out $i.aligned.fa > /dev/null 2>&1;"
                  # " rm muscle.out;"
                  " pal2nal.pl $i.aligned.fa $i.fna -output fasta > $i.codonalign.fa;"
                  " done" % DIR)

        # BUILDING CONTROL FILES
        print("preparing for codeml analysis")
        DIR = args.outdir + "/dnds-analysis"
        codealign = os.listdir(DIR)
        count = 0
        for file in codealign:
            if re.findall(r'codonalign', file):
                count += 1
                clu = file.split(".faa")[0]
                setup = open(ctl)
                print("%s/%s.ctl" % (DIR, str(clu)))
                out = open("%s/%s.ctl" % (DIR, str(clu)), "w")

                for i in setup:
                    if re.findall('seqfile', i):
                        out.write(
                            '' + ' ' + '' + ' ' + '' + ' ' + '' + ' ' + '' + ' ' + '' + ' ' + 'seqfile' + ' ' + '= ' + args.outdir + '/dnds-analysis/' + file + ' ' + '*' + ' ' + 'sequence' + ' ' + 'data' + ' ' + 'filename\n')

                    elif re.findall(r'outfile', i):
                        out.write(
                            '' + ' ' + '' + ' ' + '' + ' ' + '' + ' ' + '' + ' ' + '' + ' ' + 'outfile' + ' ' + '=' + ' '
                            + args.outdir + '/dnds-analysis/mlcTree_' + str(
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
                os.system("codeml %s/dnds-analysis/%s > /dev/null 2>&1" % (args.outdir, file))
                print("codeml %s/dnds-analysis/%s" % (args.outdir, file))
                print("\n\n")
                # os.system("rm codeml.out")

    # PARSING CODEML OUTPUT

    cwd = os.getcwd()
    DIR = args.outdir + "/dnds-analysis"

    aaiDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    blast = open("%s/pseudogene.blast" % args.outdir)
    for i in blast:
        ls = i.rstrip().split("\t")
        aaiDict[ls[0]] = ls[2]

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
        if float(dndsDict[i]["dn"]) <= args.max_ds and float(dndsDict[i]["ds"]) <= args.max_ds and float(
                dndsDict[i]["ds"]) >= args.min_ds and float(dndsDict[i]["dn"]) >= args.min_ds:
            orf = dndsDict[i]["orf"]
            dn = dndsDict[i]["dn"]
            ds = dndsDict[i]["ds"]
            dnds = float(dndsDict[i]["dn"]) / float(dndsDict[i]["ds"])
            dndsDict2[orf].append(i)
            dndsList.append(dnds)

    dsList = []
    dndsList = []
    total = 0
    out = open(args.outdir + "/dnds-summary.csv", "w")
    out.write("ORF" + "," + "referenceMatch" + "," + "amino_acid_identity" + "," + "dN" + "," + "dS" + "," + "dN/dS" + "," + "PG" + "\n")
    for i in sorted(dndsDict.keys()):
        total += 1

        if float(dndsDict[i]["ds"]) < float(args.max_ds) and float(dndsDict[i]["ds"]) > float(args.min_ds):

            dnds = float(dndsDict[i]["dn"]) / float(dndsDict[i]["ds"])

            dsList.append(float(dndsDict[i]["ds"]))
            dndsList.append(dnds)

            if dnds < args.max_dnds:
                pg = "N"
            else:
                pg = "Y"

            out.write(i + "," + str(dndsDict[i]["orf"]) + "," + str(aaiDict[i]) + "," + str(dndsDict[i]["dn"]) + "," + str(dndsDict[i]["ds"]) + "," +

                      str(dnds) + "," + str(pg) + "\n")

    out.close()

    if len(dsList) > 0:

        print("")
        print("Done!")
        print("********************************************************************************")
        print("Identified a total of %d orthologs between the query and reference datasets," % total)
        print("with %d orthologs between the dS range of %d-%d" % (len(dsList), args.min_ds, args.max_ds))
        print("")
        print("Average dS among orthologs within the specified dS range: " + str(ave(dsList)))
        print("Average dN/dS among orthologs within the specified dS range: " + str(ave(dndsList)))
        print("")
        print("See %s/dnds-summay.csv for detailed results." % args.outdir)
        print("********************************************************************************")
        print("Thanks for using!")

        # os.system("rm -f %s/pseudogene.blast" % args.out)
        os.system("rm -f 2NG.t 2NG.dN 2NG.dS rst1 rst 2ML.t 2ML.dN 2ML.dS 4fold.nuc rub")

    else:
        print("")
        print("Done!")
        print("********************************************************************************")
        print("Identified a total of %d orthologs between the query and reference datasets," % total)
        print("However, there were not orthologs between the dS range of %d-%d" % (args.m, args.M))
        print("Please consider increasing the allowable dS range, and re-running the program")
        print("Thanks for using!")


def full(args, file_dict, log_file_dict=None, skip=False):
    ref = args.reference
    nucOrfs = file_dict['cds_filename']
    pepORFs = file_dict['proteome_filename']
    referenceNucOrfs = file_dict['ref_cds_filename']
    referencePepOrfs = file_dict['ref_proteome_filename']
    c = file_dict['ctl']
    dndsLimit = args.max_dnds
    M = args.max_ds
    m = args.min_ds

    out = file_dict['dnds_out']

    try:
        if args.diamond:
            search = "diamond"
        else:
            search = "blast"
    except AttributeError:
        pass

    if skip == True:
        pass
    else:
        threads = args.threads
        os.system("echo ${ctl} > ctl.txt")
        file = open("ctl.txt")
        for i in file:
            ctl = (i.rstrip())
            ctl = ctl[0:len(ctl)]
        os.system("rm ctl.txt")

        if not re.findall(r'codeml', ctl):
            if not c:
                print("\nCannot locate the necessary \ncodeml control file. Please double check that "
                      "installation\n with the \'setup.sh\' script worked without any errors. If issues \npersist, please "
                      "report an issue on GitHub\n")
                raise SystemExit
            else:
                ctl = c

        faa = pepORFs
        fna = nucOrfs

        if referencePepOrfs != "NA" and referenceNucOrfs != "NA":
            refFaa = referencePepOrfs
            refFna = referenceNucOrfs
        else:
            if ref:

                os.system("prodigal -i %s -a %s-proteins.faa -d %s-proteins.fna" % (ref, ref, ref))

                refFna = ref + "-proteins.fna"
                refFaa = ref + "-proteins.faa"
            else:
                print(
                    "Did not find reference datasets. Please provide these using the \'-ra\' and \'rn\', or \'r\' flag")

        #print("Starting dN/dS analysis pipeline...")

        os.system("mkdir -p " + out)
        os.system("mkdir -p " + out + "/dnds-analysis")

        if search == "blast":
            print("Running BLAST")

            os.system("makeblastdb -dbtype prot -in %s -out %s > /dev/null 2>&1" % (refFaa, refFaa))
            # os.system("rm makeblastdb.out")
            os.system("blastp -query %s -db %s "
                      "-outfmt 6 -out %s/pseudogene.blast -evalue 1E-6 -num_threads %s -max_target_seqs 1  > /dev/null 2>&1" % (
                          faa, refFaa, out, threads))

            os.system("rm %s.psq" % refFaa)
            os.system("rm %s.phr" % refFaa)
            os.system("rm %s.pin" % refFaa)

        elif search == "diamond":
            print("Running DIAMOND")
            os.system(
                "diamond makedb --in %s -d %s" % (refFaa, refFaa))

            os.system("diamond blastp --db %s.dmnd --query %s --outfmt 6 --out %s/pseudogene.blast "
                      "--max-target-seqs 1 --evalue 1E-6 --threads %s > /dev/null 2>&1" % (refFaa, faa, out, threads))

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
        aaiDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
        blast = open("%s/pseudogene.blast" % out)
        for i in blast:
            ls = i.rstrip().split("\t")
            if ls[0] not in prescreened:
                alnLengthDict[ls[0]] = ls[3]
                aaiDict[ls[0]] = ls[2]
                outNUC = open(out + "/dnds-analysis/%s.faa.fna" % ls[0], "w")
                outNUC.write(">" + ls[1] + "\n")
                outNUC.write(refFnaDict[ls[1]] + "\n")
                outNUC.write(">" + ls[0] + "\n")
                outNUC.write(fnaDict[ls[0]] + "\n")
                outNUC.close()

                outAA = open(out + "/dnds-analysis/%s.faa" % ls[0], "w")
                outAA.write(">" + ls[1] + "\n")
                outAA.write(refFaaDict[ls[1]] + "\n")
                outAA.write(">" + ls[0] + "\n")
                outAA.write(faaDict[ls[0]] + "\n")
                outAA.close()

        # ALIGNING PROTEIN SEQUENCES AND CREATING A CODON ALIGNMENT
        print("Aligning files...")
        DIR = out + "/dnds-analysis"
        os.system("for i in %s/*faa; do"
                  " muscle -in $i -out $i.aligned.fa > /dev/null 2>&1;"
                  # " rm muscle.out;"
                  " pal2nal.pl $i.aligned.fa $i.fna -output fasta > $i.codonalign.fa;"
                  " done" % DIR)

        # BUILDING CONTROL FILES
        print("Preparing for codeml analysis")
        DIR = out + "/dnds-analysis"
        codealign = os.listdir(DIR)
        count = 0
        for file in codealign:
            if re.findall(r'codonalign', file):
                count += 1
                clu = file.split(".faa")[0]
                setup = open(ctl)
                OUT = open("%s/%s.ctl" % (DIR, str(clu)), "w")

                for i in setup:
                    if re.findall('seqfile', i):
                        OUT.write(
                            '' + ' ' + '' + ' ' + '' + ' ' + '' + ' ' + '' + ' ' + '' + ' ' + 'seqfile' + ' ' + '= ' + out + '/dnds-analysis/' + file + ' ' + '*' + ' ' + 'sequence' + ' ' + 'data' + ' ' + 'filename\n')

                    elif re.findall(r'outfile', i):
                        OUT.write(
                            '' + ' ' + '' + ' ' + '' + ' ' + '' + ' ' + '' + ' ' + '' + ' ' + 'outfile' + ' ' + '=' + ' '
                            + out + '/dnds-analysis/mlcTree_' + str(
                                clu) + ' ' + '' + ' ' + '' + ' ' + '' + ' ' + '' + ' ' + '' + ' ' + '' + ' ' + '' + ' '
                            + '' + ' ' + '' + ' ' + '' + ' ' + '*' + ' ' + 'main' + ' ' + 'result' + ' ' + 'file' + ' ' + 'name\n')

                    else:
                        OUT.write(i)
                OUT.close()
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

                os.system("codeml %s/dnds-analysis/%s > /dev/null 2>&1" % (out, file))
                # os.system("rm codeml.out")

    # PARSING CODEML OUTPUT
    if log_file_dict:
        DIR = log_file_dict['dnds_out'] + "/dnds-analysis"
    else:
        DIR = out + "/dnds-analysis"

    print("Summarizing codeml output")
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
        if float(dndsDict[i]["dn"]) <= M and float(dndsDict[i]["ds"]) <= M and float(
                dndsDict[i]["ds"]) >= m and float(dndsDict[i]["dn"]) >= m:
            orf = dndsDict[i]["orf"]
            dn = dndsDict[i]["dn"]
            ds = dndsDict[i]["ds"]
            dnds = float(dndsDict[i]["dn"]) / float(dndsDict[i]["ds"])
            dndsDict2[orf].append(i)
            dndsList.append(dnds)

    dsList = []
    dndsList = []
    total = 0
    try:
        OUT = open(out + "/dnds-summary.csv", "w")
    except FileNotFoundError:
        os.makedirs(out)
        OUT = open(out + "/dnds-summary.csv", "w")

    OUT.write("ORF" + "," + "referenceMatch" + "," + "amino_acid_identity" + "," + "dN" + "," + "dS" + "," + "dN/dS" + "," + "PG" + "\n")
    for i in sorted(dndsDict.keys()):
        total += 1

        if float(dndsDict[i]["ds"]) < float(M) and float(dndsDict[i]["ds"]) > float(m):

            dnds = float(dndsDict[i]["dn"]) / float(dndsDict[i]["ds"])

            dsList.append(float(dndsDict[i]["ds"]))
            dndsList.append(dnds)

            if float(dnds) < dndsLimit:
                pg = "N"
            else:
                pg = "Y"

            OUT.write(
                i + "," + str(dndsDict[i]["orf"]) + "," + str(aaiDict[i]) + "," + str(dndsDict[i]["dn"]) + "," + str(dndsDict[i]["ds"]) + "," +

                str(dnds) + "," + str(pg) + "\n")

    OUT.close()

    if len(dsList) > 0:

        print("")
        print("Done!")
        print("********************************************************************************")
        print("Identified a total of %d orthologs between the query and reference datasets," % total)
        print("with %d orthologs between the dS range of %d-%d" % (len(dsList), m, M))
        print("")
        print("Average dS among orthologs within the specified dS range: " + str(ave(dsList)))
        print("Average dN/dS among orthologs within the specified dS range: " + str(ave(dndsList)))
        print("")
        print("See %s/dnds-summary.csv for detailed results." % out)
        print("********************************************************************************")
        print("Thanks for using!")

        # os.system("rm -f %s/pseudogene.blast" % out)
        os.system("rm -f 2NG.t 2NG.dN 2NG.dS rst1 rst 2ML.t 2ML.dN 2ML.dS 4fold.nuc rub")

    else:
        print("")
        print("Done!")
        print("********************************************************************************")
        print("Identified a total of %d orthologs between the query and reference datasets," % total)
        print("However, there were not orthologs between the dS range of %d-%d" % (m, M))
        print("Please consider increasing the allowable dS range, and re-running the program")
        print("Thanks for using!")


