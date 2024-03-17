# Copyrights 2023: Nuttachat Wisittipanit

from const import (CRName, )
from util import (createFolderAndClear, read_file_dlimit, )
from structobj import (SpacerOne, RES_BLAST, CLocus, SpacerProfile, RES_CRISPR, LocusANT, NodeInfo, LocusOne,
                       LocusUSE, )
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from collections import OrderedDict
import re
import os
import util
from copy import deepcopy

"""
amrfinder --nucleotide fasta_file --organism Salmonella -o output_file 

"""


class CrisprCas_CLS:

    def __init__(self, fileCls, inputFilePath, outFolderPath, outFolderName):
        self.fileCls = fileCls
        self.inputFilePath = inputFilePath
        self.outFolderPath = outFolderPath
        self.outFolderName = outFolderName  # BK_SAL1
        # self.spResANT ->
        # {'NODE_1': {'C1':LocusANT ins,'C2':LocusANT ins,},..}
        # }
        self.spResANT = {}
        self.spRes = {}
        # self.spRes ->
        # {'NODE_1': {'C1':['TTA','TCG','AAC',..],'C2':['TTA','TCG','AAC',..]},
        #          'NODE_2': {'C1':['TTA','TCG','AAC',..],'C2':['TTA','TCG','AAC',..]},
        # }
        self.spResObj = {}
        self.spacerOne_list = []  # list of SpacerOne ins (original C1 C2 and C3)
        self.spacerOne_Use = []  # list of adjusted SpacerOne ins (only C1 C2)
        self.spacerRecord = []  # list of SecRecord ins
        self.spProfile = SpacerProfile()  # SpacerProfile ins
        self.spProfile.genomeName = self.outFolderName
        # already have node changed (and sstrand added)
        # self.spResObj_NEW -> {'C1': list of SpacerOne,'C2': ...,'C3': ...}
        # self.spResANT_NEW -> {'C1': LocusANT,'C2': LocusANT,'C3': LocusANT}
        self.spResANT_NEW = {}
        # self.spResANT_USE -> {'C1': LocusANT,'C2': LocusANT}  final form having only two locus
        self.spResANT_USE = {}
        self.spResObj_NEW = {}
        self.spacerOne_NEW = []  # list of SpacerOne ins (already reversed)
        self.nodeOrder_temp = []  # order of nodes on LT2 [('NODE_1', 45685),('NODE_2', 45989)]
        self.nodeOrder = []  # order of nodes on LT2
        # self.locusOrder_gt -> ['C2','C1','C3']
        self.locusOrder_gt = []  # order of all original locus on LT2
        # self.locusOrder_gtall -> [('C2','minus'),('C1','minus'),('C3','plus')]
        self.locusOrder_gtall = []  # order of all original locus on LT2
        # list of LocusOne ins
        self.locusOrder_One = []
        # self.C1_oneList -> list of LocusOne ins for C1
        self.C1_oneList = []
        # self.C2_oneList -> list of LocusOne ins for C2
        self.C2_oneList = []
        # if Node1->(C2,C1) and Node2->(C3)
        # self.cutoffByNodeOrder -> [1]
        # if Node1->(C3,C4) and Node2->(C2) and Node3->(C1) -> (C3,C4,C2,C1)
        # self.cutoffByNodeOrder -> [1,2]
        self.cutoffByNodeOrder = []
        self.gapTH = 1000

    def perform(self):
        """
        Overall perform
        :return:
        """
        # print('Perform CRISPR Finder')
        self.performCrisprReport()
        crisprFound = self.readCrisprReport()
        if crisprFound is False:
            return
        self.transformCrisprRes()
        self.createCLocusFasta()
        self.performBlastSpacerGenome()
        # self.performBlastSpacerGenomeRef()
        self.prepareSpacerFasta()
        self.performBlastSpacer()
        self.performReadBlastSpacerGenome()
        self.prepareNodeFASTA()
        self.performBlastNodeGenomeRef()
        readNode = self.performReadBlastNodeGenomeRef()
        if readNode is False:
            return
        # self.performReadBlastCLocusAndSP()
        newLocus = self.performNewLocusAssignment()
        if newLocus is False:
            return
        self.prepareSpacerFasta_USE()
        self.performBlastSpacer_USE()
        self.assignNewSpacerIDandCreateProfile()
        self.appendSpacerFastaTemp()
        self.makeTemporarySpacerDB()

    def performCrisprReport(self):
        """
        Perform crt tool for each genome assembly
        :return:
        """
        # print("performCrisprReport")
        # print(self.inputFilePath)
        # print(self.outFolderPath)
        # print(self.outFolderName)
        fastaFile = self.inputFilePath
        outFileName = 'crispr_result.txt'
        outFileNamePath = os.path.join(self.outFolderPath, outFileName)
        chkExist = os.path.isfile(outFileNamePath)
        if chkExist is False:
            a = os.path.dirname(__file__)
            b = os.path.realpath(os.path.join(a, '..', 'ext', 'CRT1.2-CLI.jar'))
            cmdLine = 'java -cp %s crt %s %s' % (
                b,
                fastaFile,  # FASTA file
                outFileNamePath,  # output file
            )
            os.system(cmdLine)

    def readCrisprReport(self):
        """
        Read Crispr report
        :return:
        """
        outFileName = 'crispr_result.txt'
        # / home / nuttachat / Python / Project / WGS / out_file / BK_SAL1 / crispr_result.txt
        # onePath = os.path.join(File.outputFolder,'BK_SAL1')
        # print(onePath)
        # outFileNamePath = os.path.join(onePath, outFileName)
        outFileNamePath = os.path.join(self.outFolderPath, outFileName)
        # read text file
        self.spRes = {}
        self.spResANT = {}
        lines = util.read_file_normal(outFileNamePath)
        # Determine if CRISPR is found or not
        crisprFound = True
        for line in lines:
            if line.find('No CRISPR') != -1:
                # found 'No CRISPR' -> there is no goddamn CRISPR
                crisprFound = False
                break
        if crisprFound is False:
            return False
        n = 0
        foundBatch = False  # first, set foundBatch as False
        breakAll = False
        k = True
        while k:
            try:
                line = lines[n]
            except (IndexError,):
                break
            if breakAll is True:
                break
            if line.find('ORGANISM:') != -1:
                # found a batch of CRISPR locus
                # print('foundBatch is True')
                foundBatch = True
                batchName = line.split()[1]  # 'NODE_1_length_1012076_cov_79.7719'
                self.spRes[batchName] = {}
                self.spResANT[batchName] = {}
                n += 2
                continue
            while foundBatch:
                try:
                    line = lines[n]
                except (IndexError,):
                    break
                # print('line: ', line)
                if line.find('ORGANISM') != -1:
                    # found end of reports, break all
                    breakAll = False
                    foundBatch = False
                    break
                elif line.find('Time to find repeats:') != -1:
                    # found end of reports, break all
                    breakAll = True
                    foundBatch = False
                    break
                lineSplit = line.split()
                # print('lineSplit: ', lineSplit)
                if 'CRISPR' in lineSplit:
                    # cIndex -> 'C1','C2' or 'C3'
                    """
                    (NO Reverse for now)
                    Reverse C3&C2 -> C1 and C1 -> C2
                    If there are only C1 and C2 found, then C2 -> C1 and C1 -> C2
                    """
                    cIndex = str(lineSplit[1])
                    posStart = int(lineSplit[3])
                    posEnd = int(lineSplit[5])
                    # print('cIndex: ', cIndex)
                    # cIndexUse = 'X'
                    # if cIndex in ('1',):
                    #     cIndexUse = '2'
                    # elif cIndex in ('2',):
                    #     cIndexUse = '1'
                    # elif cIndex in ('3',):
                    #     cIndexUse = '1X'
                    # print('cIndexUse: ', cIndexUse)
                    # cLocus = 'C%s' % (cIndexUse,)  # C1 (CRISPR locus 1)
                    cLocus = 'C%s' % (cIndex,)
                    # print('cLocus: ', cLocus)
                    self.spRes[batchName][cLocus] = []
                    self.spResANT[batchName][cLocus] = LocusANT()
                    self.spResANT[batchName][cLocus].posStart = posStart
                    self.spResANT[batchName][cLocus].posEnd = posEnd
                    self.spResANT[batchName][cLocus].Node = batchName
                    n += 3
                    continue
                try:
                    spacer = lineSplit[2]
                    if spacer.find('-') == -1 and spacer.find('Average') == -1:
                        # not found
                        self.spRes[batchName][cLocus].append(spacer)
                except (IndexError,):
                    pass
                try:
                    DR_seq = lineSplit[1]
                    if DR_seq.find('-') == -1 and lineSplit[0].find('Repeats:') == -1:
                        # not found
                        self.spResANT[batchName][cLocus].DR.append(DR_seq)
                except (IndexError,):
                    pass
                n += 1
        return True

    def transformCrisprRes(self):
        """
        Transform self.spRes into
        self.spResObj ->
        {'NODE_1': [CLocus obj,CLocus obj,..],..}
        from
        {'NODE_1': {'C1':['TTA','TCG','AAC',..],'C2':['TTA','TCG','AAC',..]},
         'NODE_2': {'C1':['TTA','TCG','AAC',..],'C2':['TTA','TCG','AAC',..]},
        }
        :return:
        """
        for node in self.spRes:
            # key -> 'NODE_1'
            self.spResObj[node] = []
            for locus in self.spRes[node]:
                # locus -> 'C1'
                locusObj = CLocus()
                locusObj.name = locus
                spacerList = self.spRes[node][locus]
                for k in range(len(spacerList)):
                    sp = spacerList[k]
                    spObj = SpacerOne()
                    # spID -> 'C1_0_BK_SAL1'
                    spObj.spID = '%s_%s_%s' % (locus, k, self.outFolderName,)
                    spObj.desc = self.outFolderName  # genome name
                    spObj.spSequence = Seq(sp)
                    spObj.spLength = len(sp)
                    spObj.locusName = locus
                    locusObj.spList.append(spObj)
                self.spResObj[node].append(locusObj)
            break
        # fill self.spacerOne_NEW -> list of SpacerOne ins
        for node in self.spResObj:
            # node -> 'NODE_1'
            # self.spResObj[node] -> [CLocus obj,CLocus obj,..]
            for cs in self.spResObj[node]:
                # cs -> CLocus ins
                self.spResObj_NEW[cs.name] = deepcopy(cs.spList)
                # for so in cs.spList:
                #     # filter out only C1 and C2 locus
                #     if so.locusName not in CRName.CList:
                #         continue
                #     # keep all SpacerOne ins in the list (self.spacerOne_list)
                #     self.spacerOne_NEW.append(so)

    def combineDRandSP(self, DR, spList):
        """
        Combine DR and SP
        DR -> list of direct repeat sequences
        spList -> list of spacers
        Need to build sequence from start to finish
        sequence = DR[0]+spList[0]+DR[1]+spList[1]+DR[2]+spList[2]+DR[3]
        """
        sequence = ''
        for i, seqDR in enumerate(DR):
            # i->0,seqDR->'ACTG' of DR
            try:
                sequence = sequence + seqDR + spList[i]
            except (IndexError,):
                sequence = sequence + seqDR
        return sequence

    def createCLocusFasta(self):
        """
        Take the first spacer of each locus and create a FASTA file
        .spRes -> {'NODE_41_': {'C1': ['ACC..','ACG..'],'C2': ['ACC..','ACG..'],'C3': ['ACC..','ACG..']}
                  }
        .spResANT -> {'NODE_41_': {'C1': LocusANT,'C2': LocusANT,'C3': LocusANT}
                      }
        self.spResObj -> {'NODE_1': [CLocus obj,CLocus obj,..],..}
        LocusANT contains posStart, posEnd, DR, Node
        """
        print('createCLocusFasta')
        # eachSpaperRecord -> list of SeqRecord
        eachSpaperRecord = []
        locusFound = []  # locus found so far
        for node in self.spRes:
            print('node: ', node)
            # node -> 'NODE_41'
            # self.spRes[node] -> {'C1': ['ACC..','ACG..'],'C2': ['ACC..','ACG..'],'C3': ['ACC..','ACG..']}
            # print('spRes[node]:')
            # print(self.spRes[node])
            for locus in self.spRes[node]:
                # locus -> 'C1','C2' or 'C3'
                # self.spRes[node][locus] -> ['ACC..','ACG..'] list of spacer
                # self.spResANT[node][locus] -> LocusANT ins
                # LocusANT.DR -> list of DR sequences
                spList = self.spRes[node][locus]  # list of spacer
                L = self.spResANT[node][locus]  # LocusANT ins
                L.name = locus  # C1
                L.SP = self.spResObj_NEW[locus]
                # seqUse -> DR1+SP1+DR2+SP2
                # if locus not in CRName.CList:
                #     continue
                # print('locus: ', locus)
                # print(spList)
                # print(L)
                # print(L.DR)
                combineSeq = self.combineDRandSP(L.DR, spList)
                # seqUse = Seq(L.DR[0] + spList[0] + L.DR[1] + spList[1])
                seqUse = Seq(combineSeq)
                # print('seqUse: ', seqUse)
                if locus == CRName.C1 and CRName.C1 not in locusFound:
                    # print('axa C1')
                    record = SeqRecord(seqUse, id=CRName.C1, description='')
                    locusFound.append(CRName.C1)
                elif locus == CRName.C2 and CRName.C2 not in locusFound:
                    record = SeqRecord(seqUse, id=CRName.C2, description='')
                    locusFound.append(CRName.C2)
                elif locus == CRName.C3 and CRName.C3 not in locusFound:
                    record = SeqRecord(seqUse, id=CRName.C3, description='')
                    locusFound.append(CRName.C3)
                elif locus == CRName.C4 and CRName.C4 not in locusFound:
                    record = SeqRecord(seqUse, id=CRName.C4, description='')
                    locusFound.append(CRName.C4)
                eachSpaperRecord.append(record)
        # Write into the file
        fastaFilePath = os.path.join(self.outFolderPath, 'cLocusFasta.fasta')
        with open(fastaFilePath, "w") as output_handle:
            SeqIO.write(eachSpaperRecord, output_handle, "fasta")

    def performBlastSpacerGenome(self):
        """
        Perform Blast to assign spacer ID.
        No need to grow a database
        :return:
        """
        fastaFile = 'cLocusFasta.fasta'
        fastaFilePath = os.path.join(self.outFolderPath, fastaFile)
        # dbFolderPath -> 'maindb/DBSC/BK_SAL1_DB
        dbFolderPath = os.path.join(self.fileCls.DBSC_Folder, '%s_DB' % self.outFolderName)
        # dbFolderUsePath -> 'maindb/DBSC/BK_SAL1_DB/BK_SAL1
        dbFolderUsePath = os.path.join(dbFolderPath, self.outFolderName)
        # outputFolderPath = os.path.join(File.outBlastFolder, genomeName)
        # createFolderAndClear(outputFolderPath)
        resultFilePath = os.path.join(self.outFolderPath, 'cLocusBlastResult.txt')
        # SPI1	NODE_2_length_635729_cov_82.5337	0.0	60327	98.496	82	82	4223	38434	540680	506466 plus
        # cmdLine = 'blastn -query %s -db %s -num_threads 6 -word_size 7 -strand both -perc_identity 99 -qcov_hsp_perc 99 -out %s -outfmt "6 qacc sallseqid evalue bitscore pident qcovs qcovhsp qstart qend sstart send sstrand"' % (
        #     fastaFilePath,  # FASTA file
        #     dbFolderUsePath,  # Database file
        #     resultFilePath,)  # Output file
        cmdLine = 'blastn -query %s -db %s -num_threads 6 -word_size 7 -strand both -perc_identity 99 -qcov_hsp_perc 99 -out %s -outfmt "6 qacc sseqid evalue bitscore pident qcovs qcovhsp qstart qend sstart send sstrand"' % (
            fastaFilePath,  # FASTA file
            dbFolderUsePath,  # Database file
            resultFilePath,)  # Output file
        os.system(cmdLine)

    # def performBlastSpacerGenomeRef(self):
    #     """
    #     Perform Blast to assign spacer ID.
    #     No need to grow a database
    #     :return:
    #     """
    #     fastaFile = 'cLocusFasta.fasta'
    #     fastaFilePath = os.path.join(self.outFolderPath, fastaFile)
    #     # dbFolderPath -> 'maindb/REF/LT2DB
    #     dbFolderPath = os.path.join(self.fileCls.REFDBFolder, 'LT2DB')
    #     # dbFolderUsePath -> 'maindb/REF/LT2DB/LT2
    #     dbFolderUsePath = os.path.join(dbFolderPath, 'LT2')
    #     # outputFolderPath = os.path.join(File.outBlastFolder, genomeName)
    #     # createFolderAndClear(outputFolderPath)
    #     resultFilePath = os.path.join(self.outFolderPath, 'cLocusBlastResult_REF.txt')
    #     # SPI1	NODE_2_length_635729_cov_82.5337	0.0	60327	98.496	82	82	4223	38434	540680	506466 plus
    #     cmdLine = 'blastn -query %s -db %s -num_threads 4 -word_size 7 -strand both -perc_identity 100 -qcov_hsp_perc 80 -out %s -outfmt "6 qacc sallseqid evalue bitscore pident qcovs qcovhsp qstart qend sstart send sstrand"' % (
    #         fastaFilePath,  # FASTA file
    #         dbFolderUsePath,  # Database file
    #         resultFilePath,)  # Output file
    #     os.system(cmdLine)

    def prepareSpacerFasta(self):
        """
        Prepare fasta file of all discovered spacers for the Blast operation. So this means I need to transform all
        CLocus instances in NODE_1 (in self.spResObj) into a FASTA file.
        self.spResObj_NEW -> {'C1': list of SpacerOne,'C2': ...,'C3': ...}
        self.spResANT_NEW -> {'C1': LocusANT,'C2': LocusANT,'C3': LocusANT}
        # self.spUse -> {'C1': list of SpacerOne,'C2': ...,'C3': ...} final form
        :return:
        """
        for locus in self.spResObj_NEW:
            # cs -> 'C1'
            for so in self.spResObj_NEW[locus]:
                # print(so.spID)
                # print(so.spSequence)
                # filter out only C1 and C2 locus
                # if so.locusName not in CRName.CList:
                #     continue
                record = SeqRecord(so.spSequence, id=so.spID, description=so.desc)
                self.spacerRecord.append(record)
                # keep all SpacerOne ins in the list (self.spacerOne_list)
                self.spacerOne_list.append(so)
        # # Reverse the list
        # self.spacerRecord.reverse()
        # self.spacerOne_list.reverse()
        # Write into the file
        fastaFilePath = os.path.join(self.outFolderPath, 'spacerfound.fasta')
        with open(fastaFilePath, "w") as output_handle:
            SeqIO.write(self.spacerRecord, output_handle, "fasta")

    def performBlastSpacer(self):
        """
        Perform Blast to assign spacer ID.
        No need to grow a database
        :return:
        """
        fastaFile = 'spacerfound.fasta'
        fastaFilePath = os.path.join(self.outFolderPath, fastaFile)
        # dbFolderPath -> 'SPACER/SP_DB
        dbFolderPath = os.path.join(self.fileCls.SPACER_Folder, 'SPTEMP_DB')
        # dbFolderUsePath -> 'SPACER/SP_DB/SPALL
        dbFolderUsePath = os.path.join(dbFolderPath, 'SPALL')
        # outputFolderPath -> 'maindb/out_blast/BK_SAL1'
        # outputFolderPath = os.path.join(File.outBlastFolder, genomeName)
        # createFolderAndClear(outputFolderPath)
        resultFilePath = os.path.join(self.outFolderPath, 'spBlastResult.txt')
        # SPI1	NODE_2_length_635729_cov_82.5337	0.0	60327	98.496	82	82	4223	38434	540680	506466
        cmdLine = 'blastn -query %s -db %s -num_threads 4 -word_size 7 -strand both -perc_identity 100 -qcov_hsp_perc 80 -out %s -outfmt "6 qacc sallseqid evalue bitscore pident qcovs qcovhsp qstart qend sstart send sstrand"' % (
            fastaFilePath,  # FASTA file
            dbFolderUsePath,  # Database file
            resultFilePath,)  # Output file
        os.system(cmdLine)

    def performReadBlastSpacerGenome(self):
        """
        Read 'cLocusBlastResult.txt' and 'spBlastResult.txt'
        self.spResANT -> {'NODE_41_': {'C1': LocusANT,'C2': LocusANT,'C3': LocusANT}
                        }
        self.spResANT_NEW -> {'C1': LocusANT,'C2': LocusANT,'C3': LocusANT}
        LocusANT contains posStart, posEnd, DR, Node
        LocusANT now contains all info about one locus
        """
        # print('performReadBlastSpacerGenome')
        # print('self.spResANT')
        # print(self.spResANT)
        # copy self.spResANT to self.spResANT_NEW
        # UPDATE self.spResANT_NEW
        for node in self.spResANT:
            # node -> 'NODE_41_'
            # nsp -> {'C1': LocusANT,'C2': LocusANT,'C3': LocusANT}
            print('node: ', node)
            nsp = self.spResANT[node]
            for locus in nsp:
                print('locus: ', locus)
                # locus -> 'C1'
                # nsp[locus] -> LocusANT
                self.spResANT_NEW[locus] = deepcopy(nsp[locus])
        for locus in self.spResANT_NEW:
            self.spResANT_NEW[locus].DR_One = self.spResANT_NEW[locus].DR[0]
        blastResultFilePath = os.path.join(self.outFolderPath, 'cLocusBlastResult.txt')
        # xBlast -> list of RES_BLAST ins
        xBlast = util.readBlastResult(blastResultFilePath)
        spBlastFilePath = os.path.join(self.outFolderPath, 'spBlastResult.txt')
        # xBlast_sp -> list of RES_BLAST ins
        xBlast_sp = util.readBlastResult(spBlastFilePath)
        # get the strand of each locus
        xLocus = {}  # {'C1': 'plus',..}
        for x in xBlast_sp:
            # x.qacc -> C1_0_SAL_HC6463AA_AS
            # x.sstrand -> 'plus' or 'minus'
            locus = x.qacc.split('_')[0]  # 'C1'
            if xLocus.get(locus) is None:
                xLocus[locus] = x.sstrand
        for locus in xLocus:
            self.spResANT_NEW[locus].locusStrand = xLocus[locus]
        # bLocus -> {'C1': (NODE_12,plus),..} get only the first one
        bLocus = {}
        for cx in xBlast:
            locus = cx.qacc
            if bLocus.get(locus) is None:
                subject = cx.sallseqid  # NODE_12
                strand = cx.sstrand  # plus or minus
                bLocus[locus] = (subject, strand, cx.sstart, cx.send)
        for locus in bLocus:
            self.spResANT_NEW[locus].Node = bLocus[locus][0]
            self.spResANT_NEW[locus].NodeStrand = bLocus[locus][1]
            self.spResANT_NEW[locus].NodeStart = bLocus[locus][2]
            self.spResANT_NEW[locus].NodeEnd = bLocus[locus][3]
        # print('self.spResANT_NEW')
        # print(self.spResANT_NEW)
        # for locus in self.spResANT_NEW:
        #     L = self.spResANT_NEW[locus]
        #     print('Node: ', L.Node)
        #     print('NodeStrand: ', L.NodeStrand)
        #     print('NodeStart: ', L.NodeStart)
        #     print('NodeEnd: ', L.NodeEnd)

    def prepareNodeFASTA(self):
        """
        Create the fasta file from node names found in self.spResANT_NEW getting sequnce data from input
        """
        NodeList = []  # list of Nodes
        for locus in self.spResANT_NEW:
            # L -> LocusANT
            L = self.spResANT_NEW[locus]
            NodeList.append(L.Node)
        # print('NodeList: ', NodeList)
        # find the number of node
        # isSingleNode = False
        # nodeList_source = []  # list of node id from original input file
        # for seq_record in SeqIO.parse(self.inputFilePath, "fasta"):
        #     nodeList_source.append(seq_record.id)
        # print('nodeList_source: ', nodeList_source)
        # if len(nodeList_source) == 1:
        #     isSingleNode = True
        # print('isSingleNode: ', isSingleNode)
        srList = []
        for seq_record in SeqIO.parse(self.inputFilePath, "fasta"):
            # print('seq_record.id: ', seq_record.id)
            if seq_record.id in NodeList:
                srList.append(seq_record)
        # Write into the file
        fastaFilePath = os.path.join(self.outFolderPath, 'NodeFound.fasta')
        with open(fastaFilePath, "w") as output_handle:
            SeqIO.write(srList, output_handle, "fasta")

    def performBlastNodeGenomeRef(self):
        """
        Perform Blast to assign spacer ID.
        No need to grow a database
        :return:
        """
        fastaFile = 'NodeFound.fasta'
        fastaFilePath = os.path.join(self.outFolderPath, fastaFile)
        # dbFolderPath -> 'maindb/REF/LT2DB
        dbFolderPath = os.path.join(self.fileCls.REFDBFolder, 'LT2DB')
        # dbFolderUsePath -> 'maindb/REF/LT2DB/LT2
        dbFolderUsePath = os.path.join(dbFolderPath, 'LT2')
        # outputFolderPath = os.path.join(File.outBlastFolder, genomeName)
        # createFolderAndClear(outputFolderPath)
        resultFilePath = os.path.join(self.outFolderPath, 'NodeBlastResult_REF.txt')
        # SPI1	NODE_2_length_635729_cov_82.5337	0.0	60327	98.496	82	82	4223	38434	540680	506466 plus
        cmdLine = 'blastn -query %s -db %s -num_threads 6 -strand both -perc_identity 10 -qcov_hsp_perc 5 -out %s -outfmt "6 qacc sseqid evalue bitscore pident qcovs qcovhsp qstart qend sstart send sstrand"' % (
            fastaFilePath,  # FASTA file
            dbFolderUsePath,  # Database file
            resultFilePath,)  # Output file
        os.system(cmdLine)

    def performReadBlastNodeGenomeRef(self):
        """
        Read 'cLocusBlastResult.txt' and 'spBlastResult.txt'
        self.spResANT -> {'NODE_41_': {'C1': LocusANT,'C2': LocusANT,'C3': LocusANT}
                      }
        self.spResANT_NEW -> {'C1': LocusANT,'C2': LocusANT,'C3': LocusANT}
        LocusANT contains posStart, posEnd, DR, Node
        LocusANT now contains all info about one locus
        """
        # copy self.spResANT to self.spResANT_NEW
        # UPDATE self.spResANT_NEW
        # self.NodeInfo -> {'NODE_1': NodeInfo ins,..}
        self.nodeInfo = {}
        blastResultFilePath = os.path.join(self.outFolderPath, 'NodeBlastResult_REF.txt')
        # xBlast -> list of RES_BLAST ins
        xBlast = util.readBlastResult(blastResultFilePath)
        # select the first one
        # xNode -> {'NODE_12': (start,stop,strand),..}
        xNode = {}
        for x in xBlast:
            # x.qacc -> C1_0_SAL_HC6463AA_AS
            # x.sstrand -> 'plus' or 'minus'
            nodeName = x.qacc
            if xNode.get(nodeName) is None:
                xNode[nodeName] = (x.qstart, x.qend, x.sstart, x.send, x.sstrand)
        for node in xNode:
            self.nodeInfo[node] = NodeInfo()
            self.nodeInfo[node].name = node
            self.nodeInfo[node].nodeStart = xNode[node][0]
            self.nodeInfo[node].nodeEnd = xNode[node][1]
            self.nodeInfo[node].posStart = xNode[node][2]
            self.nodeInfo[node].posEnd = xNode[node][3]
            self.nodeInfo[node].sstrand = xNode[node][4]
        # assign nodestrand
        for locus in self.spResANT_NEW:
            # locus -> 'C1','C2'
            # search for node in self.nodeInfo
            for node in self.nodeInfo:
                if node == self.spResANT_NEW[locus].Node:
                    self.spResANT_NEW[locus].NodeStrand_use = self.nodeInfo[node].sstrand
                    break
        # find original locus order on the node ['C2','C1'] or self.locusOrder_g = []
        for node in self.nodeInfo:
            # print('node: ', node)
            # print('name: ', self.nodeInfo[node].name)
            for locus in self.spResANT_NEW:
                if node == self.spResANT_NEW[locus].Node:
                    nStart = self.spResANT_NEW[locus].NodeStart
                    self.nodeInfo[node].locusOrder_gtemp.append((locus, nStart))
        # sort self.nodeInfo[node].locusOrder_gtemp
        try:
            for node in self.nodeInfo:
                print('node: ', node)
                # print('name: ', self.nodeInfo[node].name)
                print('name: ', self.nodeInfo[node].locusOrder_gtemp)
                if self.nodeInfo[node].sstrand == 'minus':
                    self.nodeInfo[node].locusOrder_gtemp.sort(key=lambda x: int(x[1]), reverse=True)
                else:
                    self.nodeInfo[node].locusOrder_gtemp.sort(key=lambda x: int(x[1]), reverse=False)
                self.nodeInfo[node].locusOrder_g = [u[0] for u in self.nodeInfo[node].locusOrder_gtemp]
        except:
            return False
        # find node order
        self.nodeOrder_temp = []
        for node in self.nodeInfo:
            pStart = self.nodeInfo[node].posStart
            self.nodeOrder_temp.append((node, pStart))
        self.nodeOrder_temp.sort(key=lambda x: int(x[1]), reverse=False)
        self.nodeOrder = [u[0] for u in self.nodeOrder_temp]
        # find self.locusOrder_gt = []  order of all original locus on LT2
        for node in self.nodeOrder:
            # node -> 'SAL_HC6463AA_AS_NODE_12_length_59827_cov_11.211156'
            # self.nodeInfo[node].locusOrder_g -> ['C2', 'C1']
            order_g = self.nodeInfo[node].locusOrder_g
            for g in order_g:
                self.locusOrder_gt.append(g)
                self.locusOrder_gtall.append((g, self.nodeInfo[node].sstrand))
        # self.cutoffByNodeOrder
        # find list of LocusOne
        cIndex = 0
        for node in self.nodeOrder:
            # node -> NODE_21
            # order_g -> ['C2', 'C1']
            order_g = self.nodeInfo[node].locusOrder_g
            # cIndex = len(order_g) - 1
            sstrand = self.nodeInfo[node].sstrand
            nodeStart = self.nodeInfo[node].nodeStart  # 6 / 83
            nodeEnd = self.nodeInfo[node].nodeEnd  # 50225 / 17197
            posStart = self.nodeInfo[node].posStart  # 3026793 / 3094325
            posEnd = self.nodeInfo[node].posEnd  # 3077012  / 3077210
            for i, locus in enumerate(order_g):
                # cOne -> LocusOne ins
                cOne = LocusOne(locus)
                cOne.sstrand = sstrand
                nStart = self.spResANT_NEW[locus].NodeStart  # 50422 / 16253
                nEnd = self.spResANT_NEW[locus].NodeEnd  # 50878 / 17197
                cOne.seqStart = nStart  # 50422
                cOne.seqEnd = nEnd  # 50878
                # find refStart and refEnd
                if sstrand == 'plus':
                    cOne.refStart = posStart + (nStart - nodeStart)
                    cOne.refEnd = posStart + (nEnd - nodeStart)
                elif sstrand == 'minus':
                    cOne.refStart = posStart - (nEnd - nodeStart)
                    cOne.refEnd = posStart - (nStart - nodeStart)
                self.locusOrder_One.append(cOne)
                if i == len(order_g) - 1:
                    # last one
                    self.cutoffByNodeOrder.append(cIndex)
                cIndex += 1
        # remove the last element from self.cutoffByNodeOrder
        try:
            self.cutoffByNodeOrder.pop()
        except (IndexError,):
            # empty list
            pass
        return True

    def performNewLocusAssignment(self):
        """
        self.spResANT_NEW -> {'C1': LocusANT,'C2': LocusANT,'C3': LocusANT}
        LocusANT contains posStart, posEnd, DR, Node
        LocusANT now contains all info about one locus
        # self.NodeInfo -> {'C1': NodeInfo ins,..}

        Need to get this
        self.spResANT_USE -> {'C1': LocusANT,'C2': LocusANT}  final form
        """
        # print('self.spResANT_NEW')
        # print(self.spResANT_NEW)
        # for locus in self.spResANT_NEW:
        #     print('locus: ', locus)
        #     print(self.spResANT_NEW[locus].name)
        #     print(self.spResANT_NEW[locus].posStart)
        #     print(self.spResANT_NEW[locus].posEnd)
        #     print(self.spResANT_NEW[locus].DR)
        #     print(self.spResANT_NEW[locus].DR_One)
        #     print(self.spResANT_NEW[locus].SP)
        #     print(self.spResANT_NEW[locus].Node)
        #     print(self.spResANT_NEW[locus].NodeStrand)
        #     print(self.spResANT_NEW[locus].NodeStart)
        #     print(self.spResANT_NEW[locus].NodeEnd)
        #     print(self.spResANT_NEW[locus].locusStrand)
        #     print(self.spResANT_NEW[locus].NodeStrand_use)
        # print('self.cutoffByNodeOrder: ', self.cutoffByNodeOrder)
        # print('self.nodeInfo')
        # # print(self.nodeInfo)
        # for node in self.nodeInfo:
        #     print('node: ', node)
        #     print('name: ', self.nodeInfo[node].name)
        #     print('posStart: ', self.nodeInfo[node].posStart)
        #     print('posEnd: ', self.nodeInfo[node].posEnd)
        #     print('sstrand: ', self.nodeInfo[node].sstrand)
        #     print('locusOrder_gtemp: ', self.nodeInfo[node].locusOrder_gtemp)
        #     print('locusOrder_g: ', self.nodeInfo[node].locusOrder_g)
        # print('self.nodeOrder on LT2')
        # print(self.nodeOrder)
        # print('self.locusOrder_gt on LT2')
        # print(self.locusOrder_gt)
        # print('self.locusOrder_gtall on LT2')
        # print(self.locusOrder_gtall)
        # print('self.locusOrder_One')
        # for one in self.locusOrder_One:
        #     print('locus: ', one.locus)
        #     print('sstrand: ', one.sstrand)
        #     print('seqStart: ', one.seqStart)
        #     print('seqEnd: ', one.seqEnd)
        #     print('refStart: ', one.refStart)
        #     print('refEnd: ', one.refEnd)
        # group into C1 and C2
        # find gap between locus
        gapList = []
        for i, one in enumerate(self.locusOrder_One):
            a = one
            try:
                b = self.locusOrder_One[i + 1]
            except (IndexError,):
                break
            gapValue = b.refStart - a.refEnd
            gapList.append(gapValue)
        print('gapList')
        print(gapList)
        if len(gapList) == 0:
            return False
        # group locus into C1 and C2
        """
        Start assigning cutOffIndex
        if cutOffIndex is None:
            Cannot find gap cut off -> nodeinfo is messed up -> use Locus on Node instead
        2 cases: 
            (1) Locus are on different node e.g. Node1->(C2,C1),Node2->(C3) cutOffIndex = 1
            (2) Locus are on the same node e.g. Node1->(C2,C1
        self.locusOrder_gt on LT2
        ['C3', 'C4', 'C2', 'C1']
        gapList
        [14, -454, 16620]
        """
        cutOffIndex = None
        for i, gapValue in enumerate(gapList):
            if gapValue > self.gapTH:
                cutOffIndex = i  # 2
                break
        if cutOffIndex is None:
            # if finding by gap fails
            try:
                cutOffIndex = self.cutoffByNodeOrder[0]
            except:
                return False
        print('cutOffIndex: ', cutOffIndex)
        self.C1_oneList = self.locusOrder_One[0:cutOffIndex + 1]
        self.C2_oneList = self.locusOrder_One[cutOffIndex + 1:]
        # print('self.C1_oneList')
        # for one in self.C1_oneList:
        #     print('locus: ', one.locus)
        #     print('sstrand: ', one.sstrand)
        #     print('seqStart: ', one.seqStart)
        #     print('seqEnd: ', one.seqEnd)
        #     print('refStart: ', one.refStart)
        #     print('refEnd: ', one.refEnd)
        # print('##############################')
        # print('self.C2_oneList')
        # for one in self.C2_oneList:
        #     print('locus: ', one.locus)
        #     print('sstrand: ', one.sstrand)
        #     print('seqStart: ', one.seqStart)
        #     print('seqEnd: ', one.seqEnd)
        #     print('refStart: ', one.refStart)
        #     print('refEnd: ', one.refEnd)

        # self.C1_oneList -> list of LocusOne ins for C1
        # self.C2_oneList -> list of LocusOne ins for C2
        # self.spResANT_NEW -> {'C1': LocusANT,'C2': LocusANT,'C3': LocusANT}
        # self.spResANT_USE -> {'C1': LocusUSE,'C2': LocusUSE}

        for locus in CRName.CListUse:
            # locus -> 'C1' or 'C2'
            self.spResANT_USE[locus] = LocusUSE()
            self.spResANT_USE[locus].name = locus
            # self.spResANT_USE[locus].SP -> list of SpacerOne ins
            if locus == CRName.C1:
                # oneList -> LocusOne ins
                oneList = self.C1_oneList
            else:
                oneList = self.C2_oneList
            # find SP (list of SpacerOne ins) from combining all LocusOne in oneList
            SP_tempall = []
            # DR_tempall -> list of DR seq
            DR_tempall = []
            Node_tempall = []
            for one in oneList:
                # one -> LocusOne ins
                name = one.locus
                # strand -> 'plus' or 'minus'
                strand = one.sstrand
                # LA -> LocusANT ins
                LA = self.spResANT_NEW[name]
                Node_tempall.append(LA.Node)
                if strand == 'plus':
                    # no reverse
                    SP_temp = deepcopy(LA.SP)
                    DR_temp = deepcopy(LA.DR)
                else:
                    # reverse
                    SP_temp = deepcopy(LA.SP)
                    SP_temp.reverse()
                    DR_temp = deepcopy(LA.DR)
                    DR_temp.reverse()
                for k, so in enumerate(SP_temp):
                    # so.locusName_new = locus
                    # so.spID_adj = '%s_%s_%s' % (locus, k, self.outFolderName,)
                    # so.desc_new = self.outFolderName
                    SP_tempall.append(so)
                for seq in DR_temp:
                    DR_tempall.append(seq)
            for k, so in enumerate(SP_tempall):
                so.locusName_new = locus
                so.spID_adj = '%s_%s_%s' % (locus, k, self.outFolderName,)
                so.desc_new = self.outFolderName
                # SP_tempall.append(so)
            self.spResANT_USE[locus].SP = SP_tempall
            self.spResANT_USE[locus].DR = DR_tempall
            self.spResANT_USE[locus].DR_One = DR_tempall[0]
            self.spResANT_USE[locus].NodeList = list(set(Node_tempall))
            # find refStart and refEnd
            self.spResANT_USE[locus].refStart = oneList[0].refStart
            self.spResANT_USE[locus].refEnd = oneList[-1].refEnd

        # for locus in self.spResANT_USE:
        #     print('locus: ', locus)
        #     print('name: ', self.spResANT_USE[locus].name)
        #     print('SP first: ', self.spResANT_USE[locus].SP[0].spID)
        #     print('SP firstNew: ', self.spResANT_USE[locus].SP[0].spID_adj)
        #     print('SP last: ', self.spResANT_USE[locus].SP[-1].spID)
        #     print('SP lastNew: ', self.spResANT_USE[locus].SP[-1].spID_adj)
        #     print(self.spResANT_USE[locus].DR)
        #     print(self.spResANT_USE[locus].DR_One)
        #     print(self.spResANT_USE[locus].NodeList)
        #     print(self.spResANT_USE[locus].refStart)
        #     print(self.spResANT_USE[locus].refEnd)
        # print('##########################################33')
        return True

    def prepareSpacerFasta_USE(self):
        """
        Prepare fasta file of all discovered spacers for the Blast operation. So this means I need to transform all
        CLocus instances in NODE_1 (in self.spResObj) into a FASTA file.
        self.spResObj_NEW -> {'C1': list of SpacerOne,'C2': ...,'C3': ...}
        self.spResANT_NEW -> {'C1': LocusANT,'C2': LocusANT,'C3': LocusANT}
        # self.spUse -> {'C1': list of SpacerOne,'C2': ...,'C3': ...} final form
        :return:
        """
        spacerRecord = []
        for locus in self.spResANT_USE:
            # cs -> 'C1'
            for so in self.spResANT_USE[locus].SP:
                # print(so.spID)
                # print(so.spSequence)
                # filter out only C1 and C2 locus
                # if so.locusName not in CRName.CList:
                #     continue
                record = SeqRecord(so.spSequence, id=so.spID_adj, description=so.desc)
                spacerRecord.append(record)
                # keep all SpacerOne ins in the list (self.spacerOne_list)
                self.spacerOne_Use.append(so)
        # # Reverse the list
        # self.spacerRecord.reverse()
        # self.spacerOne_list.reverse()
        # Write into the file
        fastaFilePath = os.path.join(self.outFolderPath, 'spacerfound_use.fasta')
        with open(fastaFilePath, "w") as output_handle:
            SeqIO.write(spacerRecord, output_handle, "fasta")

    def performBlastSpacer_USE(self):
        """
        Perform Blast to assign spacer ID.
        No need to grow a database
        :return:
        """
        fastaFile = 'spacerfound_use.fasta'
        fastaFilePath = os.path.join(self.outFolderPath, fastaFile)
        # dbFolderPath -> 'SPACER/SP_DB
        dbFolderPath = os.path.join(self.fileCls.SPACER_Folder, 'SPTEMP_DB')
        # dbFolderUsePath -> 'SPACER/SP_DB/SPALL
        dbFolderUsePath = os.path.join(dbFolderPath, 'SPALL')
        # outputFolderPath -> 'maindb/out_blast/BK_SAL1'
        # outputFolderPath = os.path.join(File.outBlastFolder, genomeName)
        # createFolderAndClear(outputFolderPath)
        resultFilePath = os.path.join(self.outFolderPath, 'spBlastResult_use.txt')
        # SPI1	NODE_2_length_635729_cov_82.5337	0.0	60327	98.496	82	82	4223	38434	540680	506466
        cmdLine = 'blastn -query %s -db %s -num_threads 4 -word_size 7 -strand both -perc_identity 100 -qcov_hsp_perc 80 -out %s -outfmt "6 qacc sallseqid evalue bitscore pident qcovs qcovhsp qstart qend sstart send sstrand"' % (
            fastaFilePath,  # FASTA file
            dbFolderUsePath,  # Database file
            resultFilePath,)  # Output file
        os.system(cmdLine)

    def assignNewSpacerIDandCreateProfile(self):
        """
        Purely assign spacer ID using self.spacerOne_List
        (1) Determine the match for each spID (C1_0_BKSAL1). If there is one match, then use the sallseqid attribute.
            If there are > one match, choose one with the lowest e-value (usually the first one on the list).
        (2) Assign the match (sallseqid) to SpacerOne.spID_new attribute.
        (3) For the unmatched (exist in self.spacerOne_List but not in the blast result), just copy SpacerOne.spID
            to SpacerOne.spID_New.
        (4) Construct self.spProfile from self.spacerOne_List
        :return:
        """
        blastResultFilePath = os.path.join(self.outFolderPath, 'spBlastResult_use.txt')
        # obtain all IDs from blastResultFilePath
        # xBlast -> list of RES_BLAST ins
        xBlast = util.readBlastResult(blastResultFilePath)
        # [(C1_0,STMB33),(C1_1,STMB19),..]
        # match -> {'C1_0':[STMB22],'C1_12':[STMB10,STMB10var1],..}
        match = {}
        # self.spacerOne_list -> list of all SpacerOne objs found in a genome
        for resB in xBlast:
            # resB -> RES_BLAST ins
            # resB.qacc -> C1_0_SAL_BA3995AA_AS
            # resB.sallseqid -> STM1
            if match.get(resB.qacc) is None:
                match[resB.qacc] = []
            match[resB.qacc].append(resB.sallseqid)
        # print('match')
        # print(match)
        # matchSelect -> {'C1_0':'STMB22','C1_12':'STMB10',..}
        matchSelect = {}
        for spID in match:
            matchSelect[spID] = match[spID][0]
        # print('matchSelect')
        # print(matchSelect)
        # If so.spID is present in matchSelect dict, the ID exists and will be assigned the matched ID.
        # Else: assign spID_new the spID
        for so in self.spacerOne_Use:
            # so -> SpacerOne ins
            # so.spID -> 'C1_0'
            if matchSelect.get(so.spID_adj) is not None:
                so.spID_new = matchSelect[so.spID_adj]
                so.foundMatch = True
            else:
                so.spID_new = so.spID_adj
                so.foundMatch = False
        # for so in self.spacerOne_list:
        #     print('old ID ', so.spID)
        #     print('new ID ', so.spID_new)
        # Create the spacer profile based on spID_new
        for so in self.spacerOne_Use:
            if so.locusName_new == CRName.C1:
                self.spProfile.C1.append(so.spID_new)
            elif so.locusName_new == CRName.C2:
                self.spProfile.C2.append(so.spID_new)
        # print('profile C1 ', self.spProfile.C1)
        # print('profile C2 ', self.spProfile.C2)
        # write the profile to 'spProfile.txt'
        spProfileFilePath = os.path.join(self.outFolderPath, 'spProfile.txt')
        util.writeSpacerProfile(spProfileFilePath, self.spProfile)

    def appendSpacerFastaTemp(self):
        """
        Append unmatched spacers in self.spacerOne_List to 'sptemp.fasta'
        Use only so.foundMatch = False
        :return:
        """
        spUnmatchedRecordList = []  # list of unmatched SpacerOne ins
        for so in self.spacerOne_Use:
            if so.foundMatch is False:
                record = SeqRecord(so.spSequence, id=so.spID_adj, description=so.desc)
                spUnmatchedRecordList.append(record)
        if len(spUnmatchedRecordList) != 0:
            # only write if there are at least one match found
            # not found, not writing
            newRecordList = []
            spTempFilePath = os.path.join(self.fileCls.SPACER_Folder, 'sptemp.fasta')
            for seq_record in SeqIO.parse(spTempFilePath, "fasta"):
                # print(type(seq_record))
                # print(seq_record.id)
                # print(seq_record.description)
                # print(seq_record.seq)
                newRecordList.append(seq_record)
            newRecordList = newRecordList + spUnmatchedRecordList
            # wrie to 'sptemp.fasta'
            with open(spTempFilePath, "w") as output_handle:
                SeqIO.write(newRecordList, output_handle, "fasta")

    def makeTemporarySpacerDB(self):
        """
        Need to create a temporary DB called \SPTEMP_DB\SPALL
        If the SPTEMP_DB folder already exists, we do not need create it anymore.
        This DB will be deleted after each batch run (all genome runs)
        And it will be deleted at 'batch.py'. Steps for the blast operation
        If this is the FIRST genome run encountered
          (1) Make a copy of 'spacerall.fasta' file to 'sptemp.fasta' (put in the same folder)
          (2) Blast 'spacerfound.fasta' against 'SPTEMP_DB' (created at batch.py)
              - e.g. >C1_0_BK_SAL1 BK_SAL1
          (3) Assign new names to all found spacers from the blast result
          (4) Create a spacer profile for this particular genome
          (5) Add the unmatched spacers to 'sptemp.fasta' (if there are)
          (6) Create SPTEMP_DB/SPALL database using the modified 'sptemp.fasta'
        If this is SECOND genome run encountered
          (1) Blast 'spacerfound.fasta' against 'SPTEMP_DB' (modified version)
              - e.g. >C1_0_BK_SAL1 BK_SAL1
          (2) Assign new names to all found spacers from the blast result
          (3) Create a spacer profile for this particular genome
          (4) Add the unmatched spacers to 'sptemp.fasta' (if there are)
          (5) Delete and create SPTEMP_DB/SPALL using the modified 'sptemp.fasta'
        After all genomes are run
          (1) Delete SPTEMP_DB
          (2) Delete 'sptemp.fasta'
        :return:
        """
        # get the list of scaffolds files in "in_file" folder
        spacerFilePath = os.path.join(self.fileCls.SPACER_Folder, 'sptemp.fasta')
        # create main db folder -> "/SPACER/SP_DB"
        spacerDB_name = 'SPTEMP_DB'
        # dbFolder -> "/SPACER/SP_DB"
        dbFolder = os.path.join(self.fileCls.SPACER_Folder, spacerDB_name)
        createFolderAndClear(dbFolder)
        # dbFolder -> "/SPACER/SP_DB/SPALL"
        dbFolder_X = os.path.join(dbFolder, 'SPALL')
        cmdLine = 'makeblastdb -in %s -dbtype nucl -parse_seqids -out %s' % (spacerFilePath,
                                                                             dbFolder_X,
                                                                             )
        os.system(cmdLine)

    def collectResult(self):
        """
        Collect result
        self.spUse_ANT -> {'C1': LocusANT,'C2': LocusANT}  final form
        :return: res -> RES_SeqSero2 obj
        """
        res = RES_CRISPR(self.outFolderName)
        # result of CRISPR (list of spacer?)
        # res.crisprRes = '|'.join(self.crisprRes)
        res.spacerC1 = '|'.join(self.spProfile.C1)
        res.spacerC2 = '|'.join(self.spProfile.C2)
        try:
            # self.spResANT_USE[CRName.C1] -> LocusUSE ins
            res.posStartC1 = self.spResANT_USE[CRName.C1].refStart
            res.posEndC1 = self.spResANT_USE[CRName.C1].refEnd
            res.DR_OneC1 = self.spResANT_USE[CRName.C1].DR_One
            res.NodeC1 = '|'.join(self.spResANT_USE[CRName.C1].NodeList)
            res.posStartC2 = self.spResANT_USE[CRName.C2].refStart
            res.posEndC2 = self.spResANT_USE[CRName.C2].refEnd
            res.DR_OneC2 = self.spResANT_USE[CRName.C2].DR_One
            res.NodeC2 = '|'.join(self.spResANT_USE[CRName.C2].NodeList)
        except (KeyError,):
            pass
        return res

    # def readSpacerContent(self):
    #     """
    #     Perform CrisprCas finder
    #     :return:
    #     """
    #     wb = load_workbook(File.Spacer_File)
    #     ws = wb['main_spacer']
    #     # d = list(ws.columns)  # [(<Cell 'A'.A1>, <Cell 'A'.A2>),(<Cell 'B'.A1>, <Cell 'B'.A2>)]
    #     self.spc_list = []
    #     d = list(ws.rows)
    #     for row in range(len(d)):
    #         if row == 0:
    #             continue
    #         v = d[row]
    #         spc = SpacerOne()
    #         spc.sid = str(v[0].value)
    #         spc.slength = int(v[1].value)
    #         spc.sequence = str(v[2].value)
    #         self.spc_list.append(spc)
    #         # break
    #     # for ss in self.spc_list:
    #     #     print('sid: ', ss.sid)
    #     #     print('slength: ', ss.slength)
    #     #     print('sequence: ', ss.sequence)

    # def transfer2FASTA(self):
    #     """
    #     Transfer to fasta file
    #     :return:
    #     """
    #     self.spcr_list = []  # list of SeqRecord ins
    #     for ss in self.spc_list:
    #         record = SeqRecord(Seq(ss.sequence, IUPAC.IUPACUnambiguousDNA),
    #                            id=ss.sid,
    #                            description="",
    #                            )
    #         self.spcr_list.append(record)
    #     # After getting seqToWrite, write to the file
    #     with open(File.SpacerFasta_File, "w") as output_handle:
    #         SeqIO.write(self.spcr_list, output_handle, "fasta")

    # def performSpacerFinderBlast(self, genomeName):
    #     """
    #     Find number of spacers in a genome
    #     E.G., if we want to find a spacer AbonB1 - GCTAGCCTGCTCCGCATTAACCGCCTTTAATG in the BK_SAL1 genome
    #     we want to perform local alignment using the direct sequence or its reverse complement sequence
    #
    #     :return:
    #     """
    #     # fastaFileName = 'SPI%s_gene.fasta' % i
    #     # fastaFile = File.SpacerFasta_FileTest
    #     fastaFile = File.SpacerFasta_File
    #     # dbFolder -> 'maindb/DBSC/BK_SAL1_DB
    #     dbFolder = os.path.join(File.DBSC_Folder, '%s_DB' % genomeName)
    #     # dbFolder_USE -> 'maindb/DBSC/BK_SAL1_DB/BK_SAL1
    #     dbFolder_USE = os.path.join(dbFolder, genomeName)
    #     # outputFolder -> 'maindb/out_blast_spc/BK_SAL1'
    #     outputFolder = os.path.join(File.outBlastSPCFolder, genomeName)
    #     createFolderAndClear(outputFolder)
    #     resultFilePath = os.path.join(outputFolder, 'SPC_blast.txt')
    #     cmdLine_1 = 'blastn -query %s -db %s -num_threads 4 -perc_identity 75 -qcov_hsp_perc 55 -out %s -outfmt "6 qacc sallseqid evalue bitscore pident qcovs qcovhsp qstart qend sstart send"' % (
    #         fastaFile,  # FASTA file
    #         dbFolder_USE,  # Database file
    #         resultFilePath,)  # Output file
    #     os.system(cmdLine_1)

    # def performRead_SpacerBlastRes(self, genomeName):
    #     """
    #     Read blast text result after performing method "performBlastN_CMD"
    #     determine if each scaffold contains an SPI or not
    #     The output format -> "qacc sallseqid evalue bitscore pident qcovs qcovhsp qstart qend sstart send"
    #     For instance, "ygbA	NODE_2_length_635729_cov_82.5337 1.24e-178 621 99.130 100 .. 345 545449 545793"
    #     qacc: Query accesion
    #     sallseqid: Subject Seq-id
    #     evalue: E-Value
    #     bitscore: Bitscore
    #     pident: Percentage of identical matches
    #     qcovs: Query Coverage Per Subject (for all HSPs)
    #     qcovhsp: Query Coverage Per HSP
    #     qstart: means Start of alignment in query
    #     qend: means End of alignment in query
    #     sstart: means Start of alignment in subject
    #     send: means End of alignment in subject
    #     :param genomeName: name of a genome e.g. BK_SAL1
    #     :return:
    #     """
    #     # access the blast result
    #     blastFileName = 'SPC_blast.txt'
    #     # sfBlastFolder -> 'out_blast_spc/BK_SAL1'
    #     spcBlastFolder = os.path.join(File.outBlastSPCFolder, genomeName)
    #     # spiBlastFile -> 'out_blast_spc/BK_SAL1/SPC_blast.txt'
    #     spcBlastFile = os.path.join(spcBlastFolder, blastFileName)
    #     xBlast = []
    #     # read info from spcBlastFile
    #     r = read_file_dlimit(spcBlastFile)
    #     for x in r:
    #         RB = RES_BLAST(genomeName, None)
    #         RB.qacc = x[0]
    #         RB.sallseqid = x[1]
    #         RB.evalue = x[2]
    #         RB.bitscore = x[3]
    #         RB.pident = x[4]
    #         RB.qcovs = x[5]
    #         RB.qcovhsp = x[6]
    #         RB.qstart = x[7]
    #         RB.qend = x[8]
    #         RB.sstart = x[9]
    #         RB.send = x[10]
    #         if float(RB.evalue) <= Parameter.SPC.E_VALUE_THRESH:
    #             # qualified hits
    #             xBlast.append(RB)
    #     print('numSPC: ', len(xBlast))

    # def createSpacerInitialDatabase(self):
    #     """
    #     Put spacers data into databasse by using Blast using spacer1.xlsx
    #     And run suffix numbers, for examples
    #     Suppose there are 5 spacers with names
    #     (1) Aba1 (2) Aba2 (3) Aba3 (4) Aba4 (5) Aba5
    #     Two cases that might happen
    #     1. If the target spacer matches perfectly with Aba3, it will be assigned the name "Aba3"
    #     2. If the target spacer matches closely with Aba3 (but not 100% identity and coverage, by looking at the blast
    #        result), it will be assigned the name "Aba6"
    #     Main Steps
    #     (1) Create database from the file "spacerall.fasta"
    #     (2) Use the blast tool to blast each spacer to such database
    #     (3) Read and interpret the blast result
    #     :return:
    #     """
    #     spacerFastaFile = File.SpacerFasta_File  # 'spacerall.fasta'
    #     # create main db folder for spacer -> "/maindb/SPACER/TEMPSP"
    #     tempDB_Folder = 'TEMPSP'
    #     tempDB_FolderPath = os.path.join(File.SPACER_Folder, tempDB_Folder)
    #     createFolderAndClear(tempDB_FolderPath)
    #     tempDB_FolderPath_X = os.path.join(tempDB_FolderPath, 'SPALL')
    #     cmdLine = 'makeblastdb -in %s -dbtype nucl -parse_seqids -out %s' % (spacerFastaFile,
    #                                                                          tempDB_FolderPath_X,
    #                                                                          )
    #     os.system(cmdLine)

    # def performSpacerGenomeBlast(self):
    #     """
    #     Blast all spacers in the entire genome
    #     Use performSingleSpacerBlast method for a single spacer blast
    #     :return:
    #     """

    # def performSingleSpacerBlast(self, spacer):
    #     """
    #     Perform blastn using a single spacer to the database of the spacer at "/maindb/SPACER/TEMPSP"
    #     The objective is to determine which spacer in database is closest to the target spacer. And must get the
    #     name of that close spacer
    #     :param spacer: target spacer sequence
    #     :return:
    #     """
    #     # make the spacer into a fasta file (temporary one)
    #     spacerList = []
    #     record = SeqRecord(
    #         spacer,
    #         id='SP1',
    #         # name=geneName,
    #         description='Single Spacer',
    #     )
    #     spacerList.append(record)
    #     # create filePath of the fasta file (spacer sequence)
    #     fastaFilePath = 5
    #
    #     with open(fastaFolder, "w") as output_handle:
    #         SeqIO.write(geneList, output_handle, "fasta")

    # geneList.append(record)
    # fastaFileName = 'SPI%s_gene.fasta' % i
    # fastaFile = os.path.join(File.SPIFT_Folder, fastaFileName)
    # dbFolder = os.path.join(File.DBSC_Folder, '%s_DB' % genomeName)
    # # dbFolder_USE -> 'maindb/DBSC/BK_SAL1_DB/BK_SAL1
    # dbFolder_USE = os.path.join(dbFolder, genomeName)
    # # outputFolder -> 'maindb/out_blast/BK_SAL1'
    # outputFolder = os.path.join(File.outBlastFolder, genomeName)
    # createFolderAndClear(outputFolder)
    # resultFilePath = os.path.join(outputFolder, 'SPI%s_blast.txt' % i)
    # cmdLine_1 = 'blastn -query %s -db %s -num_threads 4 -perc_identity 75 -qcov_hsp_perc 55 -out %s -outfmt "6 qacc sallseqid evalue bitscore pident qcovs qcovhsp qstart qend sstart send"' % (
    # fastaFile,  # FASTA file
    # dbFolder_USE,  # Database file
    # resultFilePath,)  # Output file
    # os.system(cmdLine_1)
    # break

    # def performReadNCBI_BlastResCMD(self, genomeName):
    #     """
    #     Read blast text result after performing method "performBlastN_CMD"
    #     determine if each scaffold contains an SPI or not
    #     The output format -> "qacc sallseqid evalue bitscore pident qcovs qcovhsp qstart qend sstart send"
    #     For instance, "ygbA	NODE_2_length_635729_cov_82.5337 1.24e-178 621 99.130 100 .. 345 545449 545793"
    #     qacc: Query accesion
    #     sallseqid: Subject Seq-id
    #     evalue: E-Value
    #     bitscore: Bitscore
    #     pident: Percentage of identical matches
    #     qcovs: Query Coverage Per Subject (for all HSPs)
    #     qcovhsp: Query Coverage Per HSP
    #     qstart: means Start of alignment in query
    #     qend: means End of alignment in query
    #     sstart: means Start of alignment in subject
    #     send: means End of alignment in subject
    #
    #     The idea is that if all the genes in the FASTA file (for each SPI) can be found in the blast result,
    #     we would be sure that the scaffold contains that SPI.
    #     :param genomeName: name of a genome e.g. BK_SAL1
    #     :return:
    #     """
    #     startSPI_N = 1
    #     stopSPI_N = 14
    #     for i in range(startSPI_N, stopSPI_N + 1):
    #         # find list of genes from SPI1_gene
    #         # access the file location of the FASTA file
    #         geneListFASTA = []  # gene list in the FASTA file (all genes)
    #         matchList = []  # [('R.068P_0','OTU_29'),...]
    #         fastaFileName = 'SPI%s_gene.fasta' % i
    #         fastaFile = os.path.join(File.SPIFT_Folder, fastaFileName)
    #         for seq_record in SeqIO.parse(fastaFile, "fasta"):
    #             geneListFASTA.append(seq_record.id)
    #         print('geneListFasta: ', geneListFASTA)
    #         # access the blast result
    #         blastFileName = 'SPI%s_blast.txt' % i
    #         # sfBlastFolder -> 'out_blast/BK_SAL1'
    #         sfBlastFolder = os.path.join(File.outBlastFolder, genomeName)
    #         # spiBlastFile -> 'out_blast/BK_SAL1/SPI1_blast.txt'
    #         spiBlastFile = os.path.join(sfBlastFolder, blastFileName)
    #         xBlast = []
    #         # read info from spiBlastFile
    #         r = read_file_dlimit(spiBlastFile)
    #         # print(r)
    #         for x in r:
    #             RB = RES_BLAST(genomeName, fastaFileName)
    #             RB.qacc = x[0]
    #             RB.sallseqid = x[1]
    #             RB.evalue = x[2]
    #             RB.bitscore = x[3]
    #             RB.pident = x[4]
    #             RB.qcovs = x[5]
    #             RB.qcovhsp = x[6]
    #             RB.qstart = x[7]
    #             RB.qend = x[8]
    #             RB.sstart = x[9]
    #             RB.send = x[10]
    #             if float(RB.evalue) <= Parameter.SPI.E_VALUE_THRESH:
    #                 # qualified hits
    #                 xBlast.append(RB)
    #         # print(xBlast)
    #         # take only .qacc
    #         xy = [u.qacc for u in xBlast]
    #         # verify if xy contains all the genes in geneListFASTA or not
    #         allMatch = False
    #         numGene = len(geneListFASTA)
    #         numIns = len(list(set(xy).intersection(set(geneListFASTA))))
    #         if numGene == numIns:
    #             allMatch = True
    #         print('numGene: ', numGene)
    #         print('numIns: ', numIns)
    #         print('allMatch: ', allMatch)
    #         break

    # def performSpacerFinderLocal(self, genomeName):
    #     """
    #     Find number of spacers in a genome
    #     E.G., if we want to find a spacer AbonB1 - GCTAGCCTGCTCCGCATTAACCGCCTTTAATG in the BK_SAL1 genome
    #     we want to perform local alignment using the direct sequence or its reverse complement sequence
    #
    #     :return:
    #     """
    #     # Genome name is superb, we must define genome name
    #     score = pairwise2.align.localxx(seq_record.seq.upper(),
    #                                     sOTU.seq.upper(),
    #                                     score_only=True,
    #                                     one_alignment_only=True)
    #     identity = round((score / seqLen) * 100, 2)


if __name__ == '__main__':
    obj = CrisprCas_CLS(None, None, None, None)
    # obj.readSpacerContent()
    # obj.transfer2FASTA()
    # obj.performSpacerFinderBlast('BK_SAL1')
    # obj.performRead_SpacerBlastRes('BK_SAL1')
    # obj.testCrisprTool()
    # obj.readCrisprReport()
    # obj.createSpacerInitialDatabase()
