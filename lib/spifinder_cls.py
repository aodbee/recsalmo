# Copyrights 2022: Nuttachat Wisittipanit

from const import (File, SPIFinder_CONST, Parameter, )
from structobj import (RES_SPI, RES_BLAST,)
from util import read_file_dlimit
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
import os
import shutil
from util import createFolderAndClear
from openpyxl import load_workbook

"""
amrfinder --nucleotide fasta_file --organism Salmonella -o output_file 

"""


class SPIFinder_CLS:

    def __init__(self, fileCls, inputFilePath, outFolderPath, outFolderName):
        self.fileCls = fileCls
        self.inputFilePath = inputFilePath  # in_file/BK_SAL1.scaffolds.fasta
        self.outFolderPath = outFolderPath  # out_file/BK_SAL1
        self.outFolderName = outFolderName  # BK_SAL1
        self.SPI_foundList = None

    def perform(self):
        """
        Perform main SPI finder run here, steps include
        (1) Get a total of 17 SPI sequences ready (SPI1-SPI17)
        (2) For each SPI sequence, use that as a query to do blast against each scaffold db
        (3) If the blast result shows identity > 90% and coverage > 60%, that SPI is found in the genome.
        (4) Report all the SPIs found in SPIFinder_CONST.mainResultFile
        :return:
        """
        # self.file_input -> 'BK_CAL1.scaffolds.fasta'
        # outFilePath = os.path.join(self.outFolderPath, SPIFinder_CONST.mainResultFile)
        self.performBlastN_CMDLine(self.outFolderName)
        self.SPI_foundList = self.performReadNCBI_BlastResCMD(self.outFolderName)

    def collectResult(self):
        """
        Collect result
        :return: res -> RES_SeqSero2 obj
        """
        res = RES_SPI(self.outFolderName)
        res.spiList = '|'.join(self.SPI_foundList)
        return res

    def performBlastN_CMDLine(self, genomeName):
        """
        Perform blastn using "SPI_combine.fasta" as the query to the scaffold DB
        Determine which SPI a scaffold file contains. E.g. We want to know whether 'BK_SAL1.scaffolds.fasta' contains
        any SPI or not (SPI-1 to SPI-17)
        This means WE NEED TO perform blast operations for the FASTA file containing SPI-1 to SPI-17.
        The blast results are in 'out_blast/genomeName/' folder. That folder contains 'SPI_finderResult.txt'
        :param genomeName: genome name e.g. 'BK_SAL1'
        :return:
        """
        fastaFile = 'SPI_combine.fasta'
        fastaFilePath = os.path.join(self.fileCls.SPIFS_Folder, fastaFile)
        # dbFolderPath -> 'maindb/DBSC/BK_SAL1_DB
        dbFolderPath = os.path.join(self.fileCls.DBSC_Folder, '%s_DB' % genomeName)
        # dbFolderUsePath -> 'maindb/DBSC/BK_SAL1_DB/BK_SAL1
        dbFolderUsePath = os.path.join(dbFolderPath, genomeName)
        # outputFolderPath -> 'out_file/BK_SAL1'
        outputFolderPath = os.path.join(self.fileCls.outputFolder, genomeName)
        # createFolderAndClear(outputFolderPath)
        resultFilePath = os.path.join(outputFolderPath, 'SPI_finderResult.txt')
        # SPI1	NODE_2_length_635729_cov_82.5337	0.0	60327	98.496	82	82	4223	38434	540680	506466
        cmdLine = 'blastn -query %s -db %s -num_threads 4 -perc_identity 90 -qcov_hsp_perc 60 -out %s -outfmt "6 qacc sallseqid evalue bitscore pident qcovs qcovhsp qstart qend sstart send"' % (
            fastaFilePath,  # FASTA file
            dbFolderUsePath,  # Database file
            resultFilePath,)  # Output file
        os.system(cmdLine)

    def performReadNCBI_BlastResCMD(self, genomeName):
        """
        Read blast text result after performing method "performBlastN_CMD"
        determine if each scaffold contains an SPI or not
        The output format -> "qacc sallseqid evalue bitscore pident qcovs qcovhsp qstart qend sstart send"
        For instance, "ygbA	NODE_2_length_635729_cov_82.5337 1.24e-178 621 99.130 82 82	4223 38434 540680 506466"
        qacc: Query accesion -> 'ygbA'
        sallseqid: Subject Seq-id -> 'NODE_2_length_635729_cov_82.5337'
        evalue: E-Value -> '1.24e-178'
        bitscore: Bitscore -> '621'
        pident: Percentage of identical matches -> '99.130'
        qcovs: Query Coverage Per Subject (for all HSPs) -> '82'
        qcovhsp: Query Coverage Per HSP -> '82'
        qstart: means Start of alignment in query -> '4223'
        qend: means End of alignment in query -> '38434'
        sstart: means Start of alignment in subject -> '540680'
        send: means End of alignment in subject -> '506466'

        :param genomeName: name of a genome e.g. BK_SAL1
        :return:
        """
        # access the blast result
        blastResultFile = 'SPI_finderResult.txt'
        # sfBlastFolder -> 'out_file/BK_SAL1'
        sfBlastFolder = os.path.join(self.fileCls.outputFolder, genomeName)
        # spiBlastFile -> 'out_blast/BK_SAL1/SPI1_blast.txt'
        blastResultFilePath = os.path.join(sfBlastFolder, blastResultFile)
        xBlast = []
        # read info from spiBlastFile
        r = read_file_dlimit(blastResultFilePath)
        # print(r)
        for x in r:
            RB = RES_BLAST(genomeName)
            RB.qacc = x[0]
            RB.sallseqid = x[1]
            RB.evalue = x[2]
            RB.bitscore = x[3]
            RB.pident = x[4]
            RB.qcovs = x[5]
            RB.qcovhsp = x[6]
            RB.qstart = x[7]
            RB.qend = x[8]
            RB.sstart = x[9]
            RB.send = x[10]
            xBlast.append(RB)
            # if float(RB.evalue) <= Parameter.SPI.E_VALUE_THRESH:
            #     # qualified hits
            #     xBlast.append(RB)
        # print('xBlast')
        # print(xBlast)
        # take only .qacc
        SPI_foundList = [u.qacc for u in xBlast]
        # print('SPI_foundList')
        # print(SPI_foundList)
        return SPI_foundList

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

    # def transformGB2Fasta(self):
    #     """
    #     Parsing GenBank file to extract its sequence
    #     :return:
    #     """
    #     startNum = 1
    #     stopNum = 12
    #     for i in range(startNum, stopNum + 1):
    #         mainName = 'SPI%s' % i
    #         genBankFile = '%s_gb.fasta' % mainName  # file in Genbank record format
    #         genBankFilePath = os.path.join(File.SPIFT_Folder, genBankFile)
    #         # fastaFile -> 'SPI1_gene.fasta'
    #         fastaFile = '%s_complete.fasta' % mainName  # SPI in FASTA (complete sequence)
    #         fastaFilePath = os.path.join(File.SPIFS_Folder, fastaFile)
    #         srecList = []  # list of SeqRecord
    #         for seq_record in SeqIO.parse(genBankFilePath, "genbank"):
    #             # seq_record.seq -> Seq obj
    #             seqUse = seq_record.seq.upper()
    #             record = SeqRecord(
    #                 seqUse,
    #                 id=mainName,
    #                 # name=geneName,
    #                 description='complete',
    #             )
    #             srecList.append(record)
    #
    #         with open(fastaFilePath, "w") as output_handle:
    #             SeqIO.write(srecList, output_handle, "fasta")
    #         # break

    # def transformSPI13GBtoFasta(self):
    #     """
    #     Parsing GenBank file to extract its sequence
    #     :return:
    #     """
    #     mainName = 'SPI13'
    #     genBankFile = '%s.gb' % mainName  # file in Genbank record format
    #     genBankFilePath = os.path.join(File.SPIGB_Folder, genBankFile)
    #     fastaFile = '%s_complete.fasta' % mainName  # SPI in FASTA (complete sequence)
    #     fastaFilePath = os.path.join(File.SPIFS_Folder, fastaFile)
    #     srecList = []  # list of SeqRecord
    #     for seq_record in SeqIO.parse(genBankFilePath, "fasta"):
    #         # seq_record.seq -> Seq obj
    #         seqUse = seq_record.seq.upper()
    #         seqUse = seqUse[3141722-1:3165661]
    #         record = SeqRecord(
    #             seqUse,
    #             id=mainName,
    #             # name=geneName,
    #             description='complete',
    #         )
    #         srecList.append(record)
    #     with open(fastaFilePath, "w") as output_handle:
    #         SeqIO.write(srecList, output_handle, "fasta")

    # def transformSPI14GBtoFasta(self):
    #     """
    #     Parsing GenBank file to extract its sequence
    #     :return:
    #     """
    #     mainName = 'SPI14'
    #     genBankFile = 'LT2.gb'  # file in Genbank record format
    #     genBankFilePath = os.path.join(File.REFDBFolder, genBankFile)
    #     fastaFile = '%s_complete.fasta' % mainName  # SPI in FASTA (complete sequence)
    #     fastaFilePath = os.path.join(File.SPIFS_Folder, fastaFile)
    #     srecList = []  # list of SeqRecord
    #     for seq_record in SeqIO.parse(genBankFilePath, "genbank"):
    #         # seq_record.seq -> Seq obj
    #         seqUse = seq_record.seq.upper()
    #         # print(seqUse)
    #         seqUse = seqUse[926180-1:932916]
    #         # print(seqUse)
    #         record = SeqRecord(
    #             seqUse,
    #             id=mainName,
    #             # name=geneName,
    #             description='complete',
    #         )
    #         srecList.append(record)
    #     with open(fastaFilePath, "w") as output_handle:
    #         SeqIO.write(srecList, output_handle, "fasta")

    # def transformSPIGBtoFasta(self):
    #     """
    #     Transform SPI locations
    #     :param mainname:
    #     :param gbFilePath:
    #     :param Locations:
    #     :return:
    #     """
    #     mainName = ['SPI15','SPI16','SPI17']
    #     loc = [(3053654,3060017),(605515,609992),(2460793,2465914)]
    #     genBankFile = 'refgenome_CT18.fasta'  # file in Genbank record format
    #     genBankFilePath = os.path.join(File.REFDBFolder, genBankFile)
    #     for i in range(len(mainName)):
    #         fastaFile = '%s_complete.fasta' % mainName[i]  # SPI in FASTA (complete sequence)
    #         fastaFilePath = os.path.join(File.SPIFS_Folder, fastaFile)
    #         srecList = []  # list of SeqRecord
    #         for seq_record in SeqIO.parse(genBankFilePath, "fasta"):
    #             # seq_record.seq -> Seq obj
    #             seqUse = seq_record.seq.upper()
    #             # print(seqUse)
    #             seqUse = seqUse[loc[i][0]-1:loc[i][1]]
    #             # print(seqUse)
    #             record = SeqRecord(
    #                 seqUse,
    #                 id=mainName[i],
    #                 # name=geneName,
    #                 description='complete',
    #             )
    #             srecList.append(record)
    #         with open(fastaFilePath, "w") as output_handle:
    #             SeqIO.write(srecList, output_handle, "fasta")

    # def combineSPItoFasta(self):
    #     """
    #     Transform SPI locations
    #     :param mainname:
    #     :param gbFilePath:
    #     :param Locations:
    #     :return:
    #     """
    #     start = 1
    #     stop = 17
    #     # mainName = ['SPI15','SPI16','SPI17']
    #     # genBankFile = 'refgenome_CT18.fasta'  # file in Genbank record format
    #     # SPIsFilePath = os.path.join(File.SPIFS_Folder, genBankFile)
    #     allSPI = []
    #     for i in range(start, stop+1):
    #         fastaFile = 'SPI%s_complete.fasta' % i  # SPI in FASTA (complete sequence)
    #         fastaFilePath = os.path.join(File.SPIFS_Folder, fastaFile)
    #         # srecList = []  # list of SeqRecord
    #         for seq_record in SeqIO.parse(fastaFilePath, "fasta"):
    #             # # seq_record.seq -> Seq obj
    #             # seqUse = seq_record.seq.upper()
    #             # # print(seqUse)
    #             # seqUse = seqUse[loc[i][0]-1:loc[i][1]]
    #             # # print(seqUse)
    #             # record = SeqRecord(
    #             #     seqUse,
    #             #     id=mainName[i],
    #             #     # name=geneName,
    #             #     description='complete',
    #             # )
    #             allSPI.append(seq_record)
    #     combineSPIFilePath = os.path.join(File.SPIFS_Folder, 'SPI_combine.fasta')
    #     with open(combineSPIFilePath, "w") as output_handle:
    #         SeqIO.write(allSPI, output_handle, "fasta")

    # def runParseGenbankBatch(self):
    #     """
    #     Run 'parseGenBank' method in batch
    #     :return:
    #     """
    #     startNum = 6
    #     stopNum = 12
    #     for i in range(startNum, stopNum + 1):
    #         genBankName = 'SPI%s' % i
    #         self.parseGenBank(genBankName)
    #         # break
    #
    # def parseGenBank(self, genBankName):
    #     """
    #     Simple parsing of genbank file (Transform Genbank into FASTA file)
    #     E.g. genBankName -> 'SPI1', FASTA file -> 'SPI1_gene.fasta'
    #     :param genBankName: 'SPIX' e.g. 'SPI1'
    #     :return:
    #     """
    #     print('genBankName: ', genBankName)
    #     # genBankFile -> 'SPI1.gb'
    #     genBankNameFile = '%s.gb' % genBankName  # file in Genbank record format
    #     genBankNameFolder = os.path.join(File.SPIGB_Folder, genBankNameFile)
    #     # fastaFile -> 'SPI1_gene.fasta'
    #     fastaFile = '%s_gene.fasta' % genBankName  # Genbank record in FASTA (only sequence)
    #     fastaFolder = os.path.join(File.SPIFT_Folder, fastaFile)
    #     geneList = []  # list of SeqRecord (all genes)
    #     for seq_record in SeqIO.parse(genBankNameFolder, "genbank"):
    #         # seq_record.seq -> Seq obj
    #         # print(seq_record.id)
    #         # print(repr(seq_record.seq))
    #         # print(len(seq_record))
    #         # print(seq_record.features)
    #         for seqFeature in seq_record.features:
    #             # loop for each "gene" feature in the genbank record
    #             if seqFeature.type != 'gene':
    #                 # if seqFeature.type != 'CDS':
    #                 continue
    #             # seqFeature.type -> gene
    #             # seqFeature.location -> [ExactPosition(0):ExactPosition(345)](-)
    #             # seqFeature.location.start -> ExactPosition(0)
    #             # seqFeature.location.end -> ExactPosition(345)
    #             # seqFeature.location.strand -> -1
    #             # seqFeature.qualifiers -> OrderedDict([('gene', ['ygbA']), ('locus_tag', ['SC2793'])])
    #             # seqFeature.id -> <unknown id>
    #             # seqFeature.qualifiers['gene'][0] -> "ygbA"
    #             # print('type: ', seqFeature.type)  # should be gene
    #             # print('location: ', seqFeature.location)
    #             # print('location type: ', type(seqFeature.location))
    #             # print('location start: ', seqFeature.location.start)
    #             # print('location end: ', seqFeature.location.end)
    #             # print('location strand: ', seqFeature.location.strand)
    #             # print('qualifiers: ', seqFeature.qualifiers)
    #             # print('id: ', seqFeature.id)
    #             startLoc = seqFeature.location.start
    #             endLoc = seqFeature.location.end
    #             seqUse = seq_record.seq[startLoc:endLoc]
    #             try:
    #                 geneName = seqFeature.qualifiers['gene'][0]
    #                 print('geneName: ', geneName)
    #                 # print('geneName type: ', type(geneName))
    #                 # seq_record.id -> NC_006905_P5.1
    #                 desc = 'SPI1-%s' % seq_record.id
    #                 print('desc: ', desc)
    #                 # print(seqFeature.__dict__.keys())
    #                 # for a in seqFeature.__dict__:
    #                 record = SeqRecord(
    #                     seqUse,
    #                     id=geneName,
    #                     # name=geneName,
    #                     description=desc,
    #                 )
    #                 geneList.append(record)
    #             except (KeyError,):
    #                 pass
    #             # break
    #     # turn into fasta file
    #     # count = SeqIO.convert("cor6_6.gb", "genbank", "cor6_6.fasta", "fasta")
    #     # count = SeqIO.convert(genBankFile, "genbank", fastaFile, "fasta")
    #     # print("Converted %i records" % count)
    #     # write gene into fasta file
    #     with open(fastaFolder, "w") as output_handle:
    #         SeqIO.write(geneList, output_handle, "fasta")


if __name__ == '__main__':
    obj = SPIFinder_CLS(None, None, None)
    # obj.perform()
    # obj.performBlastN_CMDLine('BK_SAL1')
    # obj.transformGB2Fasta()
    # obj.transformSPI13GBtoFasta()
    # obj.transformSPI14GBtoFasta()
    # obj.runParseGenbankBatch()
    # obj.transformSPIGBtoFasta()
    # obj.combineSPItoFasta()
    # obj.parseGenBank()
    # obj.makeBlastDB('BK_SAL1')
    # obj.performBlastN()
    # obj.performBlastN_CMDLine('BK_SAL1')
    obj.performReadNCBI_BlastResCMD('BK_SAL1')
    # obj.performReadNCBI_BlastResCMD()
