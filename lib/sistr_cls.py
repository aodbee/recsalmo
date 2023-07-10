# Copyrights 2023: Nuttachat Wisittipanit


from const import (File, SeqSero2_CONST, )
from structobj import RES_SISTR
from util import read_file_normal
from sistr.sistr_cmd import sistr_predict
import os


class SistrCmdMockArgs:
    run_mash = True
    no_cgmlst = False
    qc = True
    use_full_cgmlst_db = False


class SISTR_CLS:

    def __init__(self, fileCls, inputFilePath, outFolderPath, outFolderName):
        self.fileCls = fileCls
        self.inputFilePath = inputFilePath
        self.outFolderPath = outFolderPath
        self.outFolderName = outFolderName  # BK_SAL1
        self.sistr_results = None
        self.allele_results = None

    def perform(self):
        """
        Perform SeqSero2 run
        :return:
        """
        # self.file_input -> 'BK_CAL1.scaffolds.fasta'
        # SeqSero2_package.py -m k -t 4 -i assembly.fasta -d outFolder
        # run SISTR serovar prediction
        genome_fasta_path = self.inputFilePath
        genome_name = self.outFolderName
        tempDir = os.path.join(self.fileCls.outputFolder, '%s_%s' % (genome_name, 'tmpSistr',))
        self.sistr_results, self.allele_results = sistr_predict(genome_fasta_path,
                                                                genome_name,
                                                                keep_tmp=False,
                                                                tmp_dir=tempDir,
                                                                args=SistrCmdMockArgs)
        # print(self.sistr_results.serovar_cgmlst)
        # print(self.sistr_results.cgmlst_matching_alleles)
        # print(self.sistr_results.h1)
        # print(self.sistr_results.h2)
        # print(self.sistr_results.serovar_antigen)
        # print(self.sistr_results.cgmlst_distance)
        # print(self.sistr_results.cgmlst_genome_match)
        # print(self.sistr_results.cgmlst_ST)
        # print(self.sistr_results.serovar)
        # print(self.sistr_results.serogroup)
        # print(self.sistr_results.cgmlst_subspecies)
        # self.resSISTR.serovar_cgmlst = sistr_results.serovar_cgmlst
        # self.resSISTR.cgmlst_matching_alleles = sistr_results.cgmlst_matching_alleles
        # self.resSISTR.h1 = sistr_results.h1
        # self.resSISTR.h2 = sistr_results.h2
        # self.resSISTR.serovar_antigen = sistr_results.serovar_antigen
        # self.resSISTR.cgmlst_distance = sistr_results.cgmlst_distance
        # self.resSISTR.cgmlst_genome_match = sistr_results.cgmlst_genome_match
        # self.resSISTR.cgmlst_ST = sistr_results.cgmlst_ST
        # self.resSISTR.serovar = sistr_results.serovar
        # self.resSISTR.serogroup = sistr_results.serogroup
        # self.resSISTR.cgmlst_subspecies = sistr_results.cgmlst_subspecies
        # print(allele_results)

    def collectResult(self):
        """
        Collect result
        :return: res -> RES_SeqSero2 obj
        """
        res = RES_SISTR(self.outFolderName)
        res.serovar_cgmlst = self.sistr_results.serovar_cgmlst
        res.cgmlst_matching_alleles = self.sistr_results.cgmlst_matching_alleles
        res.h1 = self.sistr_results.h1
        res.h2 = self.sistr_results.h2
        res.serovar_antigen = self.sistr_results.serovar_antigen
        res.cgmlst_distance = self.sistr_results.cgmlst_distance
        res.cgmlst_genome_match = self.sistr_results.cgmlst_genome_match
        res.cgmlst_ST = self.sistr_results.cgmlst_ST
        res.serovar = self.sistr_results.serovar
        res.serogroup = self.sistr_results.serogroup
        res.cgmlst_subspecies = self.sistr_results.cgmlst_subspecies
        return res
