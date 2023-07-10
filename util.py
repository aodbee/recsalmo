from copy import deepcopy
import os
import re
import shutil
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from structobj import (RES_BLAST, SpacerProfile, )
from const import File


def read_file_text(file_name):
    """
    Read pure file
    :param file_name:
    :return:
    """
    data = None
    with open(file_name, 'r') as f:
        data = f.read()
    return data


def read_file_normal(file_name):
    """
    Read file and delete empty lines
    :param file_name:
    :return:
    """
    a = []
    with open(file_name, 'r') as f:
        data = f.read()
        lines = data.split('\n')
        for line in lines:
            if line == '':
                continue
            a.append(line)
    return a


def read_file_dlimit(file_name):
    """
    Read file including spliting the lines
    :param file_name:
    :return:
    """
    a = []
    with open(file_name, 'r') as f:
        data = f.read()
        lines = data.split('\n')
        for line in lines:
            if line == '':
                continue
            s = line.split()
            a.append(s)
    return a


def read_file_tab(file_name, headInclude = False):
    """
    Read file including spliting the lines
    :param file_name:
    :return:
    """
    a = []
    with open(file_name, 'r') as f:
        data = f.read()
        lines = data.split('\n')
        for i, line in enumerate(lines):
            if line == '':
                continue
            if i == 0 and headInclude is True:
                continue
            s = line.split('\t')
            a.append(s)
    return a


def writeSpacerProfile(file_name, spProfile):
    """
    Write spacer profile to a file
    :param file_name: file to write the profile to
    :param spProfile: SpacerProfile ins
    :return:
    """
    C1_profile = ' '.join(spProfile.C1)
    C2_profile = ' '.join(spProfile.C2)
    with open(file_name, 'w') as f:
        f.write(spProfile.genomeName)
        f.write('\n')
        f.write(C1_profile)
        f.write('\n')
        f.write(C2_profile)


def readSpacerProfile(file_name):
    """
    Read 'spProfile.txt' from each genome and return
    :param file_name:
    :return:
    """
    a = read_file_dlimit(file_name)
    spProfile = SpacerProfile()
    spProfile.genomeName = a[0][0]
    try:
        spProfile.C1 = a[1]
    except (IndexError,):
        pass
    try:
        spProfile.C2 = a[2]
    except (IndexError,):
        pass
    return spProfile


def readBlastResult(source):
    """
    read and return the blast result
    :param source: blast result in txt file
    :return:
    """
    xBlast = []
    # read info from spiBlastFile
    r = read_file_dlimit(source)
    # print(r)
    for x in r:
        RB = RES_BLAST()
        RB.qacc = x[0]
        # modify sallseqid -> sometimes it is in the pattern ref|xxx| where xxx is the node name
        # sjid -> 'ref|NZ_CP050731.1|' or 'NODE_1'
        sjid = x[1]
        # findall -> ['NZ_CP050731.1']
        findall = re.findall(r'ref\|(.*)\|', sjid)
        try:
            match = findall[0]
        except (IndexError,):
            match = sjid  # in case regx cannot find anything
        RB.sallseqid = match
        RB.evalue = x[2]
        RB.bitscore = x[3]
        RB.pident = x[4]
        RB.qcovs = x[5]
        RB.qcovhsp = x[6]
        RB.qstart = int(x[7])
        RB.qend = int(x[8])
        RB.sstart = int(x[9])
        RB.send = int(x[10])
        RB.sstrand = x[11]
        xBlast.append(RB)
    return xBlast


def checkCopySpacerFasta(fileCls):
    """
    Check and copy 'spacerall.fasta'
    return True if 'spacerall.fasta' already exists
           False if 'spacerall.fasta' does not exist
    :return:
    """
    spFile = 'spacerall.fasta'  # Spacer file
    spFilePath = os.path.join(fileCls.SPACER_Folder, spFile)
    # check if 'sptemp.fasta' file exists or not
    spTempFile = 'sptemp.fasta'
    spTempFilePath = os.path.join(fileCls.SPACER_Folder, spTempFile)
    spTempFileExist = os.path.exists(spTempFilePath)
    # if 'sptemp.fasta' does not exist, copy from 'spacerall.fasta'
    if spTempFileExist is False:
        shutil.copyfile(spFilePath, spTempFilePath)
        return False
    return True


def createFolderAndClear(folderPath):
    """
    Create folderPath (if not existed or do nothing if already existed) and also clear everything inside the folder
    :param folderPath:
    :return:
    """
    try:
        os.mkdir(folderPath)
    except (FileExistsError,):
        # folder already exists
        pass
    # remove all contents in the folder first
    for filename in os.listdir(folderPath):
        file_path = os.path.join(folderPath, filename)
        try:
            if os.path.isfile(file_path) or os.path.islink(file_path):
                os.unlink(file_path)
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)
        except Exception as e:
            # print('Failed to delete %s. Reason: %s' % (file_path, e))
            pass


def createFolderNotExist(folderPath):
    """
    Create folderPath (if not existed or do nothing if already existed) and also clear everything inside the folder
    :param folderPath:
    :return:
    """
    try:
        os.mkdir(folderPath)
    except (FileExistsError,):
        # folder already exists, do absolutely nothing
        pass


def removeFolderAndContent(folderPath):
    """
    Remove the folder and all of its content
    :param folderPath: folder to be removed
    :return:
    """
    try:
        shutil.rmtree(folderPath)
    except (FileNotFoundError,):
        pass


def removeFile(filePath):
    """
    Remove the file
    :param filePath: path to the file
    :return:
    """
    if os.path.isfile(filePath):
        os.remove(filePath)
    else:
        # If it fails, inform the user.
        print("Error: %s file not found" % filePath)


def getSubFasta(fastaInput, fastaOutput, seq_id, seq_desc, start, stop):
    """
    Obtain part sequence of the fastaInput
    :param fastaInput: input file in FASTA format
    :param fastaOutput: output file in FASTA format
    :param seq_id: sequence id of the output
    :param seq_desc: sequence description of the output
    :param start: start position of the input sequence
    :param stop: stop position of the input sequence
    :return:
    """
    seqList = []
    with open(fastaInput) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            # print(record.id)
            # print(record.name)
            # print(record.description)
            # print(record.features)
            # print(record.seq)
            seqUse = record.seq[start - 1:stop]
            record = SeqRecord(
                seqUse,
                id=seq_id,
                description=seq_desc,
            )
            print('Length of the sequence: ', len(record.seq))
            seqList.append(record)
            break
    with open(fastaOutput, "w") as output_handle:
        SeqIO.write(seqList, output_handle, "fasta")


if __name__ == '__main__':
    fastaInput = '/home/nuttachat/Python/Project/WGS/ref_genome/refgenome_gallinarum28791.fasta'
    fastaOutput = '/home/nuttachat/Python/Project/WGS/maindb/SPI_FT/SPI13_fasta.fasta'
    getSubFasta(fastaInput, fastaOutput, 'SPI13', 'complete sequence', 3141722, 3165661)
    # fastaInput = '/home/nuttachat/Python/Project/WGS/maindb/SPI_FT/SPI14_gene.fasta'
    # fastaOutput = '/home/nuttachat/Python/Project/WGS/res/testresult/fasta1.fasta'
    # getSubFasta(fastaInput, fastaOutput, 1, 10)
