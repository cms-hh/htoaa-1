import os
import sys
import subprocess
import json
import uproot3 as uproot

from htoaa_Settings import * 

def calculate_lumiScale(luminosity, crossSection, sumEvents):
    lumiScale = 1
    # as crosssection is in pb and luminosity in fb
    pb_to_fb_conversionFactor = 1000
    if sumEvents != 0: lumiScale = luminosity * crossSection * pb_to_fb_conversionFactor / sumEvents
    return lumiScale


def setXRootDRedirector(fileName):
    if not fileName.startswith("/store/"):
        return fileName

    redirector_toUse = None
    for redirector in xrootd_redirectorNames:
        print(f"setXRootDRedirector():: Checking {redirector + fileName}"); sys.stdout.flush()
        try:
            file1 = uproot.open(redirector + fileName)
        except:
            print(f"setXRootDRedirector():: File open {redirector + fileName} failed"); sys.stdout.flush()
        else:
            nEntries = file1['Events'].numentries
            file1.close()
            if nEntries > 0:
                print(f"{redirector + fileName}: {nEntries}"); sys.stdout.flush()
                redirector_toUse = redirector
                break
    return redirector_toUse + fileName

def xrdcpFile(sFileName, sFileNameLocal, nTry = 3):
    command_ = "time xrdcp %s %s" % (sFileName, sFileNameLocal)
    command_list_ = command_.split(" ")
    print(f"{command_ = }")
    for iTry in range(nTry):
        process = subprocess.Popen(command_list_,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True
        )
        stdout, stderr = process.communicate()
        if 'FATAL' not in stderr: # download was successful
            return True
    return False

def get_lf(files, name):
    sInputFiles_toUse = []
    for f in files:
        if "*" in f:
            sInputFiles_toUse.extend(glob.glob(f))
        elif 'eos' not in f:
            sInputFile = setXRootDRedirector(f)
            sFileLocal = f'/tmp/snandan/inputFiles/{name}/{os.path.basename(f)}'
            if xrdcpFile(sInputFile, sFileLocal, nTry = 3):
                sInputFiles_toUse.append(sFileLocal)
            else:
                print(f"Ip file {f} failed to download \t **** ERROR ****")
                exit(1)
        else:
            sInputFiles_toUse.append(f)
    return sInputFiles_toUse

def GetDictFromJsonFile(filePath):
    # Lines starting with '#' are not read out, and also content between '/* .... */' are not read.
    # Content between " '''   ....  ''' " are not read
    # Source: https://stackoverflow.com/questions/29959191/how-to-parse-json-file-with-c-style-comments
    
    contents = ""
    fh = open(filePath)
    for line in fh:
        cleanedLine = line.split("#", 1)[0]
        if len(cleanedLine) > 0 and line.endswith("\n") and "\n" not in cleanedLine:
            cleanedLine += "\n"
        contents += cleanedLine
    fh.close
    
    while "'''" in contents:
        preComment, postComment = contents.split("'''", 1)
        contents = preComment + postComment.split("'''", 1)[1]

    dictionary =  json.loads( contents )
    return dictionary


def DfColLabel_convert_bytes_to_string(df):
    cols_rename = {}
    for col in df.columns:
        if isinstance(col, (bytes, bytearray)):
            cols_rename[col] = col.decode()
    print("DfColLabel_convert_bytes_to_string:: cols_rename: {}".format(cols_rename))
    df.rename(columns=cols_rename, inplace=True)
    return df
