# -*- coding: utf-8 -*-
import numpy as np
import fnmatch
import os
import sys
from pathlib import Path


import ecephys_spike_sorting.modules.kilosort_helper.SGLXMetaToCoords as SGLXMeta

def GetFirstTrialPath(catGT_run_name, gate_string, trigger_string, probe_string ):
    prb_list = ParseProbeStr(probe_string)
    run_folder = catGT_run_name + '_g' + gate_string
    prb_folder =  run_folder + '_imec' + prb_list[0]
    first_trig, last_trig = ParseTrigStr(trigger_string, prb_folder)
    filename = catGT_run_name + '_g' + gate_string + '_t' + \
                 str(first_trig) + '.imec' + prb_list[0] + '.ap.bin'
    firstTrialPath = os.path.join(run_folder, prb_folder, filename) 

    return firstTrialPath

def GetTrialRange(prb_folder):
    tFiles = os.listdir(prb_folder)
    minIndex =  sys.maxsize
    maxIndex = 0
    for tName in tFiles:
        if (fnmatch.fnmatch(tName,'*.ap.bin')):
            parts = tName.split('_')
            tparts = parts[2].split('.')
            tInd = int(tparts[0][1:])
            if tInd > maxIndex:
                maxIndex = tInd
            if tInd < minIndex:
                minIndex = tInd
    return minIndex, maxIndex

    
def EphysParams(ap_band_file):
    # assume metadata file is in same directory as binary, Constuct metadata path
    
    # read metadata
    
    metaName, binExt = os.path.splitext(ap_band_file)
    metaFullPath = Path(metaName + '.meta')
    meta = SGLXMeta.readMeta(metaFullPath)
    
    if 'imDatPrb_type' in meta:
        pType = (meta['imDatPrb_type'])
        if pType =='0':
            probe_type = 'NP1'
        else:
            probe_type = 'NP' + pType
    else:
        probe_type = '3A'    #3A probe
    
    sample_rate = float(meta['imSampRate'])    
    
    num_channels = int(meta['nSavedChans'])
    
    uVPerBit = Chan0_uVPerBit(meta)
      
    return(probe_type, sample_rate, num_channels, uVPerBit)

# Return gain for imec channels.
# Index into these with the original (acquired) channel IDs.
#
def Chan0_uVPerBit(meta):
    # Returns uVPerBit conversion factor for channel 0
    # If all channels have the same gain (usually set that way for 
    # 3A and NP1 probes; always true for NP2 probes), can use
    # this value for all channels.
    
    imroList = meta['imroTbl'].split(sep=')')
    # One entry for each channel plus header entry,
    # plus a final empty entry following the last ')'
    # channel zero is the 2nd element in the list

    if 'imDatPrb_dock' in meta:
        # NP 2.0; APGain = 80 for all channels
        # voltage range = 1V
        # 14 bit ADC
        uVPerBit = (1e6)*(1.0/80)/pow(2,14)
    else:
        # 3A, 3B1, 3B2 (NP 1.0)
        # voltage range = 1.2V
        # 10 bit ADC
        currList = imroList[1].split(sep=' ')   # 2nd element in list, skipping header
        APgain = float(currList[3])
        uVPerBit = (1e6)*(1.2/APgain)/pow(2,10)
        
    return(uVPerBit)

def ParseProbeStr(probe_string):
    
    str_list = probe_string.split(',')
    prb_list = []
    for substr in str_list:
        if (substr.find(':') > 0):
            # split at colon
            subsplit = substr.split(':')
            for i in range( int(subsplit[0]), int(subsplit[1]) + 1):
                prb_list.append(str(i))
        else:
            # just append this string
            prb_list.append(substr)

    return prb_list

def ParseTrigStr(trigger_string, prb_folder):
    
    str_list = trigger_string.split(',')
    first_trig_str = str_list[0]
    last_trig_str = str_list[1]
    
    if last_trig_str.find('end') >= 0 or first_trig_str.find('start') >= 0 :
        # get the full range from the directory
        minInd, maxInd = GetTrialRange(prb_folder)

    if first_trig_str.find('start') >= 0:
        first_trig = minInd
    else:
        first_trig = int(first_trig_str)
    
    if last_trig_str.find('end') >= 0:
        last_trig = maxInd
    else:
        last_trig = int(last_trig_str)
        
    # trig_array =  np.arange(first_trig, last_trig+1)

    return first_trig, last_trig


def ParseTcatName(tcat_name):
    
    parts_list = tcat_name.split('.')
    baseName = parts_list[0] + '_' + parts_list[1]
    return baseName

def GetProbeStr(tcat_name):
    parts_list = tcat_name.split('.')
    imStr = parts_list[1]
    if len(imStr) == 4:
        prbStr = ''      # 3A data, no probe index
    else:
        prbStr = imStr[4:len(imStr)]
    return prbStr


def ParseCatGTLog(logPath, run_name, gate_string, prb_list):

    gfix_str = run_name + '_g' + gate_string + ' Gfix'

    num_probe = len(prb_list)
    gfix_edits = np.zeros(num_probe, dtype='float64')

    gfound = np.zeros(num_probe)
    pfound = list()             # list of strings of probes found
    nfound = 0
    log_fullpath = logPath.replace('\\', '/') + "/CatGT.log"

    with open(log_fullpath, 'r') as reader:
        line = reader.readline()
        while line != '' and nfound < num_probe:  # The EOF char is an empty string
            gstart = line.find( gfix_str )
            if gstart  > -1:      
                # parse this line to get probe string and corrections/sec
                line_end = len(line)
                gsub = line[gstart:line_end]
                datArr = gsub.split()       
                pfound.append(datArr[3])
                gfound[nfound] = float(datArr[5])
                nfound = nfound + 1
            line = reader.readline()   
    
    # order the returned gfix_edits matching the probe order specified 
    # in prb_list
    for i in range(0,len(prb_list)):
        gfix_edits[i] = gfound[pfound.index(prb_list[i])]
     
    return gfix_edits


