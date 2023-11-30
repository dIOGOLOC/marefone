#!/usr/bin/env python
# coding: utf-8

import os
import glob
import shutil
from tqdm import tqdm
import time
import sys
import subprocess

# ======
# Config
# ======

FOLDER_INPUT = '/medata02/HDs_01.04.22/'

FOLDER_OUTPUT = '/media/sysop/8d2362fc-3b46-49a7-a864-19b2a6ad097b/diogoloc/dados_posdoc/gliders_project/OUTPUT/'

# ========
# Function
# ========

def copyWithSubprocess(cmd):        
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

# =======
# Program
# =======

print("Searching PDF files...")
start_time = time.time()

pdf_files = []
for root, dirs, files in os.walk(FOLDER_INPUT):
    for file in files:
        if file.endswith('.pdf'):
            pdf_files.append(os.path.join(root, file))

pdf_files_S = sorted(pdf_files)
print("--- %.2f execution time (min) ---" % ((time.time() - start_time)/60))

# ================
# Copying pdf file
# ================

max_ = len(pdf_files_S)
for i in tqdm(pdf_files_S,total=max_, desc='Copying files'):
    cmd=None
    source = i
    dest = FOLDER_OUTPUT+'PDF/'+i.split('/')[-1]
    cmd=['cp', source, dest]
    copyWithSubprocess(cmd)

print(str(len(pdf_files_S))+' files copied!')
print('\n')