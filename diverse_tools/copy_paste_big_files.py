#!/usr/bin/env python
# coding: utf-8

import os


# =======
# Program
# =======

file_list_log = sorted([root+'/'+file for root, dirs, files in os.walk('/mnt/medata02/HDs_01.04.22/') for file in files if file.endswith('.log')])

file_list_nc = sorted([root+'/'+file for root, dirs, files in os.walk('/mnt/medata02/HDs_01.04.22/') for file in files if file.endswith('.nc')])

file_list_wav = sorted([root+'/'+file for root, dirs, files in os.walk('/mnt/medata02/HDs_01.04.22/') for file in files if file.endswith('.wav')])

print('Number of NC files:',len(file_list_log))

hds = []
for i in file_list_log:
    if not 'RECYCLE.BIN' in i:
        hds.append(i.split('/')[4])

print(sorted(list(set(hds))))
