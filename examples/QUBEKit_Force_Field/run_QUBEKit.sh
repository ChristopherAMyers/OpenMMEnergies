#!/bin/bash

pdb=moh_acceptor.pdb
xml=QUBEKit.xml
xyz=frames_moh_acceptor.xyz

python $code -pdb $pdb -xyz $xyz -xml $xml -ipt input
