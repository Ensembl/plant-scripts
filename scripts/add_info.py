#!/usr/bin/env python
# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import os
import re

file=open("merged_genes.tsv","r")    
lines=file.readlines()

output=open("merged_genes.gff", "w")
list_id=[]
list_mrna=[]
main_id=""


for line in lines:
    line=line.strip()
    if line.startswith('Chr'):
        last_column=line.split("\t")[8]
        list_id=line.split("\t")[8].split(",")
        flag1=flag2=flag3=""
        for i in list_id:
            if i.startswith('ID=Os'):
                flag1=i
            elif i.startswith('ID=gene:'):
                flag2=i
            elif i.startswith('ID=LOC'):
                flag3=i
        if flag1:
            main_id=flag1

        if flag1=="":
            if flag2:
                main_id=flag2
            else:
                main_id=flag3
        main_id=main_id.replace("ID=","")

        new_last_column="ID="+str(main_id)+";Name="+str(main_id)
        new_gene=line.replace(last_column,new_last_column)
        output.write(new_gene)
        output.write("\n")


        for gene_id in list_id:
            gene_id=gene_id.replace("ID=","")
            if gene_id.startswith('Os'):
                cmd="grep -P '\smRNA\s.+Parent=" + gene_id + "' oryza_sativa_RAPDB.gff"
            if gene_id.startswith('LOC_Os'):
                cmd="grep -P '\smRNA\s.+Parent=" + gene_id + "' oryza_sativa_MSU.gff"
            if gene_id.startswith('gene:Osativa'):
                cmd="grep -P '\smRNA\s.+Parent=" + gene_id + "' oryza_sativa_gramene.gff"
            match=os.popen(cmd).read()
            list_mrna=re.split("\n",match)
            list_mrna.remove(list_mrna[-1])
            for mrna in list_mrna:
                mrna_id=mrna.split("ID=")[1].split(";Name=")[0].split(";Parent=")[0]
                parent_id=mrna.split("Parent=")[1]
                old_mrna_parent="Parent=" +str(parent_id)
                new_mrna_parent="Parent=" + str(main_id)
                new_mrna=mrna.replace(old_mrna_parent,new_mrna_parent)
                if mrna_id.startswith('Os'):
                    cmd2="grep -w Parent=" + mrna_id + " oryza_sativa_RAPDB.gff"
                if mrna_id.startswith('LOC_Os'):
                    cmd2="grep -w Parent=" + mrna_id + " oryza_sativa_MSU.gff"
                if mrna_id.startswith('transcript:Osativa'):
                    cmd2="grep -w  Parent=" + mrna_id + " oryza_sativa_gramene.gff"
                exons_match=os.popen(cmd2).read()
                output.write(new_mrna)
                output.write("\n")
                output.write(exons_match)
