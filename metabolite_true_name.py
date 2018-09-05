#script to produce a list of actual metabolite games
#author: Noah A. Legall
#created: August 23rd 2018
##########

#Take in arguements and save to variable filename
import sys # use to access arguments
import os # use in order to call commands from terminal script is called in

#read in file name you want to parse
filename = sys.argv[1] 

#name of the output file we wish to create, then activly create it
output = sys.argv[2] 
os.system('touch '+output)

#open csv with pairing info
#open csv with enrichment info
metabolite_file = open("metabolite_names.csv", "r")
enrichment_file = open(filename,"r")
#create a dictionary which is a simple map
metabolite_name_converter = dict()

#create an array that can be easier to iterate
new_file = []

for line in metabolite_file.readlines():
    new_line1 = line.replace(':','_')
    new_line2 = new_line1.replace('-','_')
    new_file.append(new_line2.rstrip())

#fill the dictionary with the pairings
for i in range(0,len(new_file)-1,2):
    metabolite_name_converter[new_file[i]] = new_file[i+1]  
#go through the enriched list of metabolites and get their true names.
#but first a little processing
enriched_metabolites = []

for line in enrichment_file:
    #could be a better way to do this, just didn't want to use regex
    new_line1 = line.replace('x_','')
    new_line2 = new_line1.replace('_,',',')
    new_line3 = new_line2.replace(':','_')
    new_line4 = new_line3.replace('-','_')
    extract = new_line4.split(",")
    if(extract[8] == "enriched\n"):
        enriched_metabolites.append(extract[0])
    else:
        continue
true_name = []

#use the metabolite_name_converter dictionary created earlier to name the enriched metabolites
for name in enriched_metabolites:
    if name in metabolite_name_converter:
        true_name.append(metabolite_name_converter[name])
    else:
        continue
    
output_file = open(output,"w")
for metabolite in true_name:
    output_file.write(metabolite+"\n")
    

