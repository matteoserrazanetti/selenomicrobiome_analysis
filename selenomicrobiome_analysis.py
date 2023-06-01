#! /usr/bin/env python3

import sys, os, shlex, subprocess
from easyterm import *
import shutil
import multiprocessing

help_msg="""This pipeline performs different kind of analysis on a dataset of metagenomes assemblies to characterize the 
metagenomes (markers for Se usage and selenoprotein families)

gives aboslute paths in all folder requirements

example of command line:
python3 selenomicrobiome_analysis.py -i input_directory -o output_directory -db database_file 

Options
-db database folder
-i input_directory (fna files)
-o output_directory (if not present will be created)
-S False to not perform selenoprofiles step
-T False to not perform secmarker step
-I False to not perform infernal step
-e e-value for BLAST (default 1e-5)
-n cpu cores (default 8)                                                                                                  
-p profiles that selenoprofiles4 need to search (default prokarya)
-gff output from selenoprofiles in gff form (default True)
-gtf output from selenoprpfiles in gtf form (default True)
----------------------------------------
-cm cm file for RnaseP search with infernal (depending on how to distributw it but i can make it default)
-sp4_cfg .selenoprofile_config.txt file
-----------------------------------------
"""

def_opt= {'i':'', 'o':'','db':'undefined_DB',
          'S':True, 'T':True, 'I':True,
          'e':1e-5,  'n':8, 'p':'prokarya','gtf':True, 'gff':True,
          'cm':'/users/rg/mserrazanetti/wgs_data/cm/all_RF00.c.cm', 'sp4_cfg':'/users/rg/mserrazanetti/.selenoprofiles_config.txt'}

#########################   REMEMBER TO CHANGE PATH OF CM AND SP4_CONFIG ##########################

########################  catching up from MMlib #########################
  
def bbash(cmnd, print_it=0,  out=None, err=None, tolerate=False, shell=False):
  if out is None and err is None:
    out=subprocess.PIPE
    err=subprocess.STDOUT
  elif out is None or err is None:
    raise Exception('bbash ERROR if you specify out then specify err too')
  p=subprocess.run( shlex.split(cmnd), stdout=out, stderr=err, check=False, shell=shell)
  o=p.stdout
  if not tolerate and p.returncode!=0:
    if o:
      raise Exception("bbash ERROR command: {} error: {}".format(cmnd, o.decode(encoding = 'utf-8')))
    else:
      raise Exception("bbash ERROR command: {} error --> check {}".format(cmnd, o.__file__))
  if o:
    return o.decode(encoding = 'utf-8')
  else:
    return o

###########################################################################

def create_folder(path):
  if not os.path.exists(path):
    os.mkdir(path)
  if os.path.exists(path): # in this part i remove and re-create the folder if already present, this way i remove possible confusions
    shutil.rmtree(path) # remove
  os.mkdir(path) # re-create
   
#########################    Multicores main function    ####################################

def main_function(arg_list):
  
  filename = arg_list[0] 
  opt = arg_list[1]
  marker_list = arg_list[2]
  fasta_ext = arg_list[3]
  id=filename.split('/')[-1].split('.')[0]

  output_folder=opt['o']
  if output_folder[-1]=='/':
    output_folder=output_folder[:-1]
  if not os.path.exists(output_folder):
    os.mkdir(output_folder)

  if opt['S']:
    print('Starting selenoprofiles4 for '+id)
    sp_fam_list=create_fam_list(opt['sp4_cfg'], marker_list)
    selenoprofiles(filename, output_folder, opt['p'], id)
    count_selenoprotein_families(output_folder, id, marker_list, sp_fam_list)
  if opt['T']:
    print('Starting secmarker for '+id)
    with open(output_folder+'/all_secmarker.txt','w') as outfile:
      outfile.write('ID trnasec')
    secmarker(filename, output_folder, id)
    secmarker_count(filename, output_folder, id)
  if opt['I']:
    print('Starting infernal for '+id)
    infernal(filename, output_folder, opt['cm'], id)
    clean_infernal_output(output_folder, id)


#########################    Selenoprofiles part  ####################################	

# important that file option needs to take a absoulte path for the file

def create_fam_list(sp4_cfg, marker_list):
  sp_fam_list = []
  with open(sp4_cfg,'r') as infile:
    for line in infile:
      if line[0]=='f':
        line=line.split()
        if line[0]=='families_set.prokarya':
          l = line[2].split(',') 
          for el in l:
            if el not in marker_list:
              sp_fam_list.append(el)
  return(sp_fam_list)

def selenoprofiles(file, output_folder, profile, id):
  #bbash('conda activate sp4')
  sp4_folder=output_folder+'/sp4'
  if not os.path.exists(sp4_folder):
    os.mkdir(sp4_folder)
  if os.path.exists(sp4_folder+'/'+id): # in this part i remove and re-create the folder if already present, this way i remove possible confusions
    shutil.rmtree(sp4_folder+'/'+id) # remove     
  os.mkdir(sp4_folder+'/'+id) # re-create
  new_output_folder = sp4_folder+'/'+id+'/'+id # create a new folder inside the id_folder, this is done because the databases created by sp4 need to be separated
  gff_output = ''
  gtf_output = ''
  if def_opt['gff']:
    gff_output = ' -output_gff_file '+new_output_folder+'/'+'sp4_results_'+id+'.gff '
  if def_opt['gtf']:
    gtf_output = ' -output_gtf_file '+new_output_folder+'/'+'sp4_results_'+id+'.gtf '
  print('selenoprofiles -o '+new_output_folder+' -t '+file+' -s '+id+' -p '+profile+' -no_splice '+gff_output+gtf_output)
  start_folder=os.getcwd()
  os.chdir(sp4_folder+'/'+id)
  bbash('selenoprofiles -o '+new_output_folder+' -t '+file+' -s '+id+' -p '+profile+' -no_splice '+gff_output+gtf_output)
  os.chdir(start_folder)
  #bbash('conda deactivate')

def count_selenoprotein_families(output_folder, id, marker_list, sp_fam_list):
  for dir, subdirs, files in os.walk(output_folder+'/sp4'):
    for file in files:
      if file[:11]=='sp4_results' and file[-3:]=='gff':
        with open(dir+'/'+file, 'r') as infile:
          count=0
          sp_fam_dict_partial = {}
          sp_fam_dict_partial['seld_sec']=0
          for el in sp_fam_list: #initialze the partial dictionary with selenoproteins genes
            sp_fam_dict_partial[el]=0
          for el in marker_list: #initialize the partial dictionary with marker genes (actually genes that do not use selenocysteine)
            sp_fam_dict_partial[el]=0
          for line in infile:
            line = line.split()
            sp_fam=line[8].split('.')[0]
            if sp_fam[:4]=='Sec1' or sp_fam[:4]=='Sec2' or sp_fam[:4]=='Sec3': #correct the name of the family when is predicted with SeC
              sp_fam = sp_fam.split(':')[1]
              label=line[8].split('.')[-1]
            if sp_fam == 'seld' and label == 'selenocysteine':
              sp_fam_dict_partial['seld_sec']+=1
            if (sp_fam in sp_fam_list and label == 'selenocysteine') or (sp_fam in marker_list and sp_fam != 'seld' and (label == 'homologue' or label == 'pseudo')) or (sp_fam == 'seld'):
              print(sp_fam+'\t'+str(line[8].split('.')))
              count+=1
            if sp_fam not in sp_fam_dict_partial.keys():
              sp_fam_dict_partial[sp_fam]=1
            else:
              sp_fam_dict_partial[sp_fam]+=1

        with open(dir+'/sp4_family_count_'+id+'.txt','w') as outfile:
          outfile.write('  '+id+'\n')
          outfile.write('total '+str(count)+'\n')
          for x,y in sp_fam_dict_partial.items():
            outfile.write(str(x)+' '+str(y)+'\n')
  
         

#######################         Secmarker part     ####################################

def secmarker(file, output_folder, id):
  #bbash('conda activate secmarker')
  #if not os.path.exists(output_folder+'/secmarker'):
  #  os.mkdir(output_folder+'/secmarker')
  id=file.split('/')[-1].split('.')[0]
  if os.path.exists(output_folder+'/secmarker/'+id): # in this part i remove and re-create the folder if already present, this way i remove possible confusions
    shutil.rmtree(output_folder+'/secmarker/'+id) # remove     
  os.mkdir(output_folder+'/secmarker/'+id) # re-create
  print('secmarker -t '+file+' -o '+output_folder+'/secmarker/'+id+' -cpu 1')
  bbash('secmarker -t '+file+' -o '+output_folder+'/secmarker/'+id+' -cpu 1')
  #bbash('conda deactivate')

def secmarker_count(file, output_folder, id):
  with open(output_folder+'/secmarker/'+id+'/trnasec.gff','r') as infile:
    count=0
    for line in infile:
      count+=1
  with open(output_folder+'/all_secmarker.txt','w') as outfile:
    outfile.write(id+' '+str(count)+'\n')

########################       Infernal for rnaseP part    ###########################

def infernal(file, output_folder, cm_file, id):
  #bbash('conda activate infernal')
  #if not os.path.exists(output_folder+'/infernal'):
  #  os.mkdir(output_folder+'/infernal')
  if os.path.exists(output_folder+'/infernal/'+id): # in this part i remove and re-create the folder if already present, this way i remove possible confusions
    shutil.rmtree(output_folder+'/infernal/'+id) # remove     
  os.mkdir(output_folder+'/infernal/'+id) # re-create
  print(cm_file)
  print(id)
  print(output_folder)
  print(file)
  print('cmscan --rfam -E 1e-5 --cpu 1 --tblout '+output_folder+'/infernal/'+id+'/infernal_output_'+id+'.tbl '+cm_file+' '+file)
  bbash('cmscan --rfam -E 1e-5 --cpu 1 --tblout '+output_folder+'/infernal/'+id+'/infernal_output_'+id+'.tbl '+cm_file+' '+file)
  #bbash('conda deactivate')

def clean_infernal_output(output_folder, id):
  open(output_folder+'/infernal/'+id+'/infernal_output_'+id+'_cleaned.tbl', 'w').close() #empty the file before writing again on it
  with open(output_folder+'/infernal/'+id+'/infernal_output_'+id+'.tbl', 'r') as infile:
    positive_dict={}
    positive_single_dict={}
    for line in infile:
      if line[0]!='#':
        line=line.split()
        if line[2] not in positive_single_dict:
          positive_single_dict[line[2]]=(line[0],line[14])
        else:
          if float(positive_single_dict[line[2]][1]) < float(line[14]):
            positive_single_dict[line[2]]=(line[0],line[14])
      if '[ok]' in line:
        with open(output_folder+'/infernal/'+id+'/infernal_output_'+id+'.tbl', 'a') as outfile:
          outfile.write(str(id)+' '+str(len(positive_single_dict))+'\n')
          outfile.write('Contig/Ch/Other Class Score\n')
          for key, value in positive_single_dict.items():
            outfile.write(str(key)+' '+str(value[0])+' '+str(value[1])+'\n')
        with open(output_folder+'/all_rnaseP.txt','w') as outfile_all:
          outfile_all.write(str(id)+' '+str(len(positive_single_dict))+'\n')

###################################      Main     ###################################

def main(args={}):

  fasta_ext=['fna','fa','faa','frn','ffn','fasta']
  marker_list = ['tigr04348', 'egtB', 'seld', 'yqeb', 'yqec', 'ybbb']

  print('Starting all the jobs')

  #arguments reading
  if not args:
    opt=command_line_options(def_opt, help_msg)
  else:
    opt=args

  input_folder=opt['i']
  if input_folder[-1]=='/':
    input_folder=input_folder[:-1]

  output_folder=opt['o']
  if output_folder[-1]=='/':
    output_folder=output_folder[:-1]
  if not os.path.exists(output_folder):
    os.mkdir(output_folder)

  if multiprocessing.cpu_count() < int(opt['n']):
    cpu=multiprocessing.cpu_count()
  else:
    cpu=int(opt['n'])

  sp4_folder=output_folder+'/sp4'
  if not os.path.exists(sp4_folder):
    os.mkdir(sp4_folder)  
  if not os.path.exists(output_folder+'/secmarker'):
    os.mkdir(output_folder+'/secmarker')
  if not os.path.exists(output_folder+'/infernal'):
    os.mkdir(output_folder+'/infernal')


  #Orrible solution to create a list of list with many things repeated but it works
  arg_list=[]
  for dir, subdirs, files in os.walk(input_folder):
    for file in files:
      if file.split('.')[-1] in fasta_ext:
        arg_list.append([dir+'/'+file,opt,marker_list,fasta_ext])


  pool = multiprocessing.Pool(processes=cpu)
  pool.map(main_function, arg_list)
  pool.close()

if __name__ == "__main__":
  main()



