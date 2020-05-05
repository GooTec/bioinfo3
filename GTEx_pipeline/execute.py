import subprocess
import os
input_dir = 'input_jsons/'
filelist = os.listdir('input_jsons')

for input_file in filelist :
    print(input_file)
    process = os.system('sudo java -jar cromwell-49.jar run rna-seq_pipeline.wdl -i '+input_dir+input_file)
