import subprocess
import os

def run_runCIRCUST_R(tsv_file_path):

  # Build the command line
  path = os.getcwd()
  print("Current path: ", path)
  dir_list = os.listdir(path)
 
  print("Files and directories in '", path, "' :")
 
  # prints all files
  print(dir_list)
  cir_path = "C:\\app\\CIRCUST-main\\runCIRCUST.R"
  print(cir_path)
  command = ["Rscript", cir_path]

  # Run the R script and capture the output
  output = subprocess.check_output(command)

  # Print the output
  print(output)

  return output
