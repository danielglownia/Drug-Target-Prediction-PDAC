import subprocess
import os

def run_runCIRCUST_R(tsv_file_path):
  """
  Runs the "runCIRCUST.R" command line with the provided TSV file path.

  Args:
    tsv_file_path: The path to the TSV file.

  Returns:
    A tuple containing the standard output and standard error of the command as strings.

  Raises:
    FileNotFoundError: If the provided TSV file path is not found.
  """

  # Check if the provided TSV file exists
  if not os.path.exists(tsv_file_path):
    raise FileNotFoundError(f"File not found: {tsv_file_path}")

  # Build the command line
  command = ["Rscript", "./CIRCUST-main/runCIRCUST.R", "--tsv", tsv_file_path]

  # Run the command line and capture the output and error
  process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
  output, error = process.communicate()

  # Decode the output and error from byte strings to text strings
  output = output.decode("utf-8")
  error = error.decode("utf-8")

  return output, error
