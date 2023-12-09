# Use the official Python image for Windows Server Core 2019
FROM python:3.9-windowsservercore

# Set the working directory to /app
WORKDIR /app

# Copy the requirements file into the container at /app
COPY /* /app/

# Install any needed packages specified in requirements.txt
RUN pip install --upgrade pip --user
RUN pip install --no-cache-dir -r req.txt --user

# Install R in the Docker image using PowerShell
SHELL ["powershell", "-Command", "$ErrorActionPreference = 'Stop'; $ProgressPreference = 'SilentlyContinue';"]

# Download and install R
RUN Invoke-WebRequest -Uri https://cran.r-project.org/bin/windows/base/R-4.3.2-win.exe -OutFile R.exe
RUN Start-Process -Wait -FilePath .\R.exe -ArgumentList '/VERYSILENT', '/SUPPRESSMSGBOXES', '/DIR=C:\R'

# Set R_HOME environment variable
ENV R_HOME C:\\R\\bin\\x64

# Clean up
RUN Remove-Item -Force R.exe

CMD ["python", "app.py"]