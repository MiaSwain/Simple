# SIMPLE — A SIMPLE pipeline for mapping point mutations

## SYSTEM REQUIREMENTS
- Runs on mac (10.11.6) with X11 and Linux (centOS 6.7)
- Requires Java 1.7 (7u79; http://www.oracle.com/technetwork/java/javase/downloads/jdk7-downloads-1880260.html)
  - A new Simple version for Java1.8 is available. Please try it first if you have Java1.8
  - the link to it is: https://github.com/wacguy/Simple-1.8.1
- Requires R and the following packages: ggplot2, reshape2 and ggrepel
- requires Command Line Tools in Mac OS X
- Internet connection

## SIMPLE SETUP
1. Clone this repository containing the modified Simple package in your home directory
```bash
git clone https://github.com/MiaSwain/Simple.git
```

2. Create the simple conda environment
```bash
conda env create --file conda/simple.yml --name simple
```

3. Make the `simple.sh` script executable
```bash
chmod +x ./scripts/simple.sh
```

## RUNNING SIMPLE
1. Rename your fastq files as follows:
  - For the mutant and WT bulk, the names should start with mut. and wt. respectively (note the dot)
  - For single or paired-end you should then have R1 or R1 and R2, respectively and end with .fastq
    - For example, if your mutant bulk was sequence in a paired-end format and the WT as single-end you should rename the three files as follow: mut.R1.fastq, mut.R2.fastq and wt.R1.fastq.

2. If you would like to have the name of your output files be specific to the line you are mapping, open the simple_variables.sh file and change the line variable from “EMS” to your line name.
  - Letters and underscores only, please. This name will be the prefix to all of your output files (this line should look like: line=linename)

3. Place the renamed fastq files in the fastq folder located in the Simple folder.

4. Open the folder scripts inside Simple; open the data_base.txt file. Locate your species in the first column and copy it. 

5. Open the file simple_variables.sh inside the folder scripts with a text editor and paste the species name you've just copied to replace Arabidopsis_thaliana as the species name (e.g., this line should look like: my_species=Arabidopsis_thaliana or my_species=Oryza_sativa_Japonica). save the file.

6. If the mutation you are mapping is dominant, in the simple_variables.sh file, change the mutation from recessive to dominant. This line should look like: mutation=dominant. save the file.

7. Open the Terminal application.

8. Navigate to the Simple directory
```bash
cd ~/Simple
```

9. Begin a screen session so your job can run in the background
```bash
screen -S <session_name>
```

10. Activate the simple_screen conda environment
```bash
conda activate simple_screen
```

11. Activate the simple conda environment
```bash
conda activate simple
```

12. Execute the program
```bash
./scripts/simple.sh
```

13. The script will run for a few hours up to a couple of days, depending on the size of your fastq files and the size of the genome you are working with. You will know it finished once the prompt shows the following colorful text: “Simple is done”.

## Notes on using screen
- To exit a screen: Ctrl + A + D
- To return to a screen: `screen -r <screen_name>`



