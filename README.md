![](https://img.shields.io/badge/build-passing-green.svg)
![](https://img.shields.io/badge/version-1.0.0-blue.svg)
![](https://img.shields.io/badge/htsjdk-2.24.1-blue.svg)
![](https://img.shields.io/badge/java-1.8-red.svg)

# FastReadCounter
Tool for fast counting of reads from BAM file into different formatted genomic ranges (GTF, BED, VCF)

## Download software
FastReadCounter (FRC) command-line tools are provided as a [single executable jar file](../master/releases/FastReadCounter-1.0.jar?raw=true).
The .jar file contains all required materials and can be run on any terminal.

## Dependencies
### Java version
For the tools to run properly, you must have Java >=1.8 installed. 

To check your java version, open your terminal application and run the following command:

```bash
java -version
```

If the output looks something like java version "1.8.x" or above, you are good to go. 
If not, you may need to update your version; see the [Oracle Java website](http://www.oracle.com/technetwork/java/javase/downloads/) to download the latest JRE (for users) or JDK (for developers).

### htsjdk Java library
The software relies on the [htsjdk library](https://github.com/samtools/htsjdk) for reading SAM/BAM files, but the JAR is already embedded in the released JAR, so no need to install it yourself.

## How to run
To check that FastReadCounter is working properly, run the following command:

```bash
java -jar FastReadCounter.jar
```
As shown in the output of this command, FastReadCounter has the following options

```
Options:
-bam %s         [Required] Path of BAM file [do not need to be sorted or indexed].
-gtf %s         [Required (or -bed or -vcf)]Path of GTF file.
-bed %s         [Required (or -gtf or -vcf)]Path of BED file.
-vcf %s         [Required (or -gtf or -bed)]Path of VCF file.
-s %s           Do you want to count only reads falling on same strand than feature? [no, yes, reverse] [default = no].
-q %i           Minimum quality required for a read to be counted [default = 10].
-o %s           Output folder [default = folder of BAM file]
```

## Comparison with existing tools
FastReadCounter can be compared with a variety of tools that all perform the same task on only one file format (bed or gtf, or vcf) such as:
* **htseq-count** for counting reads in GTF file. FastReadCounter has stricly identical output than the union parameter in htseq-count. It runs ~20-100 times faster.
* **featureCounts** for counting reads in GTF file. FastReadCounter has identical output and runs only slightly slower. However FastReadCounter can read more file formats (BED, VCF) and can read multiplexed BAM files (single-cell RNA-seq or multiplexed bulk RNA-seq experiments).
* **bedtools** for counting reads in BED files (multicov tool). FastReadCounter has very similar output than HOMER. It does not require to transform the bam into a tag folder.
* **homer** for counting reads in BED files (annotatePeaks.pl tool). FastReadCounter has very similar output than HOMER. It does not require to transform the bam into a tag folder.
* **freebayes** for counting reads in VCF files with the --variant-input option. However freebayes current version fails with this option, and only calculate chr1.

## Directory content
* **releases**: jar releases of FastReadCounter
* **src** and **lib** folders will come soon

## Author
Vincent Gardeux - vincent.gardeux@epfl.ch
