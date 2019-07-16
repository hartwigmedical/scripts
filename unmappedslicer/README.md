# Unmapped Slicer

The unmapped slicer can be used to slice all unmapped reads from a BAM file that lives remotely and is accessible through a https url. This tool has been developed by Roel Janssen, UMC Utrecht.  

The steps to run are as follows:
 - Compile unmapped_slicer.c on the machine you run the tool on
 
 ```
 gcc -g `pkg-config htslib libcurl --libs --cflags` ./unmapped_slicer.c -o unmapped_slicer
 
```
 - Run the tool on the BAM
 
 ```
./unmapped_slicer ${bam_url} ${bai_url} ${raw_unmapped_output_bam}
```

The tool includes some mapped reads at the start of the file (in addition to all unmapped reads). These need to be filtered out as a post-process step to retain only unmapped reads.

```
sambamba view -f bam -F 'unmapped' ${raw_ummapped_output_bam} > ${final_unmapped_output_bam}
```