Walts test 2
# goleft-wdl
Collection of miscellanous WDLized [goleft](https://github.com/brentp/goleft) functions. Included are the covstats and indexcov functions.

The result of covstats is a text file, reports.txt, that prints the filename, read length, and coverage of every input file, then the average coverage and read length for the entire array of inputs. However, all covstats entries are also reported as this.txt in each shard's working directory, which can assist in gaining more statistics or debugging.

The result of indexcov is all of the output files associated with indexcov, which include interactive HTML files and text files. See [indexcov's github page](https://github.com/brentp/goleft/tree/master/indexcov#indexcov) for more information.

## Inputs
### Input array
`inputBamsOrCrams` is your array of input files. As the name implies it can contain CRAM files, BAM files, or a mixture of both.

If `inputBamsOrCrams` contains at least one CRAM file:
* You ***must*** also include `refGenome` 
* That reference genome ***must*** be the same one that was used to create the CRAM file(s)
* Mixing CRAM files made from different ref genomes is not supported 

If `inputBamsOrCrams` contains at least one BAM file:
* You should also include the respective index file for your BAM(s) (`inputIndexes` in the input JSON) but it is not a hard requirement  
* `samtools index` will be run on every BAM for which an index is not found, which can really slow things down, so if you have indicies lying around, include them  
* foo.bam's index will only be found if it is called foo.bam.bai or foo.bai
* To get accurate output from indexcov, your BAM files should be whole-genome, but this is not strictly necessary  

### Other inputs
`forceIndexcov` is a boolean which indicates whether to run indexcov even if there is not an index file for every input (determined by comparing size of the arrays). **Please read the efficiency tips section before setting this.**

Additionally, you can optionally set [runtime attributes](https://cromwell.readthedocs.io/en/stable/RuntimeAttributes/) for the tasks in this workflow. Due to how Cromwell works, they only actually do anything in the cloud. [See this document for a full list of runtime attributes](https://github.com/aofarrel/goleft-wdl/blob/main/README_runtime_attributes.md). 

## Efficiency tips
It is *highly* recommended that you include indexes for all of your inputs files in order to skip indexing with samtools.  

If you do not include indexes for all of your inputs, in the best case scenario (all cram files, skipping indexcov) there is no change in execution time. In the worst case scenario (all full-size bam files, forceIndexcov=True) `samtools index` has to be run for every input twice.  

### Authorship
This repo: [Ash O'Farrell](https://github.com/aofarrel)  
Original [goleft](https://github.com/brentp/goleft): [Brent Penderson](https://github.com/brentp)  



 *🎶 As he goleft, and you stay right 🎶 -- The Fray, kind of*
