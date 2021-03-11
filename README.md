# goleft-wdl
Collection of miscellanous WDLized [goleft](https://github.com/brentp/goleft) functions. Included are the covstats and indexcov functions.

The result of covstats is a text file, reports.txt, that prints the filename, read length, and coverage of every input file, then the average coverage and read length for the entire array of inputs. However, all covstats entries are also reported as `this.txt` in each shard's working directory, which can assist in gaining more statistics or debugging.

The result of indexcov is all of the output files associated with indexcov, which include interactive HTML files and text files. See [indexcov's github page](https://github.com/brentp/goleft/tree/master/indexcov#indexcov) for more information.

## Inputs
* `inputBamsOrCrams` is your array of input files. As the name implies it can contain CRAM files, BAM files, or a mixture of both. Please be aware that indexcov will ignore CRAM files as it currently only works on BAMs, while covstats supports both.
* If `inputBamsOrCrams` contains at least one CRAM file, you ***must*** also include a reference genome (`refGenome` in the input JSON). That reference genome ***must*** be the same one that was used to create the CRAM file.
* If `inputBamsOrCrams` contains at least BAM file, you *should* also include the respective index file for your BAM(s) (`inputIndexes` in the input JSON) but it is not a hard requirement
	* `samtools index` will be run on every BAM for which an index is not found, which can really slow things down, so if you have indicies lying around, include them
	* They should be close in name to your BAM file(s) to be found. For instance, foo.bam should have its index be named foo.bam.bai or foo.bai
	* For the purposes of accurate output from indexcov, your BAM files should be whole-genome but this is not strictly necessary if you are using this workflow as a test.
* You can optionally set [runtime attributes](https://cromwell.readthedocs.io/en/stable/RuntimeAttributes/) for the tasks in this workflow as integers. They are as follows:
		* covstatsAddlDisk
			* default: 0
			* additional storage size beyond size of inputs, in GB, to allocate to covstats
			* does nothing on AWS or local runs as they scale size automagically
		* covstatsMem
			* default: 8
			* memory to request for the covstats task, in GB, when running on cloud platforms
		Int covstatsPreempt
			* default: 1
			* number of times to attempt to run covstats on a [preemptible instance](https://cloud.google.com/compute/docs/instances/preemptible) to save money
			* does nothing unless executing on GCS (or Terra, which uses GCS)
		Int indexcovAddlDisk
			* default: 0
			* additional storage size beyond size of inputs, in GB, to allocate to indexcov
			* does nothing on AWS or local runs as they scale size automagically
		Int indexcovMemory
			* default: 2
			* memory to request for the indexcov task, in GB, when running on cloud platforms
		Int indexcovPrempt
			* default: 1
			* number of times to attempt to run covstats on a [preemptible instance](https://cloud.google.com/compute/docs/instances/preemptible) to save money
			* does nothing unless executing on GCS (or Terra, which uses GCS)
		Int reportMem
			* default: 2
			* memory to request for the report task, in GB, when running on cloud platforms
		Int reportPreemptible
			* default: 2
			* number of times to attempt to run covstats on a [preemptible instance](https://cloud.google.com/compute/docs/instances/preemptible) to save money
There's no disk size options for the report task because if you end up with over a gigabyte's worth of tiny text files, something's probably gone horribly wrong.

### Why doesn't this use the existing Biocontainer for goleft? Are there any differences in functionality?
TL;DR -- Security reasons. The main difference is that my image, due to having a different version of samtools, runs faster than the Biocontainers one at the cost of being less able to handle *incorrect* user inputs. If the user's inputs are valid -- ie, if CRAM files are being input with the same reference genome they were made with -- then handling should be equivalent.

Slightly longer explanation: Eagle-eyed users may be aware that there exists [another Docker image](https://quay.io/repository/biocontainers/goleft?tab=tags) for goleft. As it hasn't been updated in about two years and cannot be easily scanned for security purposes, I decided to make my own image. The images' differences in samtools appear to make my image run faster, but the legacy image can better handle the above cram/ref mismatch situation, likely due to the legacy image running an extra indexing step. So if you are running on dozens of crams that may have been built with different reference genomes, you may better served by [covstats-wdl](https://github.com/aofarrel/covstats-wdl/) as an alternative, as you can select which container to use there. I did not make selecting the legacy container an option in goleft-wdl, as it does not meet the standards for Biodata Catalyst Silver Workflow.

### Authorship
This repo: [Ash O'Farrell](https://github.com/aofarrel)
Original [goleft](https://github.com/brentp/goleft): [Brent Penderson](https://github.com/brentp)

 *🎶 As he goleft, and you stay right 🎶 -- The Fray, kind of*
