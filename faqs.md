# FAQs

### Why does indexing occur twice, once in indexcov and once in covstats, if an index file for a given input cannot be found?
covstats in particular expects an index file to be in the same file as the cram/bam file itself. This is never the case when an input is passed in from another task in Cromwell, which is what happens if I index in another task and then try to pass that input in. There is a possible workaround involving symlinks, but it only works on local runs, not the cloud, so I haven't included it.  

### Why are there two versions of indexcov?
TL;DR -- It is much less complicated for both me and users than making one version of indexcov. If I made it one task, I'd need to make this README even longer with a detailed explanation of the very precise needs of indexcov based on whether you are using bams, crams, or a mixture of both.  

Longer explanation: Cromwell lacks a concept of mutual exclusivity or `else`, so something like this:
```
if(a):
	call foo {input: bar = bar.bai}
if(!a):
	call foo {input bar = bar.crai}
```
Throws an error. This means you can only have one task called foo, and foo's inputs must be valid regardless of the truth of a. So that means foo must somehow account for the a and the !a situation, while still following all of the rules about how input files must be localized prior to a task beginning, how Cromwell handles optional inputs, etc. Doing this is technically possible, but has extremely specific input requirements that just aren't friendly to the user.

### Why doesn't this use the existing Biocontainer for goleft? Are there any differences in functionality?
TL;DR -- Security reasons. The main difference is that my image, due to having a different version of samtools, runs faster than the Biocontainers one at the cost of being less able to handle *incorrect* user inputs. If the user's inputs are valid -- ie, if CRAM files are being input with the same reference genome they were made with -- then handling should be equivalent.

Longer explanation: Eagle-eyed users may be aware that there exists [another Docker image](https://quay.io/repository/biocontainers/goleft?tab=tags) for goleft. As it hasn't been updated in about two years and cannot be easily scanned for security purposes, I decided to make my own image. The images' differences in samtools appear to make my image run faster, but the legacy image can better handle the above cram/ref mismatch situation, likely due to the legacy image running an extra indexing step. So if you are running on dozens of crams that may have been built with different reference genomes, you may better served by [covstats-wdl](https://github.com/aofarrel/covstats-wdl/) as an alternative, as you can select which container to use there. I did not make selecting the legacy container an option in goleft-wdl, as it does not meet the standards for Biodata Catalyst Silver Workflow.