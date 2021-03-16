All of these are integers, and as per how Cromwell works, they only actually do anything when goleft-wdl is run from the cloud.

There's no disk size options for the report task because GCS gives you about one gigabyte if you don't specify how much storage you're using for a given task. If you end up with over a gigabyte's worth of tiny text files, something's probably gone horribly wrong.

Note that all storage related inputs do nothing on AWS, as AWS scales size automagically. Also note that preemptible related inputs only work on GCS-derived platforms, including Terra.

* covstatsAddlDisk  
	* default: 0  
	* additional storage size beyond size of inputs, in GB, to allocate to covstats  
* covstatsMem  
	* default: 8  
	* memory to request for the covstats task, in GB, when running on cloud platforms  
* covstatsPreempt  
	* default: 1  
	* number of times to attempt to run covstats on a [preemptible instance](https://cloud.google.com/compute/docs/instances/preemptible) to save money  
* indexcovAddlDisk  
  * default: 0  
  * additional storage size beyond size of inputs, in GB, to allocate to indexcov  
* indexcovMemory  
  * default: 2  
  * memory to request for the indexcov task, in GB, when running on cloud platforms  
* indexcovPrempt  
  * default: 1  
  * number of times to attempt to run covstats on a [preemptible instance](https://cloud.google.com/compute/docs/instances/preemptible) to save money  
* indexrefAddlDisk
  * default: 0
  * additional storage size beyond size of inputs, in GB, to allocate to indexing the ref genome
* indexrefMem
  * default: 2
  * memory to request for the indexing of the ref genome, in GB, when running on cloud platforms  
* indexrefPreempt
  * default: 1
  * number of times to attempt to run covstats on a [preemptible instance](https://cloud.google.com/compute/docs/instances/preemptible) to save money  
* reportMem  
  * default: 2  
  * memory to request for the report task, in GB, when running on cloud platforms  
* reportPreemptible  
  * default: 2  
  * number of times to attempt to run covstats on a [preemptible instance](https://cloud.google.com/compute/docs/instances/preemptible) to save money  
