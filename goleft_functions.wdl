version 1.0

task indexRefGenome {
	input {
		# Not actually optional as only called if refGenome is defined
		File? refGenome

		# runtime attributes with defaults
		Int indexrefMem = 2
		Int indexrefPreempt = 1
		Int indexrefAddlDisk = 0
	}
	command <<<
		# samtools faidx tilde-curlyL-refGenome-curlyR somehow puts the fai in inputs folder
		# instead of the execution folder, which results in Cromwell's failure to find output.
		# Also, -o "whatever.fai" argument in samtools is ignored. Furthermore, basename
		# does not work on File? but we can't make refGenome required or else the overall
		# workflow will require refGenome, which we don't want.
		#
		# Hence, this disgusting workaround.
		# This basically doubles the size of disk space needed... ew.
		cp ~{refGenome} .
		mv "~{basename(select_first([refGenome, 'dummy']))}" "ref_copy.fa"
		samtools faidx "ref_copy.fa"
	>>>
	output {
		# select_first needed as basename does not work on File? types
		File refIndex = "ref_copy.fa.fai"
	}

	# Estimate disk size required
	Int refSize = ceil(size(refGenome, "GB"))
	Int finalDiskSize = 2*refSize + indexrefAddlDisk
	runtime {
		docker: "quay.io/aofarrel/goleft-covstats:circleci-push"
		preemptible: indexrefPreempt
		disks: "local-disk " + finalDiskSize + " HDD"
		memory: indexrefMem + "G"
	}
}



task indexcovCRAM {
	input {
		File inputCram
		Array[File] allInputIndexes
		# Not actually optional
		File? refGenomeIndex

		# runtime attributes
		Int indexcovMemory = 2
		Int indexcovPrempt = 1
		Int indexcovAddlDisk = 0
	}

	command <<<

		set -eux -o pipefail

		# Double-check this is actually a CRAM file
		AMIACRAM=$(echo ~{inputCram} | sed 's/\.[^.]*$//')
		if [ -f ${AMIACRAM}.bam ]; then
			>&2 echo "Somehow a bam file got into the cram function!"
			>&2 echo "This shouldn't happen, please report to the dev."
			exit 1
		
		else
			if [ -f ~{inputCram}.crai ]; then
				echo "Bai file already exists with pattern *.cram.crai"
			elif [ -f ${AMIACRAM}.crai ]; then
				echo "Bai file already exists with pattern *.crai"
				mv ~{inputCram}.bai ${AMIACRAM}.cram.crai
			else
				echo "Input crai file not found. We searched for:"
				echo "--------------------"
				echo "  ~{inputCram}.crai"
				echo "--------------------"
				echo "  ${AMIACRAM}.crai"
				echo "--------------------"
				echo "Finding neither, we will index with samtools."
				samtools index ~{inputCram} ~{inputCram}.crai
			fi

			INPUTCRAI=$(echo ~{inputCram}.crai)
			mkdir indexDir
			ln -s ~{inputCram} indexDir~{basename(inputCram)}
			ln -s ${INPUTCRAI} indexDir~{basename(inputCram)}.crai
			goleft indexcov --extranormalize -d indexDir/ --fai ~{refGenomeIndex} ~{inputCram}.crai
		fi

	>>>
	# Estimate disk size required
	Int indexSize = ceil(size(allInputIndexes, "GB"))
	Int thisAmSize = ceil(size(inputCram, "GB"))
	Int finalDiskSize = indexSize + thisAmSize + indexcovAddlDisk
	output {
		# Crams end up with "chr" before numbers on output filenames
		File bed = "indexDir/indexDir-indexcov.bed.gz"
		File ped = "indexDir/indexDir-indexcov.ped"
		File roc = "indexDir/indexDir-indexcov.roc"
		File html = "indexDir/index.html"
	}
	runtime {
		docker: "quay.io/aofarrel/goleft-covstats:circleci-push"
		preemptible: indexcovPrempt
		disks: "local-disk " + finalDiskSize + " HDD"
		memory: indexcovMemory + "G"
	}
}

task indexcovBAM {
	input {
		File inputBamOrCram
		Array[File] allInputIndexes

		# runtime attributes
		Int indexcovMemory = 2
		Int indexcovPrempt = 1
		Int indexcovAddlDisk = 0
	}

	command <<<

		set -eux -o pipefail

		# Double-check this is actually a BAM file
		AMIACRAM=$(echo ~{inputBamOrCram} | sed 's/\.[^.]*$//')
		if [ -f ${AMIACRAM}.cram ]; then
			>&2 echo "Cram file detected in the bam task!"
			>&2 echo "This shouldn't happen, please report to the dev."
			exit 1
		
		else
			if [ -f ~{inputBamOrCram}.bai ]; then
				echo "Bai file already exists with pattern *.bam.bai"
			elif [ -f ${AMIACRAM}.bai ]; then
				echo "Bai file already exists with pattern *.bai"
				mv ~{inputBamOrCram}.bai ${AMIACRAM}.bam.bai
			else
				echo "Input bai file not found. We searched for:"
				echo "--------------------"
				echo "  ~{inputBamOrCram}.bai"
				echo "--------------------"
				echo "  ${AMIACRAM}.bai"
				echo "--------------------"
				echo "Finding neither, we will index with samtools."
				samtools index ~{inputBamOrCram} ~{inputBamOrCram}.bai
			fi

			INPUTBAI=$(echo ~{inputBamOrCram}.bai)
			mkdir indexDir
			ln -s ~{inputBamOrCram} indexDir~{basename(inputBamOrCram)}
			ln -s ${INPUTBAI} indexDir~{basename(inputBamOrCram)}.bai
			goleft indexcov --directory indexDir/ *.bam
		fi

	>>>
	# Estimate disk size required
	Int indexSize = ceil(size(allInputIndexes, "GB"))
	Int thisAmSize = ceil(size(inputBamOrCram, "GB"))
	Int finalDiskSize = indexSize + thisAmSize + indexcovAddlDisk
	output {
		# Bams do NOT end up with "chr" before numbers on output filenames
		File bed = "indexDir/indexDir-indexcov.bed.gz"
		File ped = "indexDir/indexDir-indexcov.ped"
		File roc = "indexDir/indexDir-indexcov.roc"
		File html = "indexDir/index.html"
	}
	runtime {
		docker: "quay.io/aofarrel/goleft-covstats:circleci-push"
		preemptible: indexcovPrempt
		disks: "local-disk " + finalDiskSize + " HDD"
		memory: indexcovMemory + "G"
	}
}

task covstats {
	input {
		File inputBamOrCram
		Array[File] allInputIndexes
		File? refGenome

		# runtime attributes with defaults
		Int covstatsMem = 2
		Int covstatsPreempt = 1
		Int covstatsAddlDisk = 0
	}

	command <<<

		echo "Be aware that this does not use the goleft container provided by Biocontainers, which may have some implications for debugging. See README.md on Github for more information."

		start=$SECONDS
		set -eux -o pipefail

		# Detect if inputBamOrCram is a BAM or a CRAM file
		# A CRAM file requires a reference genome but not an index file
		# A BAM file requires an index file but not a reference genome
		AMIACRAM=$(echo ~{inputBamOrCram} | sed 's/\.[^.]*$//')

		if [ -f ${AMIACRAM}.cram ]; then
			echo "Cram file detected"

			# Check if reference genome exists
			if [ "~{refGenome}" != '' ]; then
				goleft covstats -f ~{refGenome} ~{inputBamOrCram} >> this.txt
				# Sometimes this.txt seems to be missing the header... investigate
				COVOUT=$(tail -n +2 this.txt)
				read -a COVARRAY <<< "$COVOUT"
				echo ${COVARRAY[0]} > thisCoverage
				echo ${COVARRAY[11]} > thisReadLength
				BASHFILENAME=$(basename ~{inputBamOrCram})
				echo "'${BASHFILENAME}'" > thisFilename
			
			# Cram file but no reference genome
			else
				>&2 echo "Cram detected but cannot find reference genome."
				>&2 echo "A reference genome is required for cram inputs."
				exit 1
			fi

		else

			# We now know it's a BAM file and must search for an index file
			# or make one ourselves with samtools
			OTHERPOSSIBILITY=$(echo ~{inputBamOrCram} | sed 's/\.[^.]*$//')

			if [ -f ~{inputBamOrCram}.bai ]; then
				echo "Bai file already exists with pattern *.bam.bai"
			elif [ -f ${OTHERPOSSIBILITY}.bai ]; then
				echo "Bai file already exists with pattern *.bai"
			else
				echo "Input bai file not found. We searched for:"
				echo "--------------------"
				echo "  ~{inputBamOrCram}.bai"
				echo "--------------------"
				echo "  ${OTHERPOSSIBILITY}.bai"
				echo "--------------------"
				echo "Finding neither, we will index with samtools."
				samtools index ~{inputBamOrCram} ~{inputBamOrCram}.bai
			fi

			# Not a typo; we don't input the index directly into the call,
			# it just needs to be in the directory
			goleft covstats ~{inputBamOrCram} >> this.txt

			COVOUT=$(tail -n +2 this.txt)
			read -a COVARRAY <<< "$COVOUT"
			echo ${COVARRAY[0]} > thisCoverage
			echo ${COVARRAY[11]} > thisReadLength
			BASHFILENAME=$(basename ~{inputBamOrCram})
			echo "'${BASHFILENAME}'" > thisFilename
		fi

		duration=$(( SECONDS - start ))
		echo ${duration} > duration

	>>>

	# Estimate disk size required
	Int refSize = ceil(size(refGenome, "GB"))
	Int indexSize = ceil(size(allInputIndexes, "GB"))
	Int thisAmSize = ceil(size(inputBamOrCram, "GB"))

	# If input is a cram, it will get samtools'ed into a bam,
	# so we need to at least double its size for the disk
	# calculation.
	Int finalDiskSize = refSize + indexSize + (2*thisAmSize) + covstatsAddlDisk

	output {
		Int outReadLength = read_int("thisReadLength")
		Float outCoverage = read_float("thisCoverage")
		String outFilenames = read_string("thisFilename")
		Int duration = read_int("duration")
	}
	runtime {
		docker: "quay.io/aofarrel/goleft-covstats:circleci-push"
		preemptible: covstatsPreempt
		disks: "local-disk " + finalDiskSize + " HDD"
		memory: covstatsMem + "G"
	}
}

task report {
	input {
		Array[Int] readLengths
		Array[Float] coverages
		Array[String] filenames
		Int lenReads = length(readLengths)
		Int lenCov = length(coverages)
		# runtime attributes
		Int reportMemSize = 2
		Int reportPreempt = 2
	}

	command <<<
	set -eux -o pipefail
	python << CODE
	f = open("reports.txt", "a")
	i = 0

	# if there was just one input, these will not be arrays
	pyReadLengths = ~{sep="," readLengths} # array of ints OR int
	pyCoverages = ~{sep="," coverages} # array of floats OR float
	pyFilenames = ~{sep="," filenames} # array of strings OR string
	
	if (type(pyReadLengths) == int):
		# only one input
		f.write("Filename\tRead length\tCoverage\n")
		f.write("{}\t{}\t{}\n".format(pyFilenames, pyReadLengths, pyCoverages))
		f.close()
	else:
		# print "table" with each inputs' read length and coverage
		f.write("Filename\tRead length\tCoverage\n")
		while i < len(pyReadLengths):
			f.write("{}\t{}\t{}\n".format(pyFilenames[i], pyReadLengths[i], pyCoverages[i]))
			i += 1
		# print average read length
		avgRL = sum(pyReadLengths) / ~{lenReads}
		f.write("Average read length: {}\n".format(avgRL))
		avgCv = sum(pyCoverages) / ~{lenCov}
		f.write("Average coverage: {}\n".format(avgCv))
		f.close()

	CODE
	>>>

	output {
		File finalOut = "reports.txt"
	}

	runtime {
		docker: "python:3.8-slim"
		preemptible: reportPreempt
		memory: reportMemSize + "G"
	}
}

workflow goleft_functions {
	input {
		Array[File] inputBamsOrCrams
		Array[File]? inputIndexes
		File? refGenome
		File? refGenomeIndex # currently unused

		Boolean forceIndexcov = true

		# runtime attributes with defaults
		Int covstatsAddlDisk = 0
		Int covstatsMem = 8
		Int covstatsPreempt = 1
		Int indexcovAddlDisk = 0
		Int indexcovMemory = 2
		Int indexcovPrempt = 1
		Int indexrefAddlDisk = 0
		Int indexrefMem = 2
		Int indexrefPreempt = 1
		Int reportMem = 2
		Int reportPreemptible = 2
	}

	# Weird way to fallback if no indicies are defined, necessary due to how
	# cloud localization works. There may another way to do this, however.
	Array[String] wholeLottaNada = []

	# doesn't work -- is there a WDL builtin that can print?
	#if (length(allIndexes) != length(inputBamsOrCrams)) {
			#echo "Warning: Not all files have an index. If all inputs are BAMs expect a slowdown in covstats."
			#echo "In addition, indexcov will be skipped."
	#}

	if(defined(refGenome)) {
		call indexRefGenome { input: refGenome = refGenome }
	}

	scatter(oneBamOrCram in inputBamsOrCrams) {

		Array[String] allOrNoIndexes = select_first([inputIndexes, wholeLottaNada])

		if (forceIndexcov || length(allOrNoIndexes) == length(inputBamsOrCrams)) {

			String thisFilename = "${basename(oneBamOrCram)}"
			String longerIfACram = sub(thisFilename, "\\.cram", "foobarbizbuzz")
			if (thisFilename == longerIfACram) {
				# only true if running on a BAM with an index, or if forceIndexcov
				call indexcovBAM {
					input:
						inputBamOrCram = oneBamOrCram,
						allInputIndexes = allOrNoIndexes
				}
			}

			if (thisFilename != longerIfACram) {	
				call indexcovCRAM {
					input:
						inputCram = oneBamOrCram,
						allInputIndexes = allOrNoIndexes,
						refGenomeIndex = indexRefGenome.refIndex
				}
			}
		}

		call covstats as scatteredGetStats {
			input:
				inputBamOrCram = oneBamOrCram,
				refGenome = refGenome,
				allInputIndexes = allOrNoIndexes,
				covstatsMem = covstatsMem,
				covstatsAddlDisk = covstatsAddlDisk,
				covstatsPreempt = covstatsPreempt
		}
	}

	call report {
		input:
			readLengths = scatteredGetStats.outReadLength,
			coverages = scatteredGetStats.outCoverage,
			filenames = scatteredGetStats.outFilenames,
			reportMemSize = reportMem,
			reportPreempt = reportPreemptible
	}

	meta {
		author: "Ash O'Farrell"
		email: "aofarrel@ucsc.edu"
		description: "Runs the goleft functions covstats and indexcov, then parses covstats output to extract read length and coverage. If running on Terra download indexDir to read indexcov's HTML output. \nSee README.md on Github for details."
    }
}

