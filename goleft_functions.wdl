version 1.0

task indexRefGenome {
	# This task only gets called if the refGenome is defined, but we cannot
	# make the refGenome non-optional here unless we want the refGenome to
	# be non-optional for the entire pipeline. We don't want that because
	# the refGenome isn't actually needed at all for some situations.
	input {
		File? refGenome

		# runtime attributes
		Int indexrefMem = 4
		Int indexrefPreempt = 1
		Int indexrefAddlDisk = 1
	}
	# Estimate disk size required
	Int refSize = ceil(size(refGenome, "GB"))
	Int finalDiskSize = 2*refSize + indexrefAddlDisk
	
	command <<<
		ln -s ~{refGenome} .
		samtools faidx ~{basename(select_first([refGenome, 'dummy']))}
	>>>
	
	output {
		File refIndex = glob("*.fai")[0]
	}

	runtime {
		docker: "quay.io/aofarrel/goleft-covstats:circleci-push"
		preemptible: indexrefPreempt
		disks: "local-disk " + finalDiskSize + " HDD"
		memory: indexrefMem + "G"
	}
}

task indexcovCRAM {
	# This task is only called if the user either input a ref genome index or
	# we created one earlier in the indexRefGenome task, so again, we have
	# an "optional" file here that is always going to be defined.
	input {
		File inputCram
		Array[File] allInputIndexes
		File? refGenomeIndex

		# runtime attributes
		Int indexcovMemory = 4
		Int indexcovPrempt = 1
		Int indexcovAddlDisk = 2
	}
	# Estimate disk size required
	Int indexSize = ceil(size(allInputIndexes, "GB"))
	Int thisAmSize = ceil(size(inputCram, "GB"))
	Int finalDiskSize = indexSize + thisAmSize + indexcovAddlDisk

	command <<<
		set -eux -o pipefail

		# Double-check this is actually a cram file
		FILE_EXT=$(echo ~{inputCram} | sed 's/.*\.//')
		FILE_BASE=$(echo ~{inputCram} | sed 's/\.[^.]*$//')
		if [ "$FILE_EXTENSION" = "cram" || "$FILE_EXTENSION" = "CRAM" ]; then
			# Check if an index file for the cram input exists
			if [ -f ~{inputCram}.crai ]; then
				echo "Crai file already exists with pattern *.cram.crai"
			elif [ -f ${FILE_BASE}.crai ]; then
				echo "Crai file already exists with pattern *.crai"
				mv ${FILE_BASE}.crai ${FILE_BASE}.cram.crai
			else
				echo "Input crai file not found. We searched for:"
				echo "--------------------"
				echo "  ~{inputCram}.crai"
				echo "--------------------"
				echo "  ${FILE_BASE}.crai"
				echo "--------------------"
				echo "Finding neither, we will index with samtools."
				samtools index ~{inputCram} ~{inputCram}.crai
			fi

			INPUTCRAI=$(echo ~{inputCram}.crai)
			mkdir indexDir
			ln -s ~{inputCram} indexDir~{basename(inputCram)}
			ln -s ${INPUTCRAI} indexDir~{basename(inputCram)}.crai
			goleft indexcov --extranormalize -d indexDir/ --fai ~{refGenomeIndex} ~{inputCram}.crai

		elif [ -f ${FILE_BASE}.bam ]; then
			>&2 echo "Somehow a bam file got into the cram function!"
			>&2 echo "This shouldn't happen, please report to the dev."
			exit 1
		else
			>&2 echo "Unknown file input, please report to the dev."
		fi

	>>>

	output {
		# Crams end up with "chr" before numbers on output filenames
		Array[File] indexout = glob("indexDir/*")
	}
	
	runtime {
		docker: "quay.io/aofarrel/goleft-covstats:circleci-push"
		preemptible: indexcovPrempt
		disks: "local-disk " + finalDiskSize + " HDD"
		memory: indexcovMemory + "G"
	}
}

task indexcovBAM {
	# Indexcov, when run on bams, doesn't need a refGenome index, but it does need
	# each and every bam to have an index file.
	input {
		File inputBamOrCram
		Array[File] allInputIndexes

		# runtime attributes
		Int indexcovMemory = 4
		Int indexcovPrempt = 1
		Int indexcovAddlDisk = 2
	}

	# Estimate disk size required
	Int indexSize = ceil(size(allInputIndexes, "GB"))
	Int thisAmSize = ceil(size(inputBamOrCram, "GB"))
	Int finalDiskSize = indexSize + thisAmSize + indexcovAddlDisk

	command <<<

		set -eux -o pipefail

		# Double-check this is actually a bam file
		FILE_EXT=$(echo ~{inputCram} | sed 's/.*\.//')
		FILE_BASE=$(echo ~{inputBamOrCram} | sed 's/\.[^.]*$//')
		if [ "$FILE_EXTENSION" = "bam" || "$FILE_EXTENSION" = "BAM" ]; then
			if [ -f ~{inputBamOrCram}.bai ]; then
				echo "Bai file already exists with pattern *.bam.bai"
			elif [ -f ${FILE_BASE}.bai ]; then
				echo "Bai file already exists with pattern *.bai"
				mv ${FILE_BASE}.bai ${FILE_BASE}.bam.bai
			else
				echo "Input bai file not found. We searched for:"
				echo "--------------------"
				echo "  ~{inputBamOrCram}.bai"
				echo "--------------------"
				echo "  ${FILE_BASE}.bai"
				echo "--------------------"
				echo "Finding neither, we will index with samtools."
				samtools index ~{inputBamOrCram} ~{inputBamOrCram}.bai
			fi

			INPUTBAI=$(echo ~{inputBamOrCram}.bai)
			mkdir indexDir
			ln -s ~{inputBamOrCram} indexDir~{basename(inputBamOrCram)}
			ln -s ${INPUTBAI} indexDir~{basename(inputBamOrCram)}.bai
			goleft indexcov --directory indexDir/ *.bam

		elif [ -f ${FILE_BASE}.cram ]; then
			>&2 echo "Cram file detected in the bam task!"
			>&2 echo "This shouldn't happen, please report to the dev."
			exit 1
		else
			>&2 echo "Unknown file input, please report to the dev."
		fi

	>>>
	
	output {
		# Bams do NOT end up with "chr" before numbers on output filenames
		Array[File] indexout = glob("indexDir/*")
	}
	
	runtime {
		docker: "quay.io/aofarrel/goleft-covstats:circleci-push"
		preemptible: indexcovPrempt
		disks: "local-disk " + finalDiskSize + " HDD"
		memory: indexcovMemory + "G"
	}
}

task covstats {
	# The complexity of this task is due to the complexity of input requirements.
	#
	# If the input is a CRAM:
	# * requires a reference genome and WILL fail if not provided
	# * the reference genome must be passed into goleft
	# * does NOT require a cram index file
	#
	# If the input is a BAM:
	# * requires a bam index file but will NOT fail if not provided
	# * the bam index file must be in the working directory
	# * does NOT require a reference genome
	#
	# Furthermore, indexing bams take a long time, so we pass in every index file
	# we have in hopes that one of them will be a matching bam index. If we don't
	# find one then we will index the bam ourselves. This is why we don't fail out.
	# In the cram case, we can't make a reference genome out of thin air, so without
	# a reference genome we must bail. Side note -- neither bams nor crams need a 
	# reference genome index. That's only necessary in indexcov.
	input {
		File inputBamOrCram
		Array[File] allInputIndexes
		File? refGenome

		# runtime attributes
		Int covstatsMem = 4
		Int covstatsPreempt = 2
		Int covstatsAddlDisk = 2
	}
	# Estimate disk size required
	Int refSize = ceil(size(refGenome, "GB"))
	Int indexSize = ceil(size(allInputIndexes, "GB"))
	Int thisAmSize = ceil(size(inputBamOrCram, "GB"))

	# If input is a cram, it will get samtools'ed into a bam, so we need to 
	# account for that. Thankfully a cram should always be smaller than a bam.
	Int finalDiskSize = refSize + indexSize + (2*thisAmSize) + covstatsAddlDisk

	command <<<

		echo "Be aware that this does not use the goleft container provided by Biocontainers,"
		echo "which may have implications for debugging. See README.md on Github for more info."

		start=$SECONDS
		set -eux -o pipefail

		# Detect if inputBamOrCram is a bam or a cram file
		FILE_BASE=$(echo ~{inputBamOrCram} | sed 's/\.[^.]*$//')

		if [ -f ${FILE_BASE}.cram ]; then
			echo "Cram file detected"

			# We have a cram, now check if reference genome exists
			if [ "~{refGenome}" != '' ]; then

				goleft covstats -f ~{refGenome} ~{inputBamOrCram} >> this.txt

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

			# We now know it's a bam file and must search for an index file
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
		Int reportMemSize = 3
		Int reportPreempt = 2
		Int reportDiskSize = 1
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
		disks: "local-disk " + reportDiskSize + " HDD"
		docker: "python:3.8-slim"
		preemptible: reportPreempt
		memory: reportMemSize + "G"
	}
}

workflow goleft_functions {
	input {
		Boolean forceIndexcov = true
		Array[File] inputBamsOrCrams
		Array[File]? inputIndexes
		File? refGenome
		File? refGenomeIndex # currently unused
	}

	Array[String] emptyArray = []

	if(defined(refGenome)) {
		call indexRefGenome { input: refGenome = refGenome }
	}

	scatter(oneBamOrCram in inputBamsOrCrams) {

		Array[String] allOrNoIndexes = select_first([inputIndexes, emptyArray])

		if (forceIndexcov || length(allOrNoIndexes) == length(inputBamsOrCrams)) {

			String thisFilename = "${basename(oneBamOrCram)}"
			String longerIfACram = sub(thisFilename, "\\.cram", "foobarbizbuzz")
			
			if (thisFilename == longerIfACram) {
				# This rings true in the following situations:
				# * This is a bam, and forceIndexCov is true
				# * This is a bam, and we have one index file per input bam/cram input
				# We are hoping that the second case means that the bam has an index file
				# and we won't have to index it ourselves, but this isn't certain
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
				allInputIndexes = allOrNoIndexes
		}
	}

	call report {
		input:
			readLengths = scatteredGetStats.outReadLength,
			coverages = scatteredGetStats.outCoverage,
			filenames = scatteredGetStats.outFilenames
	}

	meta {
		author: "Ash O'Farrell"
		email: "aofarrel@ucsc.edu"
		description: "Runs the goleft functions covstats and indexcov, then parses covstats output to extract read length and coverage. If running on Terra download indexDir to read indexcov's HTML output. See README.md on Github for details."
    }
}

