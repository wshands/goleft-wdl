version 1.0

task indexcov {
	input {
		File inputBamOrCram
		Array[File] allInputIndexes
		File? refGenomeIndex # currently unused

		# runtime attributes
		Int indexcovMemory = 2
		Int indexcovPrempt = 1
		Int indexcovAddlDisk = 0
	}

	command <<<

		set -eux -o pipefail

		AMIACRAM=$(echo ~{inputBamOrCram} | sed 's/\.[^.]*$//')

		if [ -f ${AMIACRAM}.cram ]; then
			>&2 echo "Cram file detected, but crams are currently not supported."
			exit(1)
		else
			if [ -f ~{inputBamOrCram}.bai ]; then
				echo "Bai file already exists with pattern *.bam.bai"
			elif [ -f ${AMIACRAM}.bai ]; then
				echo "Bai file already exists with pattern *.bai"
				mv ~{inputBamOrCram}.bai ${AMIACRAM}.bam.bai
			else
				echo "Input bai file not found. We searched for:"
				echo "  ~{inputBamOrCram}.bai"
				echo "  ${AMIACRAM}.bai"
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
		File depthHTML = "indexDir/indexDir-indexcov-depth-20.html"
		File depthPNG = "indexDir/indexDir-indexcov-depth-20.png"
		File rocHTML = "indexDir/indexDir-indexcov-roc-20.html"
		File rocPNG = "indexDir/indexDir-indexcov-roc-20.png"
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
			goleft covstats -f ~{refGenome} ~{inputBamOrCram} >> this.txt

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

		# runtime attributes with defaults
		Int covstatsAddlDisk = 0
		Int covstatsMem = 8
		Int covstatsPreempt = 1
		Int indexcovAddlDisk = 0
		Int indexcovMemory = 2
		Int indexcovPrempt = 1
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

	scatter(oneBamOrCram in inputBamsOrCrams) {

		Array[String] allOrNoIndexes = select_first([inputIndexes, wholeLottaNada])

		if (length(allOrNoIndexes) == length(inputBamsOrCrams)) {
			String thisFilename = "${basename(oneBamOrCram)}"
			String longerIfACram = sub(thisFilename, "\\.cram", "foobarbizbuzz")
			if (thisFilename == longerIfACram) {
				# Only true if running on a BAM with an index
				call indexcov {
					input:
						inputBamOrCram = oneBamOrCram,
						allInputIndexes = allOrNoIndexes
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
		description: "Runs some goleft functions, namely covstats and indexcov, then parses covstats output to extract read length and coverage."
    }
}

