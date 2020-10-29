version 1.0

# Currently assumes only one bam is the input

task index {
	input {
		File inputBam
		String outputBaiString = "${basename(inputBam)}.bai"
	}

	command <<<
		samtools index ~{inputBam} ~{outputBaiString}
	>>>

	output {
		File outputBai = outputBaiString
	}

	runtime {
        docker: "quay.io/biocontainers/goleft:0.2.0--0"
    }
}

task getReadLength {
	input {
		File inputBam
		File inputBai
	}

	command <<<

		# For some reason, a panic in go doesn't exit with status 1, so we
		# have to catch file not found exceptions ourselves
		if [ -f ~{inputBam} ]; then
			echo "Input bam file exists"
		else 
			echo "Input bam file (~{inputBam}) not found, panic"
			exit 1
		fi
		
		# Bai file is NEVER in the same directory as inputBam, trust me on this
		if [ -f ~{inputBai} ]; then
			echo "Input bai file exists"
		else 
			echo "Input bai file (~{inputBam}.bai) not found, panic"
			exit 1
		fi

		# goleft tries to look for the bai in the same folder as the bam, but 
		# they're never in the same folder when run via Cromwell, so we have
		# to symlink it. goleft automatically checks for both name.bam.bai and
		# name.bai so it's okay if we use either 
		inputBamDir=$(dirname ~{inputBam})
		ln -s ~{inputBai} ~{inputBam}.bai
		
		goleft covstats ~{inputBam} >> this.txt
		COVOUT=$(head -2 this.txt | tail -1 this.txt)
		read -a COVARRAY <<< "$COVOUT"
		echo ${COVARRAY[11]} >> readLength

		# clean up
		rm this.txt
	>>>
	output {
		File readLength = "readLength"
	}
	runtime {
        docker: "quay.io/biocontainers/goleft:0.2.0--0"
    }
}

task indexcov {
	input {
		File inputBam
		File inputBai
		#Bool? IncludeGl
		#String? ExcludePatt
		#String? Sex
		#String? Fai
		#String[]? sex
		#exclude
	}

	command <<<

		# For some reason, a panic in go doesn't exit with status 1, so we
		# have to catch file not found exceptions ourselves
		if [ -f ~{inputBam} ]; then
			echo "Input bam file exists"
		else 
			echo "Input bam file (~{inputBam}) not found, panic"
			exit 1
		fi
		
		# Bai file is NEVER in the same directory as inputBam, trust me on this
		if [ -f ~{inputBai} ]; then
			echo "Input bai file exists"
		else 
			echo "Input bai file (~{inputBam}.bai) not found, panic"
			exit 1
		fi

		mkdir indexDir
		ln -s ~{inputBam} indexDir~{basename(inputBam)}
		ln -s ~{inputBai} indexDir~{basename(inputBam)}.bai

		goleft indexcov --directory indexDir/ *.bam

	>>>
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
        docker: "quay.io/biocontainers/goleft:0.2.0--0"
    }
}

workflow goleftwdl {
	input {
		File inputBam
	}

	call index { input: inputBam = inputBam }
	call getReadLength { input: inputBam = inputBam, inputBai = index.outputBai }
	call indexcov { input: inputBam = inputBam, inputBai = index.outputBai }

	meta {
        author: "Ash O'Farrell"
        email: "aofarrel@ucsc.edu"
    }
}
