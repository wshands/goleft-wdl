version 1.0

# Replace the first URL here with the URL of the workflow to be checked.
import "https://raw.githubusercontent.com/aofarrel/goleft-wdl/master/covstats/goleft_functions.wdl" as goleft
import "https://raw.githubusercontent.com/aofarrel/checker-WDL-templates/v0.9.3/tasks/filecheck_task.wdl" as checker_file
import "https://raw.githubusercontent.com/aofarrel/checker-WDL-templates/v0.9.3/tasks/arraycheck_task.wdl" as checker_array

workflow checker {
	input {
		# These should match the inputs of the workflow being checked
		Boolean? forceIndexcov
		Array[File] inputBamsOrCrams
		Array[File]? inputIndexes
		File? refGenome

		# These are specific to the checker itself
		# Although goleft outputs Array[Array[File]], we're going to put them all in one truth file array
		# This unfortunately will mean that we can't have duplicated file names...
		# Maybe direct users to use prefixes for their outputs in their tasks?
		File truth_report
		Array[File] truth_indexcovBAM
		Array[File] truth_indexcovCRAM
	}

	call blank

	# Run the workflow to be checked
	call check_me.goleft_functions {
		input:
			forceIndexcov = forceIndexcov,
			inputBamsOrCrams = inputBamsOrCrams,
			inputIndexes = inputIndexes,
			refGenome = refGenome
	}

	# Check an array of files, wherein SOME of the files in that array might not be defined
	# Any files that might not be defined need to fall back on a file that does exist, which
	# can be done easily by passing in a bogus file via select_first. This bogus file will
	# not have a match in the truth array, so it won't get md5 checked.
	call checker_array.arraycheck_classic as nonscatteredChecker {
		input:
			test = [goleft_functions.wf_always, select_first([goleft_functions.wf_never, blank.bogus]), select_first([goleft_functions.wf_sometimesSingle, blank.bogus])],
			truth = arrayTruth
	}

	# Check an array of files, wherein the ENTIRE array might not be defined
	# In this example, the output of sometimesScattered is multiple files with the same name
	if (defined(goleft_functions.indexcov_of_bams)) {
		scatter (an_output_array in goleft_functions.indexcov_of_bams)
		call checker_array.arraycheck_optional as scatteredChecker {
			input:
				test = an_output_array,
				truth = arrayTruth
		}
	}

	# Here, we include a check for a single optional file. In the original workflow we implied
	# it is never created, but in practice you can find use for this in files that only created
	# if the user specifies a certain input flag. By including a defined() check beforehand,
	# we avoid actually calling the task unless the test file actually exists, which will save
	# us time and compute credits. Without this defined() check, the task will spin up as the
	# task considers test to be an optional input, but will fail upon execution if test does not
	# exist due to how the task is written.
	if (defined(goleft_functions.wf_never)) {
		call checker_file.filecheck as singleChecker {
			input:
				test = goleft_functions.wf_never,
				truth = singleTruth
		}
	}

}

task blank {
	input {
		Int? nada
	}

	command <<<
		touch "dummy_file.txt"
	>>>
	
	output {
		File bogus = "dummy_file.txt"
	}

	runtime {
		docker: "quay.io/aofarrel/goleft-covstats:circleci-push"
		preemptible: 3
		memory: 2 + "G"
	}
}
