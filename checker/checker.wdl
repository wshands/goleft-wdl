version 1.0

import "https://raw.githubusercontent.com/aofarrel/goleft-wdl/prefix_outputs/goleft_functions.wdl" as goleft
import "https://raw.githubusercontent.com/aofarrel/checker-WDL-templates/v0.9.3/tasks/filecheck_task.wdl" as checker_file
import "https://raw.githubusercontent.com/aofarrel/checker-WDL-templates/v0.9.3/tasks/arraycheck_task.wdl" as checker_array

# I don't recommend filling in truth_indexcovBAM and truth_indexcovCRAM with everything in truth_files,
# but in theory you could if you wanted to be a completionist. I only put a handful of these files in
# the public google bucket to prevent people from making ~200 pulls, but you could pull them in from
# a local directory (depending on the backend) or reupload them to your own bucket.

workflow goleft_checker {
	input {
		# These should match the inputs of the workflow being checked
		Boolean? forceIndexcov
		Array[File] inputBamsOrCrams
		Array[File]? inputIndexes
		File? refGenome

		# These are specific to the checker itself
		# Although goleft outputs Array[Array[File]], we put truths in a singular array
		File truth_report
		Array[File] truth_indexcovBAM
		Array[File] truth_indexcovCRAM
	}

	# Run the workflow to be checked
	call goleft.goleft_functions {
		input:
			forceIndexcov = forceIndexcov,
			inputBamsOrCrams = inputBamsOrCrams,
			inputIndexes = inputIndexes,
			refGenome = refGenome
	}

	# Check indexcov of bams
	if (defined(goleft_functions.indexcov_of_bams)) {
		scatter (an_output_array in goleft_functions.indexcov_of_bams) {
			call checker_array.arraycheck_optional as check_indexcov_bams {
				input:
					test = an_output_array,
					truth = truth_indexcovBAM
			}
		}
	}

	# Check indexcov of crams
	if (defined(goleft_functions.indexcov_of_crams)) {
		scatter (an_output_array in goleft_functions.indexcov_of_crams) {
			call checker_array.arraycheck_optional as check_indexcov_crams {
				input:
					test = an_output_array,
					truth = truth_indexcovCRAM
			}
		}
	}

	# Check covstats (same task is used for both crams and bams)
	call checker_file.filecheck as check_covstats {
		input:
			test = goleft_functions.covstats_report,
			truth = truth_report
	}

}