version 1.0

# Replace the first URL here with the URL of the workflow to be checked.
import "https://raw.githubusercontent.com/aofarrel/goleft-wdl/prefix_outputs/covstats/goleft_functions.wdl" as goleft
import "https://raw.githubusercontent.com/aofarrel/checker-WDL-templates/v0.9.3/tasks/filecheck_task.wdl" as checker_file
import "https://raw.githubusercontent.com/aofarrel/checker-WDL-templates/v0.9.3/tasks/arraycheck_task.wdl" as checker_array

workflow goleft_checker {
	input {
		# These should match the inputs of the workflow being checked
		Boolean? forceIndexcov
		Array[File] inputBamsOrCrams
		Array[File]? inputIndexes
		File? refGenome

		# These are specific to the checker itself
		# Although goleft outputs Array[Array[File]], we're going to put them all in one truth file array
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

	if (defined(goleft_functions.indexcov_of_bams)) {
		scatter (an_output_array in goleft_functions.indexcov_of_bams)
		call checker_array.arraycheck_optional as scatteredChecker {
			input:
				test = an_output_array,
				truth = truth_indexcovBAM
		}
	}

	if (defined(goleft_functions.indexcov_of_crams)) {
		scatter (an_output_array in goleft_functions.indexcov_of_crams)
		call checker_array.arraycheck_optional as scatteredChecker {
			input:
				test = an_output_array,
				truth = truth_indexcovCRAM
		}
	}

	call checker_file.filecheck as singleChecker {
		input:
			test = goleft_functions.covstats_report,
			truth = truth_report
	}

}