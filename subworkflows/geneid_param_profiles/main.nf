/*
*  Get parameters module.
*/

/*
 * Defining the output folders.
 */
OutputFolder = "${params.output}"

process paramSplit {

    // where to store the results and in which way
    // publishDir(params.OUTPUT, mode : 'copy', pattern : '*.gff3')

    // indicates to use as a container the value indicated in the parameter
    // container "ferriolcalvet/python-geneid-params"

    // indicates to use as a label the value indicated in the parameter
    label 'geneidx'

    // show in the log which input file is analysed

    input:
    path param_file

    output:
    path ("${param_file_name}.acceptor_profile.param"), emit: acceptor
    path ("${param_file_name}.donor_profile.param"), emit: donor
    path ("${param_file_name}.start_profile.param"), emit: start
    path ("${param_file_name}.stop_profile.param"), emit: stop


    script:
    param_file_name = param_file.getName()
    """
    #!/usr/bin/env python3
    # coding: utf-8

    import re

    pattern = r'([a-zA-Z]+)_[Pp]rofile'

    started = 0
    files_created = []

    # removing the new line characters
    with open("${param_file_name}") as f:
        for line in f:
            if started == 0:
                matches = re.match(pattern, line)
                # if there is any match, get the groups and start a new file
                if matches is not None:
                    groups = matches.groups()
                    # create the file that will contain the given profile
                    # only if it has not been created before
                    filename = "${param_file_name}." + groups[0].lower() + "_profile.param"
                    if filename not in files_created:
                        files_created.append(filename)
                        started = -1
                        fW = open(filename, "w")
                        fW.write(line)

            # we have just found the label of the profile of interest
            # here we need to know the length and order
            elif started == -1:
                # write line to file
                fW.write(line)

                # process the content of this line
                clean_line = line.strip()
                splitted_line = clean_line.split(" ")
                length_pwm = int(splitted_line[0])
                order_pwm = int(splitted_line[3])

                rep_position = 4 ** (order_pwm + 1)
                total_to_read = rep_position * length_pwm

                # define started as the total number of lines to be read and written to the output file
                started = total_to_read + 1

            # these are the lines that we are reading and adding to the output file of the profile
            elif started > 1:
                # write line to file and substract one
                fW.write(line)
                started -= 1
                # print(line, end = '')

            # this is the last line to be added, so we add it and close the file
            elif started == 1:
                # write line to file and substract one
                fW.write(line)
                started -= 1
                # print(line, end = '')

                # close the file as we finished adding that profile
                fW.close()
    """
}

/*
 * Workflow for choosing and providing the parameter file
 */

workflow geneid_param_profiles {

    // definition of input
    take:
    param_file

    main:

    param_file_outs = paramSplit(param_file)

    emit:
    acceptor_pwm = param_file_outs.acceptor
    donor_pwm = param_file_outs.donor
    start_pwm = param_file_outs.start
    stop_pwm = param_file_outs.stop

}
