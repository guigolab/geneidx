import re
import sys


param_file_name = sys.argv[1]
pattern = r'([a-zA-Z]+)_[Pp]rofile'

started = 0
files_created = []

print(param_file_name)
# removing the new line characters
with open(f"{param_file_name}.param") as f:
    for line in f:
        if started == 0:
            matches = re.match(pattern, line)
            # if there is any match, get the groups and start a new file
            if matches is not None:
                groups = matches.groups()
                # create the file that will contain the given profile
                # only if it has not been created before
                filename = f"{param_file_name}." + groups[0].lower() + "_profile.param"
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