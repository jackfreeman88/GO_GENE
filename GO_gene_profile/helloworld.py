import cmdlogtime
import sys

# pip install into readme
# stop commenting cant starta again.

COMMAND_LINE_DEF_FILE = "./helloworldCL.txt"
# COMMAND_LINE_DEF_FILE = "./readme_cl.txt"

# There is also an example python file that uses all this stuff in this directory.
def main():
    # Call cmdlogtime.begin() at the beginning of main()
    (start_time_secs, pretty_start_time, my_args, logfile) = cmdlogtime.begin(
        COMMAND_LINE_DEF_FILE, sys.argv[0]
    )

    hello_world(my_args)
    # the command line arguments will be in the my_args dictionary returned, so you can access them like this:
    # get_treats_vectors = my_args["get_treats_vectors"]

    #  Then you put all of your code here.....

    # if you want to add stuff to the logfile:
    # logfile.write("Run Skim here for each of the B terms files in B_TERMS_DIR")

    # The only functions you probably will ever need in cmdlogtime are begin(), end(), and maybe make_dir. Here's an example:
    # cmdlogtime.make_dir(intermediate_out_dir)

    # Call cmdlogtime.end() at the end of main()
    cmdlogtime.end(logfile, start_time_secs)


def get_output(echo_text, shout_text):
    return f"{echo_text}\n{shout_text.upper()}"


def hello_world(my_args):
    output = get_output(my_args["echo_text"], my_args["shout_text"])
    print(output)
    with open(my_args["out_dir"] + "/hw.out", "w") as f:
        f.write(output)


if __name__ == "__main__":
    main()
