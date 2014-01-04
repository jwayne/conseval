#!/usr/bin/python
import sys

COMMANDS = set(("single", "batch", "evaluate", "params"))

def print_help():
    print "usage: %s <command> [<args>]" % sys.argv[0]
    print ""
    cmds = list(COMMANDS)
    cmds.append("help")
    print "Available commands: %s" % ", ".join(sorted(cmds))
    sys.exit(0)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print_help()

    command = sys.argv[1]
    if command == "--help" or command == "help":
        print_help()
        
    if command not in COMMANDS:
        print "%s: '%s' is not a valid command. See '%s --help'." \
            % (sys.argv[0], sys.argv[1], sys.argv[0])
        sys.exit(0)

    module = __import__(command)
    args = [" ".join(sys.argv[0:2])]
    args += sys.argv[2:]
    sys.argv = args
    module.main()
