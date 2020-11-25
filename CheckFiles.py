'''
    Alberto Roldan

    Check if the input files are readable

    STRUCTURE:
        - Read geometry file
        - Read energy file

'''

def checkFiles(Files):
    for f in Files:
        try:
            fd = open(f)
        except IOError:
            raise Exception("   "+f+" is not a valid  file!")


