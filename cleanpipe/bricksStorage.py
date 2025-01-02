



def txt2list(s_filename):
    """
    theoutput is a list where each item is a string of one line
    """
    with open(s_filename, 'r') as file:
        l_lines = file.readlines()
    return l_lines