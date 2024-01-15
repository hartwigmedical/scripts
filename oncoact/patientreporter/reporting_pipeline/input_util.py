def multi_line_input(msg):
    """
    Allows the user to enter multiple lines for input. If the user enters nothing, it will return None.
    """
    print(msg, "(Type END when you are finished)\n")
    res = []
    while True:
        next_line = input()
        if next_line == '' and len(res) == 0:
            return None
        if next_line == 'END':
            return '\n'.join(res)
        res.append(next_line)