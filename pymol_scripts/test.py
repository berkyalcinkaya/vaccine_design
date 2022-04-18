from pymol import cmd

def test(i,s):
    print(s)

cmd.extend("test", test)
