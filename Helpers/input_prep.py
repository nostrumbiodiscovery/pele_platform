import sys


class input_fixer(object):

    def __init__(self, input):
        self.input = input
        self.fix_input()

    def fix_input(self):
        
        """
        1- Strip Waters
        2- Align Atoms
        """

        input_fix = []
        with open(self.input, 'r') as pdb:
            for line in pdb:
                if line.startswith("ATOM"):
                    #Strip waters
                    if line[17:20] != "HOH":
                        input_fix.append(line)
                elif line.startswith("HETATM"):
                    #Strip waters and align atoms
                    if line[17:20] != "HOH":
                        line = self.align_atomnames(line)
                        input_fix.append(line)
                else:
                    #Nothing
                    input_fix.append(line)

        with open(self.input, 'w') as fix_pdb:
            fix_pdb.write("".join(input_fix))

    @staticmethod
    def align_atomnames(line):
        
        """
          Align atom names to the right to match
          PlopRoTemp Template
        
        """
        if line[12]!= " " and line[15] == " ": 
            new_line = "{} {}{}".format(line[0:12], line[12:16], line[17:])
            return new_line
        else:
            return line




        

if __name__ == '__main__':
    fixer = input_fixer(sys.argv[1])