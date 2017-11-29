import sys


class input_fixer(object):

    def __init__(self, input):
        self.input = input
        self.strip_waters()

    def strip_waters(self):
        input_fix = []
        with open(self.input, 'r') as pdb:
            for line in pdb:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    if line[17:20] != "HOH":
                        input_fix.append(line)
                    else:
                        continue
                else:
                    input_fix.append(line)

        with open(self.input, 'w') as fix_pdb:
            fix_pdb.write("".join(input_fix))


if __name__ == '__main__':
    fixer = input_fixer(sys.argv[1])