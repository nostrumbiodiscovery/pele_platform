import os
from string import Template


class TemplateBuilder(object):

    def __init__(self, file, keywords):

        self.file = file
        self.keywords = {k: v for k, v in keywords.items() if v is not None}

        self.fill_in()

    def fill_in(self):
        """
        Fill the control file in
        """
        with open(os.path.join(self.file), 'r') as infile:
            confile_data = infile.read()

        confile_template = Template(confile_data)

        confile_text = confile_template.safe_substitute(self.keywords)

        with open(os.path.join(self.file), 'w') as outfile:
            outfile.write(confile_text)
