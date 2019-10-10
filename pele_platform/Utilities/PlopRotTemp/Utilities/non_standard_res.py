




def remove_capping_hidrogens(output_template_file, nstd):
        """
        Remove the capping hidrogens of the template
        given as a paramaters for the user when
        dealing with non standard residues
        """
	
        with open(output_template_file, "r") as f:
            lines = f.readlines()
        
        #Create Variables
        new_lines = lines
        atom_index = [line[0:5].strip() for line in lines if line[21:25] in nstd]
        index_to_remove = []
        fields = {
        "NAME" : False,
        "NBON" : False,
        "BOND" : False,
        "THET" : False,
        "PHI" : False,
        "IPHI" : False
        }

        #Remove lines from capping atoms
        for i, line in enumerate(lines):
            found=False
            for value in ["NAME", "NBON", "BOND", "THET", "PHI", "IPHI"]:
                if line.strip("\n") == value:
                    fields[value] = True
                    found=True
                if found:
                    found=False
                    continue

                if i<=2 and not fields["NBON"] and not fields["BOND"] and not fields["THET"] and not fields["PHI"] and not fields["IPHI"]:
                    if i == 2:
                        new_lines[i] = line[0:9] + str(int(line[9:11])-len(atom_index)) + line[11:]
                    else:
                        pass	
                elif i>2 and not fields["NBON"] and not fields["BOND"] and not fields["THET"] and not fields["PHI"] and not fields["IPHI"]:
                        if line[21:25].strip() in nstd:
                            index_to_remove.append(i)
                        else:
                            pass

                elif fields["NBON"] and not fields["BOND"] and not fields["THET"] and not fields["PHI"] and not fields["IPHI"]:
                        if line[0:6].strip() in atom_index:
                            index_to_remove.append(i)
                        else:
                            pass
                elif fields["NBON"] and fields["BOND"] and not fields["THET"] and not fields["PHI"] and not fields["IPHI"]:
                        atom_1, atom_2 = line.split()[0:2]
                        if atom_1 in atom_index or atom_2 in atom_index:
                            index_to_remove.append(i)
                        else:
                            pass
                elif fields["NBON"] and fields["BOND"] and fields["THET"] and not fields["PHI"] and not fields["IPHI"]:
                        atom_1, atom_2, atom_3 = line.split()[0:3]
                        if atom_1 in atom_index or atom_2 in atom_index or atom_3 in atom_index:
                            index_to_remove.append(i)
                        else:
                            pass
                elif fields["NBON"] and  fields["BOND"] and fields["THET"] and  fields["PHI"] and not fields["IPHI"]:
                        atom_1, atom_2, atom_3, atom_4 = line.split()[0:4]
                        if atom_1 in atom_index or atom_2 in atom_index or atom_3 in atom_index or atom_4 in atom_index:
                            index_to_remove.append(i)
                        else:
                            pass
                elif fields["NBON"] and  fields["BOND"] and fields["THET"] and  fields["PHI"] and  fields["IPHI"] and line != "END":
                        atom_1, atom_2, atom_3, atom_4 = line.split()[0:4]
                        if atom_1 in atom_index or atom_2 in atom_index or atom_3 in atom_index or atom_4 in atom_index:
                            index_to_remove.append(i)
                        else:
                            pass
        #Write back
        template = [line for i, line in enumerate(new_lines) if i not in index_to_remove]
        with open(output_template_file, "w") as f:
            f.write("".join(template))
