#!/usr/bin/python3
# -*- coding: utf-8 -*-

# Version incopatibilities
try:
    import Tkinter as tk
    import tkFileDialog as filedialog
    import ttk
except ImportError:
    import tkinter as tk
    from tkinter import filedialog
    from tkinter import ttk
import os
import re
import subprocess


# class Navbar(tk.Frame):
# #GUI Init
# tk.Frame.__init__(self, parent, *args, **kwargs)

# #Fill GUI



class Controller(object):

    def __init__(self, gui, model, *args, **kwargs):
        self.gui = gui
        self.model = model
        self.set_mvc()

    def set_mvc(self):
        self.model.run_button.config(command=self.run)

    def run(self):
        if not self.model.input:
            return

        self.command = self.model.create_command()

        print(self.command.split())
        self.subprocess = subprocess.call(self.command.split())


class Model(object):

    def __init__(self, gui, *args, **kwargs):
        self.gui = gui

    def create_command(self):
        options = [self.mae_charges,
                   self.clean,
                   self.max_torsions,
                   self.sidechains,
                   self.atom_core,
                   self.cpus,
                   self.confile,
                   self.native,
                   self.forcefield,
                   ]

        options = [option for option in options if option]
        self.command = 'bash PelePlop.sh {} {} {} {}'.format(
            ' '.join(options), self.input, self.residue, self.chain)
        return self.command

    @property
    def run_button(self):
        return self.gui.run_but

    @property
    def input(self):
        return self.gui.system.var_input_path.get()

    @property
    def residue(self):
        return self.gui.ligand.var_residues.get()

    @property
    def chain(self):
        return self.gui.ligand.var_chain.get()

    @property
    def mae_charges(self):
        mae_charges = self.gui.options.var_mae_charges.get()
        return mae_charges if mae_charges else None

    @property
    def clean(self):
        clean = self.gui.options.var_clean.get()
        return clean if clean else None

    @property
    def max_torsions(self):
        max_tors = self.gui.options.var_mtor.get()
        return "--mtor {}".format(max_tors) if max_tors else None

    @property
    def sidechains(self):
        sidechains = self.gui.options.var_sidechains.get()
        return "--n {}".format(sidechains) if sidechains else None

    @property
    def atom_core(self):
        atom_core = self.gui.options.var_core_atom_value.get()
        return "--core {}".format(atom_core)if atom_core else None

    @property
    def cpus(self):
        cpus = self.gui.options.var_cpus.get()
        return "--cpus {}".format(cpus) if cpus else None

    @property
    def confile(self):
        confile = self.gui.options.var_conf_file.get()
        return "--confile {}".format(confile) if confile else None

    @property
    def native(self):
        native = self.gui.options.var_native_path.get()
        return "--native {}".format(native) if native else None

    @property
    def forcefield(self):
        forcefield = self.gui.options.var_forcefield.get()
        return "--forcefield {}".format(forcefield) if forcefield else None
        
STYLES = {
    tk.Entry: {
        'background': 'white',
        'borderwidth': 1,
        'highlightthickness': 0,
        'width': 10,
    },
    tk.Listbox: {
        'height': '10',
        'width': '5',
        'background': 'white',

    },
    tk.Button: {
        'borderwidth': 1,
        'highlightthickness': 0,

    },
    tk.Checkbutton: {
        #'highlightbackground': chimera.tkgui.app.cget('bg'),
        #'activebackground': chimera.tkgui.app.cget('bg'),
    }
}

class PelePlopDialog(tk.Frame):

    """
    Displays main GUI and initializes models and controllers
    for the respective file format.
    """

    help = "https://github.com/miniaoshi/PelePlop_gui"
    VERSION = '0.0.1'
    EXIT = "PelePlopExited"

    def __init__(self, parent, *args, **kwargs):

        # GUI Init
        self.main_frame = ttk.Frame(parent, padding="3 20 3 20", borderwidth=2, relief="solid")
        self.main_frame.pack(expand=True, fill='both')
        self.main_frame.columnconfigure(0, weight=1)
        # self.main_frame.rowconfigure(0, weight=1)

        # main_frames config
        parent.title("PelePlop")
        self.run_but = tk.Button(self.main_frame, text="Run")
        self.run_but.grid(column=0, row=4)

        # Top windows
        self.system = System(self.main_frame)
        self.ligand = Ligand(self.main_frame)
        self.options = Options(self.main_frame)

        # Callbacks
        self.main_frame.bind('<Return>', (lambda e, b=self.run_but: b.invoke()))  # b is your button
        (self.system.var_input_path).trace_variable(
            "w", lambda x, y, z: _update_lig_resbox(
                self.system.var_input_path.get(),
                self.ligand.residue_combobox))

        (self.ligand.var_residues).trace_variable(
            "w", lambda x, y, z: _update_lig_chainbox(
                self.system.var_input_path.get(),
                self.ligand.residue_combobox.get(),
                self.ligand.chain_combobox))


class System(tk.Frame):

    """
        Display the widgets to
        recolect PELE input system.
    """

    def __init__(self, parent, *args, **kwargs):

        # High level variables
        self.var_input_path = tk.StringVar()

        # Create system frames
        system_frame = ttk.Frame(parent, borderwidth=2, relief="sunken")
        system_frame.grid(column=0, row=0, padx=5, pady=30)

        # Fill in system friend
        self.input_title = tk.Label(system_frame, text="System", font=20)
        self.input_button = tk.Label(system_frame, text="Input File ")
        self.input_entry = tk.Entry(system_frame, textvariable=self.var_input_path)
        self.input_search = tk.Button(system_frame, text='...',
                                      command=lambda: _browse_file(self.var_input_path, "*.pdb", "*.mae"))

        # Grid widgets
        self.input_title.grid(row=0, column=1, pady=20)
        self.input_button.grid(row=1, column=0, sticky="ew")
        self.input_entry.grid(row=1, column=1, sticky="ew")
        self.input_search.grid(row=1, column=2, sticky="ew")


class Ligand(tk.Frame):

    """
        Display the widgets to
        recapt Ligand Information
        to later separate protein
        from ligand.
    """

    def __init__(self, parent, *args, **kwargs):

        # High Level variables
        self.var_residues = tk.StringVar()
        self.var_chain = tk.StringVar()

        # Create Frame
        ligand_frame = ttk.Frame(parent, borderwidth=2, relief="sunken")
        ligand_frame.grid(column=0, row=1, padx=5, pady=20)

        # Widgets
        self.ligand_title = tk.Label(ligand_frame, text="Ligand", font=20)
        self.residue_label = tk.Label(ligand_frame, text="Ligand Residue")
        self.residue_combobox = ttk.Combobox(ligand_frame, textvariable=self.var_residues)
        self.chain_label = tk.Label(ligand_frame, text="Chain Residue")
        self.chain_combobox = ttk.Combobox(ligand_frame, textvariable=self.var_chain)

        # Grid Widgets
        self.ligand_title.grid(row=0, column=0, columnspan=2, pady=20)
        self.residue_label.grid(row=1, column=0, sticky="ew")
        self.residue_combobox.grid(row=1, column=1, sticky="ew")
        self.chain_label.grid(row=2, column=0, sticky="ew", pady=10)
        self.chain_combobox.grid(row=2, column=1, sticky="ew")


class Options(tk.Frame):

    """
        Display PlopRotTemp Options
    """

    def __init__(self, parent, *args, **kwargs):

        # High level Frame Variables
        self.var_mae_charges = tk.StringVar()
        self.var_clean = tk.StringVar()
        self.var_mtor = tk.StringVar()
        self.var_sidechains = tk.StringVar()
        self.var_core_atom_bool = tk.StringVar()
        self.var_core_atom_value = tk.StringVar()
        self.var_native_path = tk.StringVar()
        self.var_cpus = tk.StringVar()
        self.var_conf_file = tk.StringVar()
        self.var_forcefield = tk.StringVar()

        # Defaults
        self.var_mtor.set(4)
        self.var_sidechains.set(0)
        self.var_cpus.set(3)

        # Frame
        options_frame = ttk.Frame(parent, borderwidth=2, relief="sunken")
        options_frame.grid(column=0, row=2, padx=5, pady=20)

        # Widgets
        self.options_title = tk.Label(options_frame, text="Advanced Options", font=20)
        self.mae_charges_checkbut = ttk.Checkbutton(
            options_frame, text='Keep Charges From Mae File', variable=self.var_mae_charges,
            onvalue='--mae_charges', offvalue='')
        self.clean_checkbut = ttk.Checkbutton(
            options_frame, text='Clean Residual Files', variable=self.var_clean,
            onvalue='--clean', offvalue='')
        self.max_torsions_label = tk.Label(options_frame, text="Max Torsions")
        self.max_torsions_entry = tk.Entry(options_frame, textvariable=self.var_mtor)
        self.max_groups_label = tk.Label(options_frame, text="Max Sidechains")
        self.max_groups_entry = tk.Entry(options_frame, textvariable=self.var_sidechains)
        self.core_atom_entry = tk.Entry(options_frame, textvariable=self.var_core_atom_value, state='disabled')
        self.core_checkbut = ttk.Checkbutton(
            options_frame, text='Set Core Atom', variable=self.var_core_atom_bool,
            command=lambda: _enable(self.var_core_atom_bool, self.core_atom_entry))
        self.cpu_label = tk.Label(options_frame, text="CPUs")
        self.cpu_entry = tk.Entry(options_frame, textvariable=self.var_cpus)
        self.native_label = tk.Label(options_frame, text="Native File ")
        self.native_entry = tk.Entry(options_frame, textvariable=self.var_native_path)
        self.native_search = tk.Button(options_frame, text='...',
                                       command=lambda: _browse_file(self.var_native_path, "*.pdb", "*.mae"))
        self.conf_file_label=tk.Label(options_frame, text="Configuration File")
        self.conf_file_entry = tk.Entry(options_frame, textvariable=self.var_conf_file)
        self.conf_file_search = tk.Button(options_frame, text='...',
                                       command=lambda: _browse_file(self.var_conf_file, "*", "*.json"))
        self.forcefield_label = tk.Label(options_frame, text="Forcefield Residue")
        self.forcefield_combobox = ttk.Combobox(options_frame, textvariable=self.var_forcefield)
        self.forcefield_combobox['values'] = ["OPLS2005", "AMBER99sb", "AMBER99sbBSC0"]
        self.forcefield_combobox.current(0)

        # Grid Widgets
        self.options_title.grid(row=0, column=0, columnspan=2, pady=20)
        self.mae_charges_checkbut.grid(row=1, column=0, columnspan=2)
        self.clean_checkbut.grid(row=2, column=0, columnspan=2)
        self.max_torsions_label.grid(row=3, column=0, sticky="ew")
        self.core_checkbut.grid(row=5, column=0)
        self.max_torsions_entry.grid(row=3, column=1, sticky="ew")
        self.max_groups_label.grid(row=4, column=0, sticky="ew")
        self.max_groups_entry.grid(row=4, column=1, sticky="ew")
        self.core_atom_entry.grid(row=5, column=1, sticky="ew")
        self.cpu_label.grid(row=6, column=0, sticky="ew")
        self.cpu_entry.grid(row=6, column=1, sticky="ew")
        self.native_label.grid(row=7, column=0, sticky="ew")
        self.native_entry.grid(row=7, column=1, sticky="ew")
        self.native_search.grid(row=7, column=2, sticky="ew")
        self.conf_file_label.grid(row=8, column=0, sticky="ew")
        self.conf_file_entry.grid(row=8, column=1, sticky="ew")
        self.conf_file_search.grid(row=8, column=2, sticky="ew")
        self.forcefield_label.grid(row=9, column=0, sticky="ew")
        self.forcefield_combobox.grid(row=9, column=1, sticky="ew")


def _browse_file(var_store_path, file_type1, file_type2):
    """
    Browse file path
    Parameters
    ----------
    var= Interface entry widget where we wish insert the browse file.
    file_type1 = 1st type of file to open
    file_type2 = 2nd  type of file to open
    """

    path = filedialog.askopenfilename(initialdir='~/', filetypes=(
        (file_type1, file_type1), (file_type2, file_type2)))

    if path:
        var_store_path.set(path)


def _enable(var, *entries):
    """
    Enable or Disable several settings
    depending on other Checkbutton value
    Parameters
    ----------
    var: tk widget Checkbutton where we set an options
    onvalue: onvalue Checkbutton normally set as 1
    args...: tk widgets to enable or disabled
    """

    for entry in entries:
        entry['state'] = 'normal' if var.get() == '1' else 'disabled'


def _get_res_and_chain_from_pdb(pdbfile):
    """
        Get all residues and its chain from pdb input file
    """
    res_and_chains = []
    if os.path.isfile(pdbfile):
        with open(pdbfile, "r") as f:
            lines = f.readlines()
            lines = [line.strip('\n').strip() for line in lines]
        atoms = [line for line in lines if line.startswith("ATOM") or line.startswith("HETATM")]

        for atom in atoms:
            # if [residue,chain] not in res_and_chains
            if [atom[16:21].strip(), atom[21]] not in res_and_chains:
                res_and_chains.append([atom[16:21].strip(), atom[21]])
        return res_and_chains
    else:
        return False


def _update_lig_resbox(pdbfile, combobox):
    """
        update the combobox residues
        from the ligand section with
        the residues of the pdb selected
        as input file
    """

    res_and_chains = _get_res_and_chain_from_pdb(pdbfile)
    if res_and_chains:
        combobox['values'] = [res for res, chain in res_and_chains]
        combobox.current(0)


def _update_lig_chainbox(pdbfile, residue_selected, combobox_chain):
    """
        update the combobox chains
        from the ligand section with
        all the possible chain of the
        residue selected as ligand
    """
    res_and_chains = _get_res_and_chain_from_pdb(pdbfile)
    if res_and_chains:
        chains = []
        for res, chain in res_and_chains:
            if(chain not in chains and
                    res == residue_selected):
                chains.append(chain)
        combobox_chain['values'] = chains
        combobox_chain.current(0)


if __name__ == '__main__':
    root = tk.Tk()
    gui = PelePlopDialog(root)
    model = Model(gui)
    Controller(gui, model)
    root.mainloop()
