from dataclasses import dataclass
import logging
import pele_platform.Utilities.Helpers.center_of_mass as cm
import pele_platform.constants.constants as cs

@dataclass
class BoxSetter:
    '''
    Base class to generate the
    exploration box for the simulation
    '''
    box_center: list
    box_radius: float
    ligand_ref: str=""
    logger: logging.Logger=None

    def generate_json(self) -> str:
        # From a box center and box radius generates 
        # json string for PELE control file
        self.logger.info("Generating exploration box")
        self._set_box_center()
        self.logger.info(f"Box with center {self.box_center} \
and radius {self.box_radius} generated\n\n")
        return self.box_json

    def _set_box_center(self):
        # If not user's box center
        # set center of mass of the ligand
        if not self.box_center:
            self.box_center = cm.center_of_mass(self.ligand_ref)

    @property
    def box_json(self):
        # generate the json format for the box
        return cs.BOX.format(self.box_radius, self.box_center) if  self.box_radius else ""

