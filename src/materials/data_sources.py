"""External data source clients for compound properties."""

import os
import re
from statistics import mean, mode
from typing import Final

import requests
from bs4 import BeautifulSoup
from chemparse import parse_formula

from .elements import ElementsManager


class PubChemManager:
    """Client for PubChem REST endpoints."""

    URL_SEARCH_BY_FORMULA: Final = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastformula/"
    URL_GET_FROM_WEBPAGE: Final = "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/"

    def __init__(self):
        """Initialize in-memory page cache."""
        self._soup = None
        self._current_cid = None

    def search_cids(self, chemical_formula):
        """Return CIDs for a molecular formula."""
        response = requests.get(self.URL_SEARCH_BY_FORMULA + chemical_formula + "/cids/JSON")
        if response.status_code == 200:
            return response.json()["IdentifierList"]["CID"]
        return []

    def get_properties(self, cid):
        """Fetch compound properties by PubChem CID."""
        self._get_xml_page(cid)
        return {
            "formula": self._get_molecular_formula(),
            "name": self._get_chem_name(),
            "density": self._get_density(),
            "mol_weight": self._get_mol_weight(),
            "id": self._current_cid,
            "source": "PubChem",
        }

    def _get_xml_page(self, cid):
        """Load XML page for a CID into parser cache."""
        response = requests.get(self.URL_GET_FROM_WEBPAGE + str(cid) + "/XML")
        if response.status_code == 200:
            self._soup = BeautifulSoup(response.content, features="xml")
            self._current_cid = cid

    def _get_density(self, few_val_process="mode"):
        """Parse and aggregate density values from PubChem page."""
        str_densities = self._find_strings("Density")
        if str_densities is None:
            return None

        densities = []
        for s in str_densities:
            d = re.findall(r"[\d.]+[\d.\-]*", s.find(string=True))
            dens = list(filter(lambda x: float(x.split("-")[0]) // 100 == 0, d))[0]
            if "-" in dens:
                dens = mean(map(float, dens.split("-")))
            densities.append(float(dens))

        if few_val_process == "mode":
            return mode(densities)
        if few_val_process == "mean":
            return mean(densities)
        return densities[0]

    def _get_mol_weight(self):
        """Return molecular weight value from cached page."""
        str_mol_weight = self._find_strings("Molecular Weight")[0]
        return float(str_mol_weight.find(string=True))

    def _get_molecular_formula(self):
        """Return molecular formula from cached page."""
        return self._find_strings("Molecular Formula")[0].find(string=True)

    def _get_chem_name(self):
        """Return record title from cached page."""
        return self._soup.find("RecordTitle").find(string=True)

    def _find_strings(self, property_name: str):
        """Find all value nodes for given property heading."""
        try:
            strings = self._soup.find("TOCHeading", string=re.compile(property_name)).parent.find_all("String")
        except AttributeError:
            return None
        return strings


class SpringerMaterialsManager:
    """Client for Springer Materials web pages."""

    SM_URL = "https://materials.springer.com/"
    ISP_URL = "isp/crystallographic/docs/"
    package_directory = os.path.dirname(os.path.abspath(__file__))
    PATH_TO_ELEMENT_PROPERTIES = os.path.join(package_directory, "../../substance_properties")
    params = {"searchTerm": "", "datasourceFacet": "sm_isp"}

    def search_sds(self, chemical_formula):
        """Search Springer dataset IDs by chemical formula."""
        self.params["searchTerm"] = chemical_formula
        response = requests.get(self.SM_URL + "search", params=self.params)
        if response.status_code == 200:
            bs = BeautifulSoup(response.content, "html.parser")
            return [a.get("href").split("/")[-1] for a in bs.find_all("a", {"class": "search_result"})]
        return []

    def get_properties(self, sd):
        """Fetch formula, density and molar mass from Springer dataset page."""
        response = requests.get(self.SM_URL + self.ISP_URL + sd)
        if response.status_code == 200:
            bs = BeautifulSoup(response.content, "html.parser")
            props = bs.find("div", {"id": "substance-summary"}).find_all("li")
            dens = props[-1].find("span").text.split(" ")[2]
            st_formula = props[0].find("span").text.strip()
            return {
                "formula": st_formula,
                "name": "",
                "density": dens,
                "mol_weight": self._calc_mol_weight(st_formula),
                "id": sd,
                "source": "Springer Materials",
            }
        return None

    @staticmethod
    def _calc_mol_weight(chem_formula):
        """Calculate molar mass from stoichiometry and atomic weights."""
        parsed_f = parse_formula(chem_formula)
        sum_mol_weight = 0
        for elem in parsed_f.keys():
            sum_mol_weight += ElementsManager(elem).at_weight * parsed_f[elem]
        return sum_mol_weight

