from dataclasses import dataclass

@dataclass
class Nuclide:
    """
    Nuclide dataclass.

    Attributes
    ----------
    symbol : str
        Chemical element symbol (e.g., 'Am').
    Z : int
        Atomic number (e.g., 95).
    A : int
        Mass number (e.g., 242).
    is_fission_product : bool
        Whether the nuclide is a fission product, default False.
    is_metastable : bool
        Whether the nuclide is a metastable/isomeric state, default False.
    """
    symbol: str
    Z: int
    A: int
    is_fission_product: bool = False
    is_metastable: bool = False

    def label(self):
        metastable_str = "m" if self.is_metastable else ""
        return f"{self.symbol}{self.A}{metastable_str}"

@dataclass
class ENDFLibrary:
    name: str
    label: str
    url: str

    def get_filename_of_nuclide(self, nuclide: Nuclide):
        """
        Get the filename of a nuclide in this ENDF library.

        Parameters
        ----------
        nuclide : Nuclide
            The nuclide for which to get the filename.

        Returns
        -------
        str
            The filename of the nuclide in this ENDF library.
        """
        match self.label:
            case "endf71" | "endf80" | "endf81":
                metastable_suffix = "m1" if nuclide.is_metastable else ""
                return f"n-{nuclide.Z:03d}_{nuclide.symbol}_{nuclide.A:03d}{metastable_suffix}.endf"
            case "jeff33":
                metastable_suffix = "m" if nuclide.is_metastable else "g"
                return f"{nuclide.Z}-{nuclide.symbol}-{nuclide.A}{metastable_suffix}.jeff33"
            case "jeff40":
                metastable_suffix = "m" if nuclide.is_metastable else "g"
                return f"n_{nuclide.Z}-{nuclide.symbol}-{nuclide.A}{metastable_suffix}.jeff"
            case _:
                raise ValueError(f"Unknown library label: {self.label}")

endf71_library = ENDFLibrary(
    name="ENDF/B-VIII.0", 
    label="endf71", 
    url="https://www.nndc.bnl.gov/endf-b7.1/zips/ENDF-B-VII.1.tar.gz",
)

endf80_library = ENDFLibrary(
    name="ENDF/B-VIII.0", 
    label="endf80", 
    url="https://www.nndc.bnl.gov/endf-b8.0/zips/ENDF-B-VIII.0.zip",
)

endf81_library = ENDFLibrary(
    name="ENDF/B-VIII.1", 
    label="endf81", 
    url="https://www.nndc.bnl.gov/endf-releases/releases/B-VIII.1/ENDF-B-VIII.1.tar.gz",
)

jeff33_library = ENDFLibrary(
    name="JEFF-3.3", 
    label="jeff33", 
    url="https://www.oecd-nea.org/dbdata/jeff/jeff33/downloads/JEFF33-n.tgz",
)

jeff40_library = ENDFLibrary(
    name="JEFF-4.0", 
    label="jeff40", 
    url="https://data.oecd-nea.org/records/e9ajn-a3p20/files/JEFF40-Evaluations-Neutron-593.zip?download=1",
)


