from dataclasses import dataclass, field
from enum import Enum, auto
from pathlib import Path
import os
from typing import Literal

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
    metastable_state : int | None
        The metastable state of the nuclide, default None.
    """
    symbol: str
    Z: int
    A: int
    is_fission_product: bool = False
    metastable_state: int | None = None

    def label(self):
        metastable_str = f"m{self.metastable_state}" if self.metastable_state is not None else ""
        return f"{self.symbol}{self.A}{metastable_str}"

@dataclass
class ThermalScatteringMaterial:
    """
    Thermal scattering material dataclass.

    Attributes
    ----------
    subject : str
        Scattering atom or species in the material (e.g. 'H', 'D').
    material : str
        Host material name containing the subject (e.g. 'H2O', 'D2O', 'graphite').
    type : str
        Description of the TSL, e.g. 'hh2o' or 'dd2o'
    """
    subject: str
    material: str
    type: str

    def label(self):
        return f"{self.subject}_in_{self.material}"

class DownloadableType(Enum):
    NEUTRONS = auto()
    TSL = auto()

@dataclass
class Downlodable:
    url: str
    filename_archived: str
    filename_extracted: str
    types: list[DownloadableType]

@dataclass
class ENDFLibrary:
    name: str
    label: Literal["endf71", "endf80", "endf81", "jeff33", "jeff40"]
    base_path: Path
    neutrons_path_suffix: str
    tsl_path_suffix: str
    downloadables: list[Downlodable] = field(default_factory=list)

    @property
    def neutrons_path(self) -> Path:
        return self.base_path / self.neutrons_path_suffix

    @property
    def tsl_path(self) -> Path:
        return self.base_path / self.tsl_path_suffix

    def download_library_if_not_exists(self):
        if len(self.downloadables) == 0:
            raise ValueError("No downloadables specified for this library.")

        for downloadable in self.downloadables:
            if downloadable.types == []:
                raise ValueError("Downloadable must have at least one type specified.")
            elif downloadable.types == [DownloadableType.NEUTRONS]:
                if self.neutrons_path.exists():
                    print(f"Neutrons data for {self.name} already exists at '{self.neutrons_path_suffix}', skipping download.")
                    continue
            elif downloadable.types == [DownloadableType.TSL]:
                if self.tsl_path.exists():
                    print(f"Thermal scattering data for {self.name} already exists at '{self.tsl_path_suffix}', skipping download.")
                    continue
            elif set(downloadable.types) == {DownloadableType.NEUTRONS, DownloadableType.TSL}:
                if self.neutrons_path.exists() and self.tsl_path.exists():
                    print(f"Neutrons and thermal scattering data for {self.name} already exist at '{self.neutrons_path}' and '{self.tsl_path}', skipping download.")
                    continue
            else:
                raise ValueError(f"Unsupported combination of DownloadableType in downloadable: {downloadable.types}")

            archive_path = self.base_path / downloadable.filename_archived
            destination_path = self.base_path / downloadable.filename_extracted

            if not self.base_path.exists():
                print(f"Creating base path directory at '{self.base_path}'...")
                self.base_path.mkdir(parents=True, exist_ok=True)

            # Download the library archive
            if archive_path.exists():
                print(f"Archive '{archive_path}' already exists, skipping download.")
            else:
                print(f"Downloading '{downloadable.url}'...")
                os.system(f"wget {downloadable.url} -P {self.base_path}")
            
            # Extract the downloaded file
            print(f"Extracting '{archive_path}' to '{destination_path}'...")
            if archive_path.suffix == ".zip":
                exit_code = os.system(f"unzip -oq {archive_path} -d {destination_path}")
                if exit_code != 0:
                    raise RuntimeError(f"Failed to unzip archive '{archive_path}'.")
            elif archive_path.suffix in [".tar", ".tgz", ".gz"]:
                destination_path.mkdir(parents=True, exist_ok=True)
                exit_code = os.system(f"tar -xzf {archive_path} -C {destination_path} --no-same-owner")
                if exit_code != 0:
                    raise RuntimeError(f"Failed to extract archive '{archive_path}'.")
            else:
                raise ValueError(f"Unsupported archive format: '{archive_path.suffix}'")

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
                metastable_suffix = f"m{nuclide.metastable_state}" if nuclide.metastable_state is not None else ""
                return f"n-{nuclide.Z:03d}_{nuclide.symbol}_{nuclide.A:03d}{metastable_suffix}.endf"
            case "jeff33":
                metastable_suffix = "m" if nuclide.metastable_state is not None else "g"
                return f"{nuclide.Z}-{nuclide.symbol}-{nuclide.A}{metastable_suffix}.jeff33"
            case "jeff40":
                metastable_suffix = "m" if nuclide.metastable_state is not None else "g"
                return f"n_{nuclide.Z}-{nuclide.symbol}-{nuclide.A:03d}{metastable_suffix}.jeff"
            case _:
                raise ValueError(f"Unknown library label: {self.label}")
            
    def get_filename_of_tsl_material(self, tsl_material: ThermalScatteringMaterial):
        """
        Get the filename of a thermal scattering material in this ENDF library.

        Parameters
        ----------
        tsl_material : ThermalScatteringMaterial
            The thermal scattering material for which to get the filename.

        Returns
        -------
        str
            The filename of the thermal scattering material in this ENDF library.
        """
        match self.label:
            case "endf71" | "endf80" | "endf81":
                return f"tsl-{tsl_material.subject}in{tsl_material.material}.endf"
            case "jeff33":
                return f"tsl-{tsl_material.subject}in{tsl_material.material}.jeff33"
            case "jeff40":
                return f"tsl_{tsl_material.subject}_{tsl_material.material}.jeff"
            case _:
                raise ValueError(f"Unknown library label: {self.label}")

    # TODO: Use ENDFtk to parse available temperatures from TSL files instead of hardcoding
    def get_tsl_temperatures(self, tsl_material: ThermalScatteringMaterial) -> list[float]:
        """
        Get the available temperatures for a thermal scattering material in this ENDF library.

        Parameters
        ----------
        tsl_material : ThermalScatteringMaterial
            The thermal scattering material for which to get the temperatures.

        Returns
        -------
        list[float]
            The available temperatures for the thermal scattering material in this ENDF library.
        """
        match self.label:
            case "endf71":
                match tsl_material.label():
                    case "H_in_H2O":
                        return [
                            293.6, 350.0, 400.0, 450.0, 500.0, 550.0, 600.0, 650.0, 800.0,
                        ]
                    case "D_in_D2O":
                        return [
                            293.6, 350.0, 400.0, 450.0, 500.0, 550.0, 600.0, 650.0,
                        ]
                    case _:
                        raise ValueError(f"Unknown TSL material: {tsl_material.label()} for library {self.label}")
            case "endf80" | "endf81":
                match tsl_material.label():
                    case "H_in_H2O":
                        return [
                            283.6, 293.6, 300.0, 323.6, 350.0, 373.6, 400.0, 423.6, 450.0,
                            473.6, 500.0, 523.6, 550.0, 573.6, 600.0, 623.6, 650.0, 800.0,
                        ]
                    case "D_in_D2O":
                        return [
                            283.6, 293.6, 300.0, 323.6, 350.0, 373.6, 400.0, 423.6, 450.0,
                            473.6, 500.0, 523.6, 550.0, 573.6, 600.0, 623.6, 650.0,
                        ]
                    case _:
                        raise ValueError(f"Unknown TSL material: {tsl_material.label()} for library {self.label}")
            case "jeff33" | "jeff40":
                match tsl_material.label():
                    case "H_in_H2O": # jeff40 has additional temperatures, but using only a subset of these in order to be consistent with jeff33
                        return [
                            293.6, 323.6, 373.6, 423.6, 473.6, 523.6, 573.6, 623.6, 647.2, 
                            800.0, 1000.0
                        ]
                    case "D_in_D2O":
                        return [
                            283.6, 293.6, 300.0, 323.6, 350.0, 373.6, 400.0, 423.6, 450.0,
                            473.6, 500.0, 523.6, 550.0, 573.6, 600.0, 623.6, 650.0
                        ]
                    case _:
                        raise ValueError(f"Unknown TSL material: {tsl_material.label()} for library {self.label}")
            case _:
                raise NotImplementedError(f"TSL temperatures not implemented for library: {self.label}")

base_endf_path = Path(__file__).parent / "endf_libraries"

endf71_library = ENDFLibrary(
    name="ENDF/B-VIII.0", 
    label="endf71", 
    base_path=base_endf_path / "endf71",
    neutrons_path_suffix="ENDF-B-VII.1/ENDF-B-VII.1-neutrons/neutrons",
    tsl_path_suffix="ENDF-B-VII.1/ENDF-B-VII.1-thermal_scatt/thermal_scatt",
    downloadables=[
        Downlodable(
            url="https://www.nndc.bnl.gov/endf-b7.1/zips/ENDF-B-VII.1.tar.gz",
            filename_archived="ENDF-B-VII.1.tar.gz",
            filename_extracted="ENDF-B-VII.1",
            types=[DownloadableType.NEUTRONS, DownloadableType.TSL],
        )
    ]
)

endf80_library = ENDFLibrary(
    name="ENDF/B-VIII.0", 
    label="endf80", 
    base_path=base_endf_path / "endf80",
    neutrons_path_suffix="ENDF-B-VIII.0/ENDF-B-VIII.0/neutrons",
    tsl_path_suffix="ENDF-B-VIII.0/ENDF-B-VIII.0/thermal_scatt",
    downloadables=[
        Downlodable(
            url="https://www.nndc.bnl.gov/endf-b8.0/zips/ENDF-B-VIII.0.zip",
            filename_archived="ENDF-B-VIII.0.zip",
            filename_extracted="ENDF-B-VIII.0",
            types=[DownloadableType.NEUTRONS, DownloadableType.TSL],
        )
    ]
)

endf81_library = ENDFLibrary(
    name="ENDF/B-VIII.1", 
    label="endf81", 
    base_path=base_endf_path / "endf81",
    neutrons_path_suffix="ENDF-B-VIII.1/neutrons-version.VIII.1",
    tsl_path_suffix="ENDF-B-VIII.1/thermal_scatt-version.VIII.1",
    downloadables=[
        Downlodable(
            url="https://www.nndc.bnl.gov/endf-releases/releases/B-VIII.1/ENDF-B-VIII.1.tar.gz",
            filename_archived="ENDF-B-VIII.1.tar.gz",
            filename_extracted="ENDF-B-VIII.1",
            types=[DownloadableType.NEUTRONS, DownloadableType.TSL],
        )
    ]
)

jeff33_library = ENDFLibrary(
    name="JEFF-3.3", 
    label="jeff33",
    base_path=base_endf_path / "jeff33",
    neutrons_path_suffix="JEFF33-n/endf6",
    tsl_path_suffix="JEFF33-tsl/JEFF33-tsl",
    downloadables=[
        Downlodable(
            url="https://www.oecd-nea.org/dbdata/jeff/jeff33/downloads/JEFF33-n.tgz",
            filename_archived="JEFF33-n.tgz",
            filename_extracted="JEFF33-n",
            types=[DownloadableType.NEUTRONS],
        ),
        Downlodable(
            url="https://www.oecd-nea.org/dbdata/jeff/jeff33/downloads/JEFF33-tsl.tgz",
            filename_archived="JEFF33-tsl.tgz",
            filename_extracted="JEFF33-tsl",
            types=[DownloadableType.TSL],
        ),
    ]
)

jeff40_library = ENDFLibrary(
    name="JEFF-4.0", 
    label="jeff40", 
    base_path=base_endf_path / "jeff40",
    neutrons_path_suffix="JEFF40-Evaluations-Neutron-593",
    tsl_path_suffix="JEFF40-Evaluations-TSL",
    downloadables=[
        Downlodable(
            url="https://data.oecd-nea.org/records/e9ajn-a3p20/files/JEFF40-Evaluations-Neutron-593.zip",
            filename_archived="JEFF40-Evaluations-Neutron-593.zip",
            filename_extracted="JEFF40-Evaluations-Neutron-593",
            types=[DownloadableType.NEUTRONS],
        ),
        Downlodable(
            url="https://data.oecd-nea.org/records/scmva-nqh68/files/JEFF40-Evaluations-TSL.zip",
            filename_archived="JEFF40-Evaluations-TSL.zip",
            filename_extracted="JEFF40-Evaluations-TSL",
            types=[DownloadableType.TSL],
        ),
    ]
)


