import ROOT
import json
from abc import ABC, abstractmethod
from typing import Any

class Dataset(ABC):

    class DefaultSample:
        trees: list[str] = ["events"]
        files: list[str] = []
        metadata: dict[str, str] = {}

    _dataset: dict[str, DefaultSample] = {}

    # TODO: get_xsec_weight
    # What about the pol weight? Is the logic for this also in
    # the responsibility of this class? Maybe yes, then instead of storing
    # it we can always return 1. for FCC-ee Datasets

    @abstractmethod
    def __init__(self, *args):
        pass


    @staticmethod
    @abstractmethod
    def _parse_path(path: str) -> str:
        pass


    @abstractmethod
    def _get_meta_for_process(self, process_name):
        pass


    def _add_process(self, process_name: str):
        meta = self._get_meta_for_process(process_name)
        entry = self.DefaultSample()
        entry.metadata = meta
        self._dataset[process_name] = entry


    def _add_file(self, path) -> None:
        process_name = self._parse_path(path)
        if not process_name in self._dataset:
            self._add_process(process_name)
        self._dataset[process_name].files.append(path)


    def to_json(self):
        res = {}
        for name, dataset in self._dataset.items():
            res[name] = vars(dataset)
        return json.dumps({"samples": res})


    @abstractmethod
    def get_weight(self, int_lumi: float, *args, **kwargs) -> float:
        return 0.


# TODO: maybe give this its own file
class ILDMC2020Dataset(Dataset):
    """Implementation to create a Dataset from the ILD 2020 mass production

    Relies on the filename scheme of the production.
    """

    _meta_dict: dict[str, Any] = {}
    _name2id: dict[str, str] = {}

    def __init__(self, files: str | list[str], genmeta: str):
        def process_lines(content):
            for line in content:
                path = line.strip()
                self._add_file(path)
        if isinstance(files, str):
            with open(files) as file:
                process_lines(file)
        elif isinstance(files, list):
            process_lines(files)
        else:
            return TypeError

        self._meta_dict = json.loads(genmeta)
        for id, meta in self._meta_dict.items():
            process_name = meta["process_names"]
            e_pol = meta["polarization1"]
            p_pol = meta["polarization2"]
            self._name2id[f"{process_name}_e{e_pol}p{p_pol}"] = id


    @staticmethod
    def _parse_path(path: str) -> str:
        dir, _, fname = path.rpartition("/")
        # the mc-2020 filenames follow a certain fixed naming scheme, example:
        # rv02-02.sv02-02.mILD_l5_o1_v02.E250-SetA.I500102.P4f_sze_sl.eL.pR.n024.d_dstm_15180_122_mini-DST.edm4hep.root
        # <reco v>.<sim v>.<det model>.<machine setting>.<process id>.<process name>.<e pol>.<p pol>.<... other stuff
        parts = fname.split(".")
        process_name = parts[5].lstrip("P")
        e_pol = parts[6]
        p_pol = parts[7]
        return f"{process_name}_{e_pol}{p_pol}"


    def _get_meta_for_process(self, process_name):
        id = self._name2id[process_name]
        meta = self._meta_dict[id]
        # now decide what parts we actually want
        process_meta = {}
        process_meta["xsec_fb"] = meta["cross_section_in_fb"]
        process_meta["n_events"] = 0
        process_e_pol = -1.0 if "eL" in process_name else 1.0 if "eR" in process_name else 0.0
        process_p_pol = -1.0 if "pL" in process_name else 1.0 if "pR" in process_name else 0.0
        process_meta["e_pol"] = process_e_pol
        process_meta["p_pol"] = process_p_pol


    def _add_file(self, path):
        # does all the heavy lifting
        super()._add_file(path)
        # still need to add our events to the count
        with ROOT.TFile(path) as file:
            tree = file["events"]
            meta = self._dataset[self._parse_path(path)].metadata
            meta["n_events"] += tree.GetEntries()


    def get_weight(self, int_lumi, *args, **kwargs):
        e_pol, p_pol = args
        # TODO continue here!!