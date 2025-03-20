import json
from abc import ABC, abstractmethod
from typing import Any, Generator

class Dataset(ABC):

    class Sample:
        trees: list[str]
        files: list[str]
        metadata: dict[str, Any]

        def __init__(self):
            # FIXME: will need to make this settable by the dataset implementation
            self.trees = ["events"]
            self.files = []

    _dataset: dict[str, Sample] = {}


    @abstractmethod
    def __init__(self, *args):
        pass



    @staticmethod
    @abstractmethod
    def _parse_path(path: str) -> str:
        pass


    @abstractmethod
    def _get_meta_for_process(self, process_name) -> dict[str, Any]:
        pass


    def _add_process(self, process_name: str):
        meta = self._get_meta_for_process(process_name)
        entry = self.Sample()
        entry.metadata = meta
        self._dataset[process_name] = entry


    def _add_file(self, path) -> None:
        process_name = self._parse_path(path)
        if process_name not in self._dataset:
            self._add_process(process_name)
        self._dataset[process_name].files.append(path)


    def to_json(self, indent: int = 0):
        res = {}
        for name, dataset in self._dataset.items():
            res[name] = vars(dataset)
        return json.dumps({"samples": res}, indent=indent)


    @staticmethod
    def from_json(path: str):
        with open(path) as file:
            dataset = Dataset()
            for name, content in json.load(file)["samples"].items():
                sample = Dataset.Sample()
                sample.trees = content["trees"]
                sample.files = content["files"]
                sample.metadata = content["metadata"]
                dataset._dataset[name] = sample
        return dataset


    def get_weight(self, process_name, int_lumi: float, e_pol: float = 0., p_pol: float = 0.) -> float:
        meta = self._dataset[process_name].metadata
        process_e_pol = meta.get("e_pol", 0.)
        process_p_pol = meta.get("p_pol", 0.)
        pol_weight = 0.25 * (1.0 + e_pol * process_e_pol) * (1.0 + p_pol * process_p_pol)
        n_events = meta["n_events"]
        xsec = meta["xsec_fb"]
        lumi_weight = int_lumi / (n_events / xsec)
        return pol_weight * lumi_weight


    def get_samples(self) -> Generator[str, str, list[str]]:
        for name, sample in self._dataset.items():
            yield name, sample.trees[0], sample.files
