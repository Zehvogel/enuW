import json
from abc import ABC, abstractmethod
from typing import Any, Generator

class Dataset(ABC):

    class Sample:
        trees: list[str] = ["events"]
        files: list[str] = []
        metadata: dict[str, Any] = {}

    _dataset: dict[str, Sample] = {}


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
        entry = self.Sample()
        entry.metadata = meta
        self._dataset[process_name] = entry


    def _add_file(self, path) -> None:
        process_name = self._parse_path(path)
        if process_name not in self._dataset:
            self._add_process(process_name)
        self._dataset[process_name].files.append(path)


    def to_json(self):
        res = {}
        for name, dataset in self._dataset.items():
            res[name] = vars(dataset)
        return json.dumps({"samples": res})


    @abstractmethod
    def get_weight(self, process_name: str, int_lumi: float, *args, **kwargs) -> float:
        return 0.


    def get_samples(self) -> Generator[str, str, list[str]]:
        for name, sample in self._dataset.items():
            yield name, sample.trees[0], sample.files
