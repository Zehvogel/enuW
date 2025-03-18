from collections import namedtuple
from itertools import chain, filterfalse
import re
from typing import Any
from collections import defaultdict
import ROOT
import json
from overrides import override

Stats = namedtuple("Stats", ["int_lumi", "n_events"])


def id_genrange():
    # the generator overlords have chosen these numbers...
    higgs_inc_range = range(402001, 402014 + 1)
    # some are missing:
    missing = [500000 + i for i in [
        5, 7, 9, 11, 23, 24, 31, 32, 41, 42, 49, 50, 61, 63, 65, 67, 69,
        71, 73, 75, 77, 79, 81, 83, 85, 87, 89, 91, 93, 95, 97, 99, 109,
        111, 121, 123, 129, 133, 137
        ]]
    all2fto5f_range = filterfalse(lambda i: i in missing, range(500001, 500248 + 1))
    the6f_range     = range(402301, 402324 + 1)
    for id in chain(higgs_inc_range, all2fto5f_range, the6f_range):
        yield str(id)

# Petroff 10 color scheme backport
kP10 = [
ROOT.TColor( 63./255., 144./255., 218./255.),# "kP10Blue"),
ROOT.TColor(       1., 169./255.,  14./255.),# "kP10Yellow"),
ROOT.TColor(189./255.,  31./255.,   1./255.),# "kP10Red"),
ROOT.TColor(131./255.,  45./255., 182./255.),# "kP10Violet"),
ROOT.TColor(169./255., 107./255.,  89./255.),# "kP10Brown"),
ROOT.TColor(231./255.,  99./255.,        0.),# "kP10Orange"),
ROOT.TColor(185./255., 172./255., 112./255.),# "kP10Green"),
ROOT.TColor(113./255., 117./255., 129./255.),# "kP10Ash"),
ROOT.TColor(148./255., 164./255., 162./255.),# "kP10Gray"),
ROOT.TColor(146./255., 218./255., 221./255.),# "kP10Cyan")
]

def process_id(path: str) -> str:
    dir, _, fname = path.rpartition("/")
    # the mc-2020 filenames follow a certain fixed naming scheme, example:
    # rv02-02.sv02-02.mILD_l5_o1_v02.E250-SetA.I500102.P4f_sze_sl.eL.pR.n024.d_dstm_15180_122_mini-DST.edm4hep.root
    # <reco v>.<sim v>.<det model>.<machine setting>.<process id>.<process name>.<e pol>.<p pol>.<... other stuff
    parts = fname.split(".")
    return parts[4].lstrip("I")

class Dataset:
    """Initialise a dataset from a list of file paths in ILD mc-2020 production notation"""
    # we want to have something like {"process": {"pol: [path1, path2]"}}
    _dataset = defaultdict(lambda: defaultdict(list))

    # TODO I want to be able to supply either a list or a file :/
    def __init__(self, input_path: str | list[str]):
        def process_lines(content):
            for line in content:
                path = line.strip()
                self.add_file(path)
        if isinstance(input_path, str):
            with open(input_path) as file:
                process_lines(file)
        elif isinstance(input_path, list):
            process_lines(input_path)
        else:
            return TypeError


    def process_path(self, path: str) -> str:
        dir, _, fname = path.rpartition("/")
        # the mc-2020 filenames follow a certain fixed naming scheme, example:
        # rv02-02.sv02-02.mILD_l5_o1_v02.E250-SetA.I500102.P4f_sze_sl.eL.pR.n024.d_dstm_15180_122_mini-DST.edm4hep.root
        # <reco v>.<sim v>.<det model>.<machine setting>.<process id>.<process name>.<e pol>.<p pol>.<... other stuff
        parts = fname.split(".")
        process_name = parts[5].lstrip("P")
        e_pol = parts[6]
        p_pol = parts[7]
        return process_name, e_pol, p_pol

    def get_dataset(self):
        # convert back to regular dict
        return dict(self._dataset)

    # you probably do not want to call this yourself
    # I will probably put in stuff soon that assumes to only be called once the Dataset is complete
    # but no explicit checks...
    def add_file(self, path: str):
        process_name, e_pol, p_pol = self.process_path(path)
        self._dataset[process_name][e_pol + p_pol].append(path)

    def get_files(self, process_name: str, pol: str) -> list[str]:
        # use the getter to not have a defaultdict anymore and just let the dict throw
        return self.get_dataset()[process_name][pol]

    def get_keys(self):
        return [(process_name, pol) for process_name in self._dataset.keys() for pol in self._dataset[process_name].keys()]

    def get_merged_keys(self):
        return [f"{process_name}_{pol}" for process_name in self._dataset.keys() for pol in self._dataset[process_name].keys()]


class SnapshotDataset(Dataset):
    @override
    def process_path(self, path: str) -> str:
        # urgh here I relly made my life unnecessarily hard :(
        # my snapshot names look like this:
        # <category>_<process_name>_<pol>.snapshot.root
        # where the first two can also contain '_' and pol doesn't...
        # I can only split this under the assumption that the
        # process name starts with either '<n>f' or matches '<any>h'
        parts = path.split("_")
        process_name = parts[5].lstrip("P")
        e_pol = parts[6]
        p_pol = parts[7]
        return process_name, e_pol, p_pol




class Analysis:
    """Holds a Dataset, all the data frames and possibly also all the histograms"""

    # _df: dict[str, ROOT.RDataFrame]
    # _df: dict[str, cppyy.gbl.ROOT.RDataFrame]
    _df: dict[str, Any] = {}
    _stats: dict[str, Stats] = {}
    _counts = {}
    _histograms = {}
    _weights = {}
    _stacks = {}
    _canvases = {}
    # need a place to "park" them somewhere as the THStacks do not take ownership :(
    _scaled_histograms = {}
    _legends = {}
    _xsecs = {}
    _k2id = {}
    _reports = {}
    _categories: dict[str, list[str]] = {}
    _snapshots = {}

    def __init__(self, dataset: Dataset):
        self._dataset = dataset
        for process, pol in dataset.get_keys():
            files = dataset.get_files(process, pol)
            df = ROOT.RDataFrame("events", files)
            k = f"{process}_{pol}"
            self._df[k] = df
            # need process id to get genmeta like xsec
            id = process_id(files[0])
            self._k2id[k] = id


    def Define(self, *args):
        for k, df in self._df.items():
            self._df[k] = df.Define(*args)


    def populate_cross_sections(self, path: str):
        """Populates cross sections from ILD genmetaByID.json"""
        with open(path) as file:
            meta = json.load(file)
            for k, id in self._k2id.items():
                process_meta = meta[id]
                xsec = float(process_meta["cross_section_in_fb"])
                self._xsecs[k] = xsec


    def init_parameters(self, params: list[tuple[str, str, str]]):
        """Inits the podio generic parameters supplied as a list of (name, c++ typename, alternative value)"""
        ROOT.gInterpreter.Declare("#include <podio/GenericParameters.h>")
        self.Define("Parameters", "podio::GenericParameters par; par.loadFrom(GPDoubleKeys, GPDoubleValues); par.loadFrom(GPFloatKeys, GPFloatValues); par.loadFrom(GPIntKeys, GPIntValues); par.loadFrom(GPStringKeys, GPStringValues); return par;")
        for p_name, p_type, p_alt in params:
            self.Define(f"params_{p_name.replace('.', '_')}", f"Parameters.get<{p_type}>(\"{p_name}\").value_or({p_alt})")


    def book_statistics(self):
        # loop over all df, maybe also wrap this in the future
        for k, df in self._df.items():
            count = df.Count()
            # FIXME: refactor
            # self._stats[k] = Stats(int_lumi=0, n_events=count)
            self._counts[k] = count


    def print_stats(self):
        # for k, stats in self._stats.items():
            # n_events = stats.n_events.GetValue()
        for k, count in self._counts.items():
            xsec = self._xsecs[k]
            lumi = count.GetValue() / xsec
            print(f"process: {k} events: {count.GetValue()} int lumi: {lumi}")


    def evaluate_completeness(self, path: str):
        """Compare event counts in dataset to genmeta json file (generated amount of events)"""
        found = []
        with open(path) as file:
            meta = json.load(file)
            for k, id in self._k2id.items():
                found.append(id)
                process_meta = meta[id]
                total_generated = int(process_meta["total_number_of_events"])
                count = self._counts[k].GetValue()
                fraction_available = count / total_generated
                print(f"process: {k}, id: {id}, generated: {total_generated}, available: {count} ({fraction_available*100:.2f}%)")
            # print missing processes
            for id in id_genrange():
                if id not in found:
                    process_meta = meta[id]
                    name = process_meta["process_names"]
                    total_generated = int(process_meta["total_number_of_events"])
                    e_pol = process_meta["polarization1"]
                    p_pol = process_meta["polarization2"]
                    print(f"missing process {name}_e{e_pol}p{p_pol} with id {id} and {total_generated} generated events")


    def get_keys(self):
        return self._dataset.get_merged_keys()


    def book_histogram_1D(self, name: str, column: str, config: tuple):
        histograms = {}
        for k, df in self._df.items():
            histograms[k] = df.Histo1D(config, column)
        self._histograms[name] = histograms


    def _calculate_lumi_weight(self, k: str, int_lumi: float):
        # return int_lumi / (self._stats[k].n_events.GetValue() / self._xsecs[k])
        return int_lumi / (self._counts[k].GetValue() / self._xsecs[k])


    def _calculate_pol_weight(self, k: str, e_pol: float, p_pol: float):
        process_e_pol = -1.0 if "eL" in k else 1.0 if "eR" in k else 0.0
        process_p_pol = -1.0 if "pL" in k else 1.0 if "pR" in k else 0.0

        return 0.25 * (1.0 + e_pol * process_e_pol) * (1.0 + p_pol * process_p_pol)


    def _calculate_weight(self, k: str, int_lumi: float, e_pol: float, p_pol: float):
        lumi_weight = self._calculate_lumi_weight(k, int_lumi)
        pol_weight = self._calculate_pol_weight(k, e_pol, p_pol)
        return lumi_weight * pol_weight


    def get_weights(self, int_lumi: float, e_pol: float, p_pol: float):
            return {k: self._calculate_weight(k, int_lumi, e_pol, p_pol) for k in self._df}


    def run(self):
        """Execute all booked computations"""
        objects = list(self._counts.values())
        objects += list(self._snapshots.values())
        objects += list(self._reports.values())
        # unpacking the histograms is more involved
        for process in self._histograms.values():
            objects += process.values()
        ROOT.RDF.RunGraphs(objects)


    def draw_histogram(self, name: str, int_lumi: float = 5000, e_pol: float = 0.0, p_pol: float = 0.0):
        histograms = self._histograms[name]
        stack = ROOT.THStack()
        params = (name, int_lumi, e_pol, p_pol)
        self._scaled_histograms[params] = {}
        legend = ROOT.TLegend(0.6, 0.7, 1., 1,)
        # FIXME draw by category and use one color for each, just in the order of definition
        for i, (category_name, dataframes) in enumerate(self._categories.items()):
            scaled_cat_histograms = []
            if len(dataframes) == 0:
                # nothing to do for empty categories
                continue
            for k in dataframes:
                # get histogram, clone it, scale it, put it in a list
                h = histograms[k].Clone()
                # our crossection table only knows about the full processes
                weight_key = k.removesuffix("_signal").removesuffix("_bkg")
                # print(k, weight_key)
                weight = self._calculate_weight(weight_key, int_lumi, e_pol, p_pol)
                h.Scale(weight)
                scaled_cat_histograms.append(h)
            # sum up scaled category histograms, starting with the first one
            h = sum(scaled_cat_histograms[1:], start=scaled_cat_histograms[0])
            h.SetFillColor(kP10[i].GetNumber())
            # TODO: parse k to nice process name
            legend.AddEntry(h, category_name, "f")
            stack.Add(h)
            # store h so that it does not get deleted
            self._scaled_histograms[params][category_name] = h
        legend.SetNColumns(2)
        self._legends[params] = legend
        self._stacks[params] = stack
        canvas = ROOT.TCanvas()
        self._canvases[params] = canvas
        stack.SetTitle(f";{name}")
        stack.Draw("hist")
        legend.Draw()
        canvas.Draw()


    def add_filter(self, expression: str, name: str):
        for k, df in self._df.items():
            self._df[k] = df.Filter(expression, name)


    def book_reports(self):
        for k, df in self._df.items():
            self._reports[k] = df.Report()


    # TODO: adapt to use the new categories
    def print_reports(self, int_lumi: float = 5000, e_pol: float = 0.0, p_pol: float = 0.0):
        """Print cut-flow report for each category/prefix"""
        # all processes should have the same filters applied:
        # should still work with categories because the cut frames are defined last
        names = list(list(self._df.values())[0].GetFilterNames())
        n_filters = len(names)
        # aggregate all the information
        largest = 0
        numbers = {}
        errors2 = {}
        for category_name, frames in self._categories.items():
            nums = [0.0] * (n_filters + 1)
            errs2 = [0.0] * (n_filters + 1)
            # print(frames)
            for k in frames:
                reports = self._reports[k]
                weight_key = k
                if category_name.endswith("_signal") or category_name.endswith("_bkg"):
                    weight_key  = weight_key.removesuffix("_signal").removesuffix("_bkg")
                    # need to remove signal/background definition cut to get meaningful numbers and correct alignment of samples
                    # XXX: NO! signal cut is unnamed so not in the report!
                    # from itertools import islice
                    # reports = islice(reports, 1, None)
                weight = self._calculate_weight(weight_key, int_lumi, e_pol, p_pol)
                for i, cut_info in enumerate(reports):
                    # print(f"processing cutinfo for category: {category_name}, frame: {k}, cut_name: {cut_info.GetName()}, pass: {cut_info.GetPass()}")
                    if i == 0:
                        n = cut_info.GetAll() * weight
                        nums[0] += n
                        errs2[0] += n * weight # need to square weight
                    n = cut_info.GetPass() * weight
                    nums[i+1] += n
                    errs2[i+1] += n * weight # need to square weight
            # print(nums)
            numbers[category_name] = nums
            errors2[category_name] = errs2
            largest = max(largest, max(nums))

        # do the actual printing
        # TODO: just also give a fixed format for the yields like Teresa
        from math import log10, ceil, sqrt
        max_digits = ceil(log10(largest+1)) #+ len(" 0e-00")
        max_name = max((len(name) for name in self._categories))
        offset = max(max_digits, max_name) + 1
        line = ""
        for category_name in self._categories:
            line_piece = f"{category_name: >{offset + len(' 0e-00')}}"
            line += line_piece
        print(line)
        for i, name in enumerate(["All"] + names):
            line = ""
            for category_name in self._categories:
                num = numbers[category_name][i]
                err = sqrt(errors2[category_name][i])
                line_piece = f"{round(num): >{offset}} ({err:.0e})"
                # print(f"::{line_piece}::")
                line += line_piece
            print(line, name)
        line = ""
        for category_name in self._categories:
            eff = numbers[category_name][-1] / numbers[category_name][0] if numbers[category_name][0] > 0 else 0
            if eff >= 0.01:
                line_piece = f"{eff: >{offset}.2f}"
            else:
                line_piece = f"{eff: >{offset}.0e}"
            line += line_piece
        print(line, "efficiency")


    def set_categories(self, input: dict[str, dict[str, str]]):
        """Defines categories that can be used that group multiple samples together or just use parts of a sample. The categories are expected to be disjoint. Two signal categories can not share the same processes."""
        # example_input = {
        # "4f_sw_sl_signal": {"pattern": "4f_sw_sl", "cut": "pre_defined_signal_cut"},
        # "4f_sl_bkg": {"pattern": r"4f\w+sl", "cut": None}, # inverse signal cut is applied automatically where necessary
        # "3f": {"pattern": "ea_3f|ae_3f", "cut": None} # FIXME: cutting not even supported at the moment, add a signal category if needed :(
        #}
        # TODO: The second part could also just be a named tuple
        # TODO: at the moment signals are supposed to end in '_signal' and backgrounds that need the inverse signal cut need to end in '_bkg' but cleaner would be to set a bool
        categories: dict[str, list[str]] = {}
        signals: list[str] = []
        for category_name, category_def in input.items():
            dataframes = []
            pattern = re.compile(category_def["pattern"])
            for df_name in filter(pattern.search, self._df.keys()):
                # print(category_name, df_name)
                dataframes.append(df_name)
            if category_name.endswith("_signal"):
                signals.append(category_name)
            categories[category_name] = dataframes
        # now apply cuts and clean up the backgrounds
        for signal_name in signals:
            signal_def = input[signal_name]
            # need to apply cut to each dataframe of the signal
            cut = signal_def.get("cut")
            new_df_names = []

            # maybe a bit ugly, sorry :(
            def rename_and_cut(df_name, background=False):
                new_df = self._df[df_name]
                if cut:
                    if not background:
                        new_df = new_df.Filter(f"return {cut};")
                    else:
                        new_df = new_df.Filter(f"return !({cut});")
                new_df_name = f"{df_name}_signal" if not background else f"{df_name}_bkg"
                # print(f"creating {new_df_name} from {df_name} with background == {background}")
                self._df[new_df_name] = new_df
                return new_df_name


            for df_name in categories[signal_name]:
                new_df_names.append(rename_and_cut(df_name))
            # overwrite old df name list
            categories[signal_name] = new_df_names

            # now look for any backgrounds to rename and apply inverse cut
            pattern = re.compile(signal_def["pattern"])
            for category_name, dataframes in categories.items():
                if not category_name.endswith("_bkg"):
                    continue
                new_df_names = dataframes.copy()
                for df_name in filter(pattern.search, dataframes):
                    # print(f"signal {signal_name} has background {df_name}")
                    new_df_names.remove(df_name)
                    new_df_names.append(rename_and_cut(df_name, True))
                # overwrite old df name list
                categories[category_name] = new_df_names

        self._categories = categories


    def is_complete_categorisation(self):
        df_names = list(self._df.keys())
        categories = self._categories
        for c, frames in categories.items():
            for frame in frames:
                short_frame = frame.removesuffix("_bkg").removesuffix("_signal")
                if short_frame in df_names:
                    df_names.remove(short_frame)
                if frame in df_names:
                    df_names.remove(frame)


        return len(df_names) == 0

        # TODO: Add a method to store selected histograms to a root file

    # FIXME: need to also write-out some metadata, especially the initial counts per process!
    def book_snapshots(self, tree_name: str, out_dir: str, column_list = None, no_rvec = False):
        """Will be written out when the event loop runs"""
        # Need to avoid double counting, either write out all non-suffixed dataframes or start from the categories
        categories = self._categories
        for category_name, frames in categories.items():
            for frame in frames:
                # # we only want to write out non-empty frames because RDataFrame always creates the file
                # # and in the multi-threaded case it does not even put the empty tree inside :(
                # # We need counts after the last cut, and we don't want to recalculate them again
                # # But we can take them from the cut reports
                # report = self._reports[frame]
                # # get the last cut
                # from collections import deque
                # last_cut = deque(report, maxlen=1).pop()
                # if not last_cut.GetPass() > 0:
                #     continue
                # if we arrive here there is actually something to do
                df = self._df[frame]
                file_name = f"{out_dir}/{category_name}_{frame}.snapshot.root"
                args = [tree_name, file_name]
                if column_list:
                    args.append(column_list)
                else:
                    # append empty column selection regex -> all columns
                    selection = ""
                    args.append(selection)
                snapshot_options = ROOT.RDF.RSnapshotOptions()
                snapshot_options.fLazy = True
                snapshot_options.fCompressionLevel = 9
                snapshot_options.fCompressionAlgorithm = ROOT.RCompressionSetting.EAlgorithm.kZSTD
                if no_rvec:
                    snapshot_options.fVector2RVec = False
                args.append(snapshot_options)
                self._snapshots[frame] = df.Snapshot(*args)
