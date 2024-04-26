#!/usr/bin/env python
# coding: utf-8

import ROOT

df = ROOT.RDataFrame(10)

code = """
#include <array>
#include <vector>
struct test {
  std::array<float, 6> values{};
};

ROOT::VecOps::RVec<test> makeTestData() {
  return ROOT::VecOps::RVec<test>(10);
}
"""
ROOT.gInterpreter.Declare(code)

df = df.Define("test_col", "makeTestData()")
print(df.Describe())
df.Snapshot("test_tree", "test_file.root")
df2 = ROOT.RDataFrame("test_tree", "test_file.root")
print(df2.Describe())
