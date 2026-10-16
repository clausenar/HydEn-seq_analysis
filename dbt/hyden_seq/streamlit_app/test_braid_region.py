"""Unit tests for braid_region.py - no Snowflake or Streamlit needed.

Run from this directory:
    python -m pytest -q test_braid_region.py
"""

import numpy as np
import pandas as pd
import pytest

import braid_region as br


class TestRegionBounds:
    def test_window_is_centred_and_start_aligned_to_bins(self):
        assert br.region_bounds(10_250, 1_000, 500, 100_000) == (9_000, 11_250)

    def test_window_is_clipped_to_the_chromosome(self):
        assert br.region_bounds(300, 1_000, 100, 5_000) == (0, 1_300)
        assert br.region_bounds(4_900, 1_000, 100, 5_000) == (3_900, 5_000)

    def test_too_many_bins_is_refused(self):
        with pytest.raises(ValueError, match="too many"):
            br.region_bounds(500_000, 500_000, 1, 1_000_000)

    def test_limit_itself_is_allowed(self):
        start, end = br.region_bounds(10_000, 10_000, 1, 1_000_000)
        assert end - start == br.MAX_BINS

    @pytest.mark.parametrize("centre, half, bin_size", [(-1, 10, 1), (5_000, 10, 1), (10, 0, 1), (10, 10, 0)])
    def test_invalid_inputs_are_refused(self, centre, half, bin_size):
        with pytest.raises(ValueError):
            br.region_bounds(centre, half, bin_size, 5_000)


class TestQuery:
    def test_one_placeholder_per_parameter(self):
        samples = ["a", "b", "c"]
        sql = br.region_query(len(samples))
        params = br.region_params(samples, "chrI", 0, 1_000, 50)
        assert sql.count("?") == len(params)

    def test_values_are_bound_not_pasted_into_the_sql(self):
        sql = br.region_query(1)
        assert "chrI" not in sql and "'a'" not in sql


class TestDenseBins:
    def rows(self):
        return pd.DataFrame({
            "SAMPLE_NAME": ["s1", "s1", "s2"],
            "BIN_START": [100, 300, 200],
            "WATSON": [3.0, 0.0, 1.0],
            "CRICK": [1.0, 7.0, 0.0],
        })

    def test_every_bin_of_every_sample_is_present(self):
        dense = br.dense_bins(self.rows(), ["s1", "s2"], 100, 400, 100)
        assert len(dense) == 2 * 3
        assert dense.groupby("SAMPLE_NAME")["BIN_START"].apply(list).to_dict() == {
            "s1": [100, 200, 300], "s2": [100, 200, 300],
        }

    def test_empty_bins_count_as_zero_and_ratio_uses_pseudocount_1(self):
        dense = br.dense_bins(self.rows(), ["s1", "s2"], 100, 400, 100).set_index(["SAMPLE_NAME", "BIN_START"])
        assert dense.loc[("s1", 200), "LOG2_WATSON_CRICK"] == 0.0
        assert dense.loc[("s1", 100), "LOG2_WATSON_CRICK"] == pytest.approx(np.log2(4 / 2))
        assert dense.loc[("s1", 300), "LOG2_WATSON_CRICK"] == pytest.approx(np.log2(1 / 8))

    def test_no_signal_at_all_gives_a_flat_zero_line(self):
        empty = pd.DataFrame(columns=["SAMPLE_NAME", "BIN_START", "WATSON", "CRICK"])
        dense = br.dense_bins(empty, ["s1"], 0, 500, 100)
        assert len(dense) == 5 and (dense["LOG2_WATSON_CRICK"] == 0).all()


def test_plot_has_one_panel_per_sample():
    dense = br.dense_bins(TestDenseBins().rows(), ["s1", "s2"], 100, 400, 100)
    fig = br.plot_region(dense, ["s1", "s2"], "chrI", 100, 400, 100, origins=[250, 9_999])
    assert len(fig.axes) == 2
    assert fig.axes[0].get_xlim() == (100, 400)
