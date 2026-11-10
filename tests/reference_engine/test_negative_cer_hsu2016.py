from __future__ import annotations

import unittest

from sphingolipid_toolkit._vendor.lipidgate.ms2.models import (
    ExperimentalPeak,
    FragmentMatch,
    FragmentRecord,
    LibraryRecord,
)
from sphingolipid_toolkit._vendor.lipidgate.ms2.negative_cer_hsu2016 import (
    logical_negative_cer_fragment_types,
    normalize_negative_cer_hsu2016_fragments,
    selected_family,
)
from sphingolipid_toolkit._vendor.lipidgate.ms2.scoring import _with_negative_cer_match_bonus


def fragment(mz, name, role):
    return FragmentRecord(mz, name, role)


class NegativeCerHsu2016Tests(unittest.TestCase):
    def normalize(self, name, peaks, adduct='[M-H]-'):
        return normalize_negative_cer_hsu2016_fragments(name, 999., adduct, peaks)

    def test_only_selected_alpha_hydroxy_and_nonhydroxy_families(self):
        expected = {
            'Cer(d18:0/24:0)': 'd0_nfa',
            'Cer(d18:0/h16:0)': 'd0_alpha_hfa',
            'Cer(t18:0/20:0)': 't0_nfa',
        }
        for name, family in expected.items():
            for adduct in ('[M-H]-', '[M+HCOO]-', '[M+CH3COO]-'):
                self.assertEqual(selected_family(name, adduct), family)
        for name in ('Cer(d18:1/24:1)', 'Cer(d18:1/h16:0)',
                     'Cer(t18:0/h20:0)', 'Cerd18:0/24:0(O-16:0)'):
            self.assertIsNone(selected_family(name, '[M-H]-'))

    def test_dihydrosphinganine_nonhydroxy_removes_unsupported_routes(self):
        peaks = [
            fragment(650.6457, '[M-H]-', 'Precursor Ion'),
            fragment(620.6351, 'M-H-HCHO', 'C类碎片'),
            fragment(270.2802, 'LCB-H-HCHO', 'LCB碎片'),
            fragment(265.2537, 'M-H-H2O-RCONH', 'LCB碎片'),
            fragment(618.6195, 'M-H-H2-HCHO', 'C类碎片'),
            fragment(300.2908, 'LCB-H', 'LCB碎片'),
        ]
        result = self.normalize('Cer(d18:0/24:0)', peaks)
        self.assertEqual({f.name for f in result}, {'[M-H]-', 'M-H-H2-HCHO', 'LCB-H'})

    def test_dihydrosphinganine_alpha_hydroxy_adds_c6_and_b8(self):
        peaks = [
            fragment(554.5154, '[M-H]-', 'Precursor Ion'),
            fragment(270.2802, 'LCB-H-HCHO', 'LCB碎片'),
            fragment(300.2908, 'LCB-H', 'LCB碎片'),
        ]
        result = self.normalize('Cer(d18:0/h16:0)', peaks)
        by_name = {f.name: f for f in result}
        self.assertNotIn('LCB-H-HCHO', by_name)
        self.assertEqual(by_name['M-H-2H2O-HCHO'].mz, 488.4837)
        self.assertEqual(by_name['M-H-2H2O-HCHO'].fragment_type, 'C类碎片')
        self.assertEqual(by_name['LCB-H+CO'].mz, 328.2857)
        self.assertEqual(by_name['LCB-H+CO'].fragment_type, 'LCB碎片')

    def test_phytosphingosine_nonhydroxy_adds_a10_and_b3(self):
        peaks = [
            fragment(610.5780, '[M-H]-', 'Precursor Ion'),
            fragment(310.3115, '[RCONH]-(20:0)', 'FA类碎片'),
            fragment(316.2857, 'LCB-H', 'LCB碎片'),
        ]
        result = self.normalize('Cer(t18:0/20:0)', peaks)
        by_name = {f.name: f for f in result}
        self.assertEqual(by_name['[RCONH+C3H4O]-(20:0)'].mz, 366.3377)
        self.assertEqual(by_name['[RCONH+C3H4O]-(20:0)'].fragment_type, 'FA类碎片')
        self.assertEqual(by_name['LCB-H-H2-HCHO'].mz, 284.2595)
        self.assertEqual(by_name['LCB-H-H2-HCHO'].fragment_type, 'LCB碎片')

    def test_normalization_is_idempotent(self):
        cases = [
            ('Cer(d18:0/24:0)', [
                fragment(650.6457, '[M-H]-', 'Precursor Ion'),
                fragment(620.6351, 'M-H-HCHO', 'C类碎片'),
                fragment(270.2802, 'LCB-H-HCHO', 'LCB碎片'),
                fragment(265.2537, 'M-H-H2O-RCONH', 'LCB碎片'),
            ]),
            ('Cer(d18:0/h16:0)', [
                fragment(554.5154, '[M-H]-', 'Precursor Ion'),
                fragment(270.2802, 'LCB-H-HCHO', 'LCB碎片'),
                fragment(300.2908, 'LCB-H', 'LCB碎片'),
            ]),
            ('Cer(t18:0/20:0)', [
                fragment(610.5780, '[M-H]-', 'Precursor Ion'),
                fragment(310.3115, '[RCONH]-(20:0)', 'FA类碎片'),
                fragment(316.2857, 'LCB-H', 'LCB碎片'),
            ]),
        ]
        for name, peaks in cases:
            once = self.normalize(name, peaks)
            self.assertEqual(self.normalize(name, once), once)

    def test_exact_isobar_is_stored_once_with_two_logical_types(self):
        peaks = [
            fragment(372.3119, '[M-H]-', 'Precursor Ion'),
            fragment(172.1707, 'LCB-H-HCHO', 'LCB碎片'),
            fragment(202.1813, 'LCB-H', 'LCB碎片'),
            fragment(230.1762, '[NAE-H]-(10:0;O)', 'NAE碎片'),
        ]
        result = self.normalize('Cer(d11:0/h10:0)', peaks)
        at_mass = [f for f in result if f.mz == 230.1762]
        self.assertEqual(len(at_mass), 1)
        self.assertEqual(at_mass[0].name, '[NAE-H]-(10:0;O) | LCB-H+CO')
        self.assertEqual(logical_negative_cer_fragment_types(at_mass[0]),
                         {'NAE碎片', 'LCB碎片'})
        self.assertEqual(self.normalize('Cer(d11:0/h10:0)', result), result)

    def test_d18_1_16_0_two_isobars_are_each_stored_once(self):
        peaks = [
            fragment(237.2224, 'LCB-H-C2H7NO', 'LCB碎片'),
            fragment(237.2224, '[RCOO-H2O]-(16:0)', 'FA类碎片'),
            fragment(298.2752, 'LCB-H', 'LCB碎片'),
            fragment(298.2752, '[NAE-H]-(16:0)', 'NAE碎片'),
        ]
        result = self.normalize('Cer(d18:1/16:0)', peaks)
        self.assertEqual(len(result), 2)
        self.assertEqual(result[0].name, 'LCB-H-C2H7NO | [RCOO-H2O]-(16:0)')
        self.assertEqual(result[1].name, 'LCB-H | [NAE-H]-(16:0)')
        self.assertEqual(logical_negative_cer_fragment_types(result[0]),
                         {'LCB碎片', 'FA类碎片'})
        self.assertEqual(logical_negative_cer_fragment_types(result[1]),
                         {'LCB碎片', 'NAE碎片'})

    def test_negative_cer_ten_physical_match_bonus_caps_at_100(self):
        record = LibraryRecord(1, 'Cer', 'Cer(d34:1)', 'Cer(d18:1/16:0)',
                               536.5, '[M-H]-')
        matches = [
            FragmentMatch(fragment(float(i), str(i), 'Common'),
                          ExperimentalPeak(float(i), 1., .1), 0.)
            for i in range(10)
        ]
        self.assertEqual(_with_negative_cer_match_bonus(record, matches[:9], 70.), 70.)
        self.assertEqual(_with_negative_cer_match_bonus(record, matches, 70.), 80.)
        self.assertEqual(_with_negative_cer_match_bonus(record, matches, 95.), 100.)
        record.adduct = '[M+H]+'
        self.assertEqual(_with_negative_cer_match_bonus(record, matches, 70.), 70.)


if __name__ == '__main__':
    unittest.main()
