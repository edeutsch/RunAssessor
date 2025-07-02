import sys

sys.path.append("../lib")
from mzML_assessor import MzMLAssessor


def test_parse_filter_string_ms1():
    assessor = MzMLAssessor(mzml_file="fake_test_file.mzML.gz")
    stats = assessor.get_empty_stats()

    filter_string = 'FTMS + p NSI Full ms [350.00-1500.00]'
    result = assessor.parse_filter_string(filter_string,stats)
    assert result == None
    assert stats['n_ms1_spectra'] == 1
    assert stats['high_accuracy_precursors'] == 'true'
    assert stats['fragmentation_tag'] == '??'


def test_parse_filter_string_1():
    assessor = MzMLAssessor(mzml_file="fake_test_file.mzML.gz")
    stats = assessor.get_empty_stats()

    filter_string = 'FTMS + c NSI d Full ms2 697.37@cid32.00 [180.00-1405.00]'
    result = assessor.parse_filter_string(filter_string,stats)
    assert result == None
    assert stats['n_ms2_spectra'] == 1
    assert stats['n_HR_IT_CID_spectra'] == 1
    assert stats['fragmentation_type'] == 'HR_IT_CID'
    assert stats['fragmentation_tag'] == 'HR IT CID'


def test_parse_filter_string_2():
    assessor = MzMLAssessor(mzml_file="fake_test_file.mzML.gz")
    stats = assessor.get_empty_stats()

    filter_string = 'FTMS + c NSI d Full ms2 697.37@etd32.00 [180.00-1405.00]'
    result = assessor.parse_filter_string(filter_string,stats)
    assert result == None
    assert stats['n_ms2_spectra'] == 1
    assert stats['n_HR_IT_ETD_spectra'] == 1
    assert stats['fragmentation_type'] == 'HR_IT_ETD'
    assert stats['fragmentation_tag'] == 'HR IT ETD'


def test_parse_filter_string_3():
    assessor = MzMLAssessor(mzml_file="fake_test_file.mzML.gz")
    stats = assessor.get_empty_stats()

    filter_string = 'FTMS + c NSI d Full ms2 697.37@hcd32.00 [180.00-1405.00]'
    result = assessor.parse_filter_string(filter_string,stats)
    assert result == None
    assert stats['n_ms2_spectra'] == 1
    assert stats["n_HR_HCD_spectra"] == 1
    assert stats['fragmentation_type'] == 'HR_HCD'
    assert stats['fragmentation_tag'] == 'HR HCD'


def test_parse_filter_string_4():
    assessor = MzMLAssessor(mzml_file="fake_test_file.mzML.gz")
    stats = assessor.get_empty_stats()

    filter_string = ''
    result = assessor.parse_filter_string(filter_string,stats)
    assert result == None
    assert stats['high_accuracy_precursors'] == 'unknown'
    assert stats['fragmentation_type'] == 'unknown_fragmentation_type'
    assert stats['fragmentation_tag'] == 'unknown fragmentation type'


def test_parse_filter_string_5():
    assessor = MzMLAssessor(mzml_file="fake_test_file.mzML.gz")
    stats = assessor.get_empty_stats()

    filter_string = 'ITMS + c NSI d w Full ms2 857.62@cid30.00 [225.00-870.00]'
    result = assessor.parse_filter_string(filter_string,stats)
    assert result == None
    assert stats['n_ms2_spectra'] == 1
    assert stats['n_LR_IT_CID_spectra'] == 1
    assert stats['fragmentation_type'] == 'LR_IT_CID'
    assert stats['fragmentation_tag'] == 'LR IT CID'


def test_parse_filter_string_6():
    assessor = MzMLAssessor(mzml_file="fake_test_file.mzML.gz")
    stats = assessor.get_empty_stats()

    filter_string = 'ITMS + c NSI r d sa Full ms2 470.93@etd33.33 [50.00-1425.00]'
    result = assessor.parse_filter_string(filter_string,stats)
    assert result == None
    assert stats['n_ms2_spectra'] == 1
    assert stats['n_LR_IT_ETD_spectra'] == 1
    assert stats['fragmentation_type'] == 'LR_IT_ETD'
    assert stats['fragmentation_tag'] == 'LR IT ETD'


def test_parse_filter_string_7():
    assessor = MzMLAssessor(mzml_file="fake_test_file.mzML.gz")
    stats = assessor.get_empty_stats()
    
    filter_string = 'ITMS + c NSI r d Full ms2 509.2148@hcd30.00 [110.0000-1029.0000]'
    result = assessor.parse_filter_string(filter_string,stats)
    assert result == None
    assert stats['n_ms2_spectra'] == 1
    assert stats['n_LR_HCD_spectra'] == 1
    assert stats['fragmentation_type'] == 'LR_HCD'
    assert stats['fragmentation_tag'] == 'LR HCD'


def test_parse_filter_string_8():
    assessor = MzMLAssessor(mzml_file="fake_test_file.mzML.gz")
    stats = assessor.get_empty_stats()

    filter_string = 'FTMS + c NSI d sa Full ms2 557.2769@etd116.52@cid30.00 [120.0000-1125.0000]'
    result = assessor.parse_filter_string(filter_string,stats)
    assert result == None
    assert stats['n_ms2_spectra'] == 1
    assert stats['n_HR_ETciD_spectra'] == 1
    assert stats['fragmentation_type'] == 'HR_ETciD'
    assert stats['fragmentation_tag'] == 'HR ETciD'


def test_parse_filter_string_9():
    assessor = MzMLAssessor(mzml_file="fake_test_file.mzML.gz")
    stats = assessor.get_empty_stats()

    filter_string = 'ITMS + c NSI r d sa Full ms2 984.7690@etd95.59@hcd25.00 [110.0000-1980.0000]'
    result = assessor.parse_filter_string(filter_string,stats)
    assert result == None
    assert stats['n_ms2_spectra'] == 1
    assert stats['n_LR_EThcD_spectra'] == 1     # not ETciD
    assert stats['fragmentation_type'] == 'LR_EThcD'
    assert stats['fragmentation_tag'] == 'LR EThcD'


def test_parse_filter_string_10():
    assessor = MzMLAssessor(mzml_file = "fake_test_file.mzML.gz")
    stats = assessor.get_empty_stats()

    filter_string = 'FTMS + p NSI d sa Full ms2 982.3831@etd95.59@hcd25.00 [110.0000-1975.0000]'
    result = assessor.parse_filter_string(filter_string,stats)
    assert result == None
    assert stats['n_ms2_spectra'] == 1
    assert stats['n_HR_EThcD_spectra'] == 1
    assert stats['fragmentation_type'] == 'HR_EThcD'
    assert stats['fragmentation_tag'] == 'HR EThcD'