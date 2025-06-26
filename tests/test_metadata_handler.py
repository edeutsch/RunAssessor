import sys

sys.path.append("../lib")
from metadata_handler import MetadataHandler
from pyfakefs.fake_filesystem_unittest import Patcher

def test_init_with_various_parameters():
    metadata = MetadataHandler()
    assert metadata.metadata_filepath.endswith('.json')
    metadata = MetadataHandler('')
    assert metadata.metadata_filepath.endswith('.json')
    metadata = MetadataHandler('.')
    assert metadata.metadata_filepath.endswith('.json')

    metadata = MetadataHandler('.jso') #Should be able to make a new metadata file path if not relevant
    assert metadata.metadata_filepath.endswith('.json')

    metadata= MetadataHandler('study_metadata')
    assert metadata.metadata_filepath.endswith('.json')

    assert len(metadata.sdrf_hints) == 0


def test_create_template():
    metadata = MetadataHandler()
    assert(hasattr(metadata, 'metadata'))
    result = metadata.create_template()
    assert len(metadata.metadata) == 6
    assert len(metadata.metadata['problems']['errors']) == 3
    assert metadata.metadata['problems']['errors']['count'] == 0


#Unit test to check that function will use template file if study_metadata.txt is not availible
def test_find_txt_file(capsys):
    with Patcher() as patcher:

        #### study_metadata.txt file exists
        patcher.fs.create_file("study_metadata.txt")
        metadata = MetadataHandler("study_metadata.json", verbose=2)

        metadata.find_txt_file()
        captured = capsys.readouterr()
        assert "study_metadata.txt found" in captured.out or "study_metadata.txt found" in captured.err      

        #### study_metadata.txt file exists 
        metadata = MetadataHandler("study_metadata.json", verbose=2)
        metadata.find_txt_file()
        captured = capsys.readouterr()
        assert "study_metadata.txt found" in captured.out or "study_metadata.txt found" in captured.err   

        #### study_metadata.txt does not exist so template should be used
        patcher.fs.remove("study_metadata.txt")
        metadata = MetadataHandler("study_metadata.json", verbose=2)
        metadata.find_txt_file()
        captured = capsys.readouterr()
        assert "/study_metadata_template.txt" in captured.out or "/study_metadata_template.txt" in captured.err

        #### Should try and use study_metadata.ballon.txt (DNE) and use template
        metadata = MetadataHandler("study_metadata.ballon", verbose=2)
        metadata.find_txt_file()
        captured = capsys.readouterr()
        assert "/study_metadata_template.txt" in captured.out or "/study_metadata_template.txt" in captured.err
