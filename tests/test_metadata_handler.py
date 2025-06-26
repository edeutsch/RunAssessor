import sys
import os
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

        file = metadata.find_txt_file()
        assert file == "study_metadata.txt"    

        #### study_metadata.txt file exists 
        metadata = MetadataHandler("study_metadata.json", verbose=2)
        file = metadata.find_txt_file()
        assert file == "study_metadata.txt" 

        #### defualt study_metadata.json is used so base template should be used
        patcher.fs.remove("study_metadata.txt")
        metadata = MetadataHandler(verbose=2)
        file = metadata.find_txt_file()
        
        #### study_metadata.txt is specified by user but does not exits, SDRF table wont be made
        metadata = MetadataHandler("study_metadata.json", verbose=2)
        file = metadata.find_txt_file()
        assert file == None

        #### Should not find baloon.txt and not generate the 
        metadata = MetadataHandler("balloon.json", verbose=2)
        file = metadata.find_txt_file()
        assert file == None