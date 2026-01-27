import sys
import os
sys.path.append("../lib")
from metadata_handler import MetadataHandler


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


# Unit tests to check that study_metadata.txt will be created based on template if not available
def test_find_txt_file():

    filename = "study_metadata.txt"
    file_path = os.path.join(os.getcwd(), filename)

    ### Create an empty study_metadata.txt
    open(file_path, 'w').close()

    #### Verify that function finds it
    metadata = MetadataHandler("study_metadata.json", verbose=2)
    file = metadata.find_txt_file()
    assert file == "study_metadata.txt"    

    #### Test that a study_metadata.txt is generated from the template
    os.remove(file_path)
    metadata = MetadataHandler(verbose=2)
    file = metadata.find_txt_file()

    if os.path.exists(file_path):
        os.remove(file_path)
    else:   
        raise AssertionError(f"File {file_path} was not generated")
    
    #### Verify that a balloon.txt is created to accompany the balloon.json file
    filename = "balloon.txt"
    file_path = os.path.join(os.getcwd(), filename)
    metadata = MetadataHandler("balloon.json", verbose=2)
    file = metadata.find_txt_file()
    assert file == filename
    os.remove(file_path)

