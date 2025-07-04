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


#Unit test to check that function will use template file if study_metadata.txt is not availible
def test_find_txt_file(capsys):

    ### Creates study_metadata.txt
    filename = "study_metadata.txt"
    file_path = os.path.join(os.getcwd(), filename)
    open(file_path, 'w').close()

    #### study_metadata.txt file exists
    metadata = MetadataHandler("study_metadata.json", verbose=2)
    file = metadata.find_txt_file()
    assert file == "study_metadata.txt"    

    #### study_metadata.txt file exists 
    metadata = MetadataHandler("study_metadata.json", verbose=2)
    file = metadata.find_txt_file()
    assert file == "study_metadata.txt" 
    
    #### defualt study_metadata.json is used and no stud_metadat.txt ecists so base template should be used
    #### study_metadata.txt should be generated
    os.remove(file_path)
    metadata = MetadataHandler(verbose=2)
    file = metadata.find_txt_file()

    if os.path.exists(os.getcwd() + "/study_metadata.txt"):
        os.remove(os.getcwd() + "/study_metadata.txt")
    else:   
        raise AssertionError(f"File {file_path} was not generated:" + 
                             os.getcwd() + "/study_metadata.txt")
    
    #### study_metadata.txt is specified by user but does not exits, SDRF table wont be made
    metadata = MetadataHandler("study_metadata.json", verbose=2)
    file = metadata.find_txt_file()
    assert file == None

    #### Should not find baloon.txt and not generate the sdrf table
    metadata = MetadataHandler("balloon.json", verbose=2)
    file = metadata.find_txt_file()
    assert file == None