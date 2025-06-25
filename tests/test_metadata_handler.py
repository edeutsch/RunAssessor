import sys

sys.path.append("../lib")
from metadata_handler import MetadataHandler


def test_init_with_various_parameters():
    metadata = MetadataHandler()
    assert metadata.metadata_filepath.endswith('.json')
    metadata = MetadataHandler('')
    assert metadata.metadata_filepath.endswith('.json')
    metadata = MetadataHandler('.')
    assert metadata.metadata_filepath.endswith('.json')
    assert len(metadata.sdrf_hints) == 0


def test_create_template():
    metadata = MetadataHandler()
    assert(hasattr(metadata, 'metadata'))
    result = metadata.create_template()
    assert len(metadata.metadata) == 6
    assert len(metadata.metadata['problems']['errors']) == 3
    assert metadata.metadata['problems']['errors']['count'] == 0


