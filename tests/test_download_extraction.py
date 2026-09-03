import os
import tempfile
import yaml
import shutil
import datetime
from ismn.download import ISMNExtractor

testdata_root = os.path.join(os.path.dirname(__file__), "test_data")

def test_download_extraction():
    """
    Test the extraction process of the ISMN dataset.
    """

    with tempfile.TemporaryDirectory() as tmpdir:
        archive_path = os.path.join(testdata_root, 'archive_sample', "ismn_archive.zip")

        output_dir = os.path.join(tmpdir, "ISMN")
        nrt_networks = ['RSMN']

        # Pre-populate the output dir with the RSMN folder
        os.makedirs(output_dir, exist_ok=True)
        rsmn_src = os.path.join(testdata_root, 'archive_sample', "RSMN")
        rsmn_dst = os.path.join(output_dir, "RSMN")
        shutil.copytree(rsmn_src, rsmn_dst)

        extractor = ISMNExtractor(
            archive_path=archive_path,
            output_dir=output_dir,
            nrt_networks=nrt_networks
        )
        extractor.run()

        # The overview file should be generated in the output dir
        overview_path = os.path.join(output_dir, "overview.yml")
        assert os.path.isfile(overview_path), f"{overview_path} was not created"

        # And its contents should match the expected metadata
        with open(overview_path) as f:
            overview = yaml.safe_load(f)

        assert overview == {
            "product": "ISMN",
            "version": "v202505",
            "period_to": datetime.date(2024, 12, 31),
        }