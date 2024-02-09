import os
from click.testing import CliRunner
from ismn.cli import collect_metadata, export_geojson
from tempfile import TemporaryDirectory

testdata_root = os.path.join(os.path.dirname(__file__), "test_data")

def test_cli_meta_collect():
    with TemporaryDirectory() as tempdir:
        data_path = os.path.join(
            testdata_root, "zip_archives", "ceop",
            "Data_seperate_files_20170810_20180809.zip")
        runner = CliRunner()
        result = runner.invoke(collect_metadata,
                               [data_path, "--meta_path", tempdir, "-p"])
        assert result.exit_code == 0
        assert os.path.isfile(os.path.join(
            tempdir, "Data_seperate_files_20170810_20180809.csv"))

def test_cli_export_geojson():
    with TemporaryDirectory() as tempdir:
        data_path = os.path.join(
            testdata_root, "zip_archives", "ceop",
            "Data_seperate_files_20170810_20180809.zip")
        runner = CliRunner()
        result = runner.invoke(export_geojson,
                               [data_path, "--file_out",
                                os.path.join(tempdir, "test.geojson"),
                                "-f", "network",
                                '-f', "timerange",
                                "-f", "lc_2010",
                                "-f", "variable",
                                "-var", "soil_moisture"])
        assert result.exit_code == 0
        assert os.path.isfile(os.path.join(tempdir, "test.geojson"))
        with open(os.path.join(tempdir, "test.geojson"), "r") as f:
            content = f.readlines()
            assert "network" in content[0]
            assert "timerange_from" in content[0]
            assert "timerange_to" in content[0]
            assert "lc_2010" in content[0]
            assert 'soil_moisture' in content[0]
            assert "precipitation" not in content[0]
