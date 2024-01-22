import os
import click
from ismn.interface import ISMN_Interface

@click.command("collect_metadata", short_help="Collect all ISMN metadata.")
@click.argument('data_path', type=click.STRING)
@click.option('--meta_path', type=click.Path(writable=True), default=None,
              help="Directory where the metadata should be stored. The file"
                   "will be created automatically. Existing metadata in this"
                   "directory will be replaced! If not specified, "
                   "we use DATA_PATH.")
@click.option('--parallel', '-p', is_flag=True, show_default=True,
              default=False,
              help="Pass this flag to activate parallel metadata collection "
                   "(recommended for large archives). Deactivated by default."
              )
def collect_metadata(data_path, meta_path, parallel):
    """
    Command line program to initialise ISMN metadata collection.
    THIS WILL OVERWRITE ANY EXISTING METADATA!

    \b
    DATA_PATH: string
        Path where the downloaded ISMN archive is stored. This is either
        - The downloaded ISMN ZIP archive or
        - A directory with network folders extracted from the ZIP archive.
        ISMN data can be downloaded from https://ismn.earth after registration.
    """
    # The docstring above is slightly different to the normal python one to
    # display it properly on the command line.
    if not os.path.exists(data_path):
        raise ValueError("The passed DATA_PATH does not exist.")
    if meta_path is not None:
        os.makedirs(meta_path, exist_ok=True)
    _ = ISMN_Interface(data_path, force_metadata_collection=True,
                       meta_path=meta_path, parallel=parallel)

@click.command("export_geojson", short_help="Export ISMN sensors to geojson.")
@click.argument('data_path', type=click.STRING)
@click.option('--file_out',
              type=click.STRING, default=None,
              help="Path to the json file that should be created. "
                   "If the file already exists it will be overwritten. "
                   "If not specified this is a file called "
                   "`ismn_sensors.json` and stored in the DATA_PATH.")
@click.option('--markercolor', '-m',
              type=click.STRING, default='"#00aa00"', show_default=True,
              help='Hex color (USE QUOTES!, e.g. "#00aa00") to assign to '
                   'markers in json file. The default color is green.')
def export_geojson(data_path, file_out, markercolor):
    """
    Calls
    Command line program to initialise ISMN metadata collection. THIS WILL
    OVERWRITE ANY EXISTING METADATA!

    \b
    Parameters
    ----------
    DATA_PATH: string
        Path where the downloaded ISMN archive is stored. This is either
        - The downloaded ISMN ZIP archive or
        - A directory with network folders extracted from the ZIP archive.
        ISMN data can be downloaded from https://ismn.earth after registration.
    """
    # The docstring above is slightly different to the normal python one to
    # display it properly on the command line.
    markercolor = str(markercolor.replace('"', '').replace("'", ""))
    if not os.path.exists(data_path):
        raise ValueError("The passed DATA_PATH does not exist.")
    ds = ISMN_Interface(data_path)
    if file_out is None:
        file_out = os.path.join(ds.root.root_dir, 'ismn_sensors.json')
    os.makedirs(os.path.dirname(file_out), exist_ok=True)
    print(f"Exporting geojson to: {file_out}")
    ds.collection.export_geojson(file_out, markercolor=markercolor)


@click.group(short_help="ISMN Command Line Programs.")
def ismn():
    pass

ismn.add_command(collect_metadata)
ismn.add_command(export_geojson)
