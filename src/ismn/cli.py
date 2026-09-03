import os
import click
from ismn.interface import ISMN_Interface
from ismn.download import ISMNDownloader, ISMNExtractor
from ismn.const import nrt_networks

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
@click.option('--file_out', '-o',
              type=click.STRING, default=None,
              help="Path to the json file that should be created. "
                   "If the file already exists it will be overwritten. "
                   "If not specified this is a file called "
                   "`ismn_sensors.json` and stored in the DATA_PATH.")
@click.option('--field', '-f', multiple=True,
              help="Fields to include. This option can be called multiple times"
                   "with different fields. Allowed are: "
                   "network, station, sensor, depth, timerange. \n "
                   "Or any sensor properties (also custom ones) that have a "
                   "value.")
@click.option('--variable', '-var', multiple=True,
              help="To include only the metadata for a certain variable (e.g."
                   "soil_moisture) pass the name here. This option is allowed"
                   "multiple times.")
def export_geojson(data_path, file_out, field, variable):
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
    if not os.path.exists(data_path):
        raise ValueError("The passed DATA_PATH does not exist.")
    ds = ISMN_Interface(data_path)
    if file_out is None:
        file_out = os.path.join(ds.root.root_dir, 'ismn_sensors.json')
    os.makedirs(os.path.dirname(file_out), exist_ok=True)
    print(f"Exporting geojson to: {file_out}")
    print(f"Include fields: {field}")
    print(f"Filter for variables: {variable}")

    kwargs = {}
    field = [f.lower() for f in field]
    for opt in ['network', 'station', 'sensor', 'depth', 'timerange']:
        if opt in field:
            field.pop(field.index(opt))
            kwargs[opt] = True
        else:
            kwargs[opt] = False

    kwargs['extra_props'] = field if len(field) > 0 else None

    if len(variable) > 0:
        kwargs['filter_kwargs'] = {'variable': variable}

    ds.collection.export_geojson(file_out, **kwargs)

@click.group(short_help="ISMN Command Line Programs.")
def ismn():
    pass

@click.command("nrt-download", short_help="Download the ISMN archive from ismn.earth.")
@click.argument(
    "output_path",
    type=click.Path(dir_okay=False, writable=True),
    default="/tmp/ismn_archive.zip",
    required=False,
)
@click.option(
    "--username",
    default=lambda: os.environ.get("ISMN_USERNAME"),
    help="ISMN account username. Falls back to $ISMN_USERNAME.",
)
@click.option(
    "--password",
    default=lambda: os.environ.get("ISMN_PASSWORD"),
    help="ISMN account password. Falls back to $ISMN_PASSWORD.",
)
def nrt_download(output_path, username, password):
    """Download the ISMN archive from ismn.earth.

    OUTPUT_PATH is where to write the downloaded archive
    (default: /tmp/ismn_archive.zip).
    """
    if not username or not password:
        raise click.UsageError(
            "username and password are required "
            "(pass --username/--password or set $ISMN_USERNAME/$ISMN_PASSWORD)."
        )

    ISMNDownloader(
        username=username,
        password=password,
        output_path=output_path,
    ).run()


@click.command("nrt-extract", short_help="Extract hardcoded NRT networks from an ISMN archive.")
@click.argument("archive_path", type=click.Path(exists=True, dir_okay=False))
@click.argument("output_dir", type=click.Path(file_okay=False))
def nrt_extract(archive_path, output_dir):
    """Extract hardcoded NRT networks from an ISMN archive.

    ARCHIVE_PATH is the path to the ISMN archive zip.
    OUTPUT_DIR is where to write extracted output.
    """
    ISMNExtractor(
        archive_path=archive_path,
        nrt_networks=nrt_networks,
        output_dir=output_dir,
    ).run()


ismn.add_command(collect_metadata)
ismn.add_command(export_geojson)
ismn.add_command(nrt_download, name="nrt-download")
ismn.add_command(nrt_extract, name="nrt-extract")