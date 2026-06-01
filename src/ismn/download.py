import requests
import zipfile
import shutil
from pathlib import Path
from ismn.misc import collect_stm_cov, write_overview
from ismn.const import nrt_networks
import re

class ISMNDownloader:
    BASE_URL = "https://ismn.earth"
    LOGIN_URL = f"{BASE_URL}/en/accounts/login/"
    DOWNLOAD_URL = f"{BASE_URL}/en/dataviewer/api/download_archive"

    def __init__(self,
                 username: str,
                 password: str,
                 output_path: str = "ismn_archive.zip"):
        self.username = username
        self.password = password
        self.output_path = Path(output_path).expanduser().resolve()
        self.session = requests.Session()

    def _get_csrf_token(self) -> str:
        print("Step 1: Fetching login page and CSRF token...")
        self.session.get(self.LOGIN_URL)
        csrf_token = self.session.cookies["csrftoken"]
        print(f"  ✓ CSRF token obtained: {csrf_token[:10]}...")
        return csrf_token

    def _login(self, csrf_token: str) -> None:
        print("\nStep 2: Logging in...")
        response = self.session.post(
            self.LOGIN_URL,
            data={
                "csrfmiddlewaretoken": csrf_token,
                "login": self.username,
                "password": self.password,
            },
            headers={"Referer": self.LOGIN_URL},
        )
        if response.ok:
            print(f"  ✓ Login successful (status {response.status_code})")
        else:
            raise RuntimeError(f"Login failed (status {response.status_code})")

    def _download(self) -> None:
        print("\nStep 3: Starting archive download...")
        self.output_path.parent.mkdir(parents=True, exist_ok=True)
        response = self.session.get(self.DOWNLOAD_URL, stream=True)

        if not response.ok:
            raise RuntimeError(
                f"Download request failed (status {response.status_code})")

        total_size = int(response.headers.get("content-length", 0))
        if total_size:
            print(f"  ✓ File size: {total_size / (1024**3):.2f} GB")
        else:
            print("  ! File size unknown (no Content-Length header)")

        chunk_size = 1024 * 1024  # 1 MB
        downloaded = 0

        with open(self.output_path, "wb") as f:
            for chunk in response.iter_content(chunk_size=chunk_size):
                if chunk:
                    f.write(chunk)
                    downloaded += len(chunk)
                    if total_size:
                        percent = downloaded / total_size * 100
                        downloaded_gb = downloaded / (1024**3)
                        total_gb = total_size / (1024**3)
                        print(
                            f"  Downloading... {downloaded_gb:.2f} / {total_gb:.2f} GB ({percent:.1f}%)",
                            end="\r")
                    else:
                        print(
                            f"  Downloaded: {downloaded / (1024**3):.2f} GB",
                            end="\r")

        print(
            f"\n  ✓ Download complete: {self.output_path} ({downloaded / (1024**3):.2f} GB)"
        )

    def run(self) -> None:
        csrf_token = self._get_csrf_token()
        self._login(csrf_token)
        self._download()


class ISMNExtractor:
    # Matches "<prefix>_<YYYYMMDD>.stm" -> group(1) is everything except the final date
    _STM_RE = re.compile(r"^(.*)_(\d{8})\.stm$")

    def __init__(self, archive_path: str, nrt_networks: list[str], output_dir: str = "ISMN"):
        self.archive_path = Path(archive_path)
        self.nrt_networks = nrt_networks
        self.output_dir = Path(output_dir)

    @classmethod
    def _stm_prefix(cls, name: str) -> str | None:
        """Return the filename with the trailing '_YYYYMMDD.stm' stripped, or None if it doesn't match."""
        m = cls._STM_RE.match(name)
        return m.group(1) if m else None

    def _index_existing_stm(self, network: str) -> dict[tuple[Path, str], Path]:
        """Map (parent_dir, prefix) -> existing .stm file on disk for the given network."""
        index: dict[tuple[Path, str], Path] = {}
        net_dir = self.output_dir / network
        if not net_dir.exists():
            return index

        for path in net_dir.rglob("*.stm"):
            prefix = self._stm_prefix(path.name)
            if prefix is None:
                continue
            index[(path.parent, prefix)] = path
        return index

    def run(self) -> None:
        self.output_dir.mkdir(parents=True, exist_ok=True)
        print(f"Updating {len(self.nrt_networks)} NRT networks from {self.archive_path}...")
        print(f"Output directory: {self.output_dir.resolve()}\n")

        totals = {"overwritten": 0, "no_match": 0, "missing_networks": 0}

        with zipfile.ZipFile(self.archive_path, "r") as zf:
            all_entries = zf.namelist()

            for network in self.nrt_networks:
                matching = [e for e in all_entries if e.startswith(f"{network}/")]

                if not matching:
                    print(f"  ! Network not found in archive: {network}")
                    totals["missing_networks"] += 1
                    continue

                # Build an index of the .stm files already present for this network
                existing_index = self._index_existing_stm(network)

                print(f"\n  Processing {network} ({len(matching)} entries)...")
                net_overwritten = net_no_match = 0

                for entry in matching:
                    # Skip directory entries and anything that isn't a .stm file
                    if entry.endswith("/") or not entry.endswith(".stm"):
                        continue

                    archive_name = Path(entry).name
                    prefix = self._stm_prefix(archive_name)
                    if prefix is None:
                        continue  # not a dated .stm file we care about

                    target_dir = (self.output_dir / entry).parent
                    existing = existing_index.get((target_dir, prefix))

                    if existing is None:
                        # No file with the same prefix on disk -> nothing to update
                        net_no_match += 1
                        print(f"    - no match {archive_name}")
                        continue

                    # Overwrite the existing file, keeping its (old) name
                    with zf.open(entry) as src, open(existing, "wb") as dst:
                        shutil.copyfileobj(src, dst)

                    net_overwritten += 1
                    print(f"    ~ updated  {existing.name}  (from {archive_name})")

                totals["overwritten"] += net_overwritten
                totals["no_match"] += net_no_match
                print(f"  ✓ {network}: {net_overwritten} overwritten, {net_no_match} without a match")

        print(
            f"\n✓ All done. "
            f"{totals['overwritten']} overwritten, {totals['no_match']} without a match."
        )
        if totals["missing_networks"]:
            print(f"  ! {totals['missing_networks']} requested network(s) not found in archive.")

        print("Collecting coverage information from .stm files...")
        period_to = collect_stm_cov(str(self.output_dir), n_proc=4)
        print("Writing overview file...")
        write_overview(self.output_dir, period_to=period_to)


if __name__ == "__main__":
    ISMNExtractor(
        archive_path="ismn_archive.zip",
        output_dir="/ISMN",
        nrt_networks=nrt_networks).run()
